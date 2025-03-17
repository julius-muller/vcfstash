
"""
VEP Annotation Cache

This script manages a database of genetic variants in BCF/VCF format,
providing functionality to initialize, add to, and annotate variant databases.

Key features:
- Supports BCF/VCF format (uncompressed/compressed)
- Requires pre-indexed input files (CSI/TBI index)
  -> This ensures input is properly sorted and valid
- Maintains database integrity through MD5 checksums
- Provides versioned annotation workflow support
- Includes detailed logging of all operations

Author: Julius Müller, PhD
Organization: GHGA - German Human Genome-Phenome Archive
Date: 16-03-2025

"""

import tarfile
import tempfile
from datetime import datetime
import json
import os
import sys
import subprocess
import argparse
from pathlib import Path
from datetime import datetime
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional


def create_annotation_archive(run_dir: Path) -> Path:
    """Create a single archive file from annotation run directory"""
    archive_name = f"{run_dir.name}.vepstash"
    archive_path = run_dir.parent / archive_name

    # Create archive with metadata
    with tarfile.open(archive_path, "w:gz") as tar:
        # Add all files from run directory
        tar.add(run_dir, arcname="")

        # Add metadata file
        metadata = {
            "created": datetime.now().isoformat(),
            "format_version": "1.0",
            "run_id": run_dir.name
        }
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tf:
            json.dump(metadata, tf, indent=2)
            tf_path = Path(tf.name)

        tar.add(tf_path, arcname="metadata.json")
        tf_path.unlink()

    # Remove original directory after successful archive creation
    shutil.rmtree(run_dir)
    return archive_path

def read_annotation_archive(archive_path: Path, extract_path: Optional[Path] = None) -> Dict:
    """Read metadata from annotation archive without extracting"""
    with tarfile.open(archive_path, "r:gz") as tar:
        try:
            meta_file = tar.extractfile("metadata.json")
            if meta_file:
                metadata = json.loads(meta_file.read().decode())
                if extract_path:
                    tar.extractall(path=extract_path)
                return metadata
        except KeyError:
            raise ValueError("Invalid annotation archive: metadata.json not found")
    return {}


def log_message(log_file: Path, message: str, level: str = "INFO") -> None:
    """Log a message with timestamp and level"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    formatted_message = f"[{timestamp}] {level}: {message}"
    with open(log_file, "a") as log:
        log.write(formatted_message + "\n")

def check_duplicate_md5(info_file: Path, new_md5: str) -> bool:
    """Check if a file with the same MD5 was already added"""
    try:
        with open(info_file, "r") as f:
            for line in f:
                if "Input file MD5:" in line and new_md5 in line:
                    return True
    except FileNotFoundError:
        return False
    return False

def ensure_indexed(file_path: Path) -> None:
    """Ensure input file has an index file (CSI or TBI)"""
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    csi_index = Path(str(file_path) + ".csi")
    tbi_index = Path(str(file_path) + ".tbi")

    if not (csi_index.exists() or tbi_index.exists()):
        raise RuntimeError(
            f"No index found for {file_path}. Use bcftools index for BCF/compressed VCF "
            "or tabix for compressed VCF files."
        )

def validate_bcf_header(bcf_path: Path, norm: bool = True) -> Tuple[bool, Optional[str]]:
    """
    Validate BCF header for required normalization command and contig format.
    Returns tuple (is_valid, error_message).
    """
    try:
        header = subprocess.run(
            ["bcftools", "view", "-h", bcf_path],
            check=True,
            capture_output=True,
            text=True
        ).stdout

        if norm:
            # Check normalization command
            norm_lines = [line for line in header.splitlines()
                         if line.startswith("##bcftools_normCommand")]

            if not norm_lines:
                return False, "Missing bcftools_normCommand in header"

            norm_cmd = norm_lines[0]
            required_options = ["norm", "-c x", "-m-"]
            missing_options = [opt for opt in required_options if opt not in norm_cmd]

            if missing_options:
                return False, f"Missing required normalization options: {', '.join(missing_options)}"

        # Check contig format
        contig_lines = [line for line in header.splitlines()
                       if line.startswith("##contig=")]

        if not contig_lines:
            return False, "No contig lines found in header"

        invalid_contigs = [line for line in contig_lines
                          if not line.startswith("##contig=<ID=chr")]

        if invalid_contigs:
            example = invalid_contigs[0] if invalid_contigs else ""
            return False, f"Invalid contig format (should start with 'chr'): {example}"

        return True, None

    except subprocess.CalledProcessError as e:
        return False, f"Error reading BCF header: {e}"

def get_bcf_stats(bcf_path: Path) -> Dict[str, str]:
    """Get statistics from BCF file using bcftools stats"""
    try:
        result = subprocess.run(
            ["bcftools", "stats", bcf_path],
            capture_output=True,
            text=True,
            check=True
        )
        stats = {}
        for line in result.stdout.splitlines():
            if line.startswith("SN"):
                parts = line.split("\t")
                if len(parts) >= 4:
                    key = parts[2].strip(":")
                    value = parts[3]
                    stats[key] = value
        return stats
    except subprocess.CalledProcessError as e:
        return {"error": f"Failed to get statistics: {e}"}


def compute_md5(file_path: Path) -> str:
    try:
        result = subprocess.run(
            ["md5sum", file_path],
            check=True,
            capture_output=True,
            text=True
        )
        return result.stdout.split()[0]  # The MD5 hash is the first word in the output
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error computing MD5 checksum: {e}")


def check_bcftools_installed() -> None:
    try:
        subprocess.run(["bcftools", "--version"], check=True, capture_output=True)
    except FileNotFoundError:
        sys.exit("Error: bcftools is not installed or not in PATH.")


def log_script_command(info_file: Path) -> None:
    """Log the exact command used to execute the script"""
    command = f"python3 {sys.argv[0]} {' '.join(sys.argv[1:])}"
    log_message(info_file, f"Script command: {command}")

def store_workflow_dag(run_dir: Path, cmd: List[str]) -> None:
    """Generate and store workflow DAG visualization"""

    try:
        del cmd[cmd.index("-with-trace")]
    except ValueError:
        pass  # Element not in list

    try:
        # Generate DAG in dot format
        dag_cmd = cmd + [ "-preview", "-with-dag", str(run_dir / "flowchart.html") ]

        subprocess.run(dag_cmd, check=True, capture_output=True)


    except subprocess.CalledProcessError as e:
        log_message(run_dir / "annotation.info", f"Warning: Failed to generate workflow DAG: {e}", level="WARN")


def init_mode(args: argparse.Namespace) -> None:
    ensure_indexed(args.i)
    is_valid, error = validate_bcf_header(args.i, norm=False)
    if not is_valid:
        sys.exit(f"Invalid input BCF: {error}")

    output_dir = Path(args.output) / args.name
    output_dir.mkdir(parents=True, exist_ok=True)
    output_bcf = output_dir / "blueprint" / "variants.bcf"
    info_file = output_dir / "blueprint" / "variants.info"

    if not Path(args.i).exists():
        sys.exit("Error: Input BCF file does not exist.")
    if output_bcf.exists():
        sys.exit("Error: Output database already exists.")

    threads = str(max(int(args.t), 1))
    input_md5 = compute_md5(args.i)

    try:
        cmd = [
            "bcftools", "view", "-Ou", "--threads", threads, args.i,
            "|", "bcftools", "annotate", "-x", "INFO", "--threads", threads,
            "|", "bcftools", "norm", "-m-", "-f", args.fasta, "-c", "x", "--threads", threads, "--rm-dup", "all",
            "-Ob", "--write-index", "-o", str(output_bcf)
        ]

        start_time = datetime.now()
        subprocess.run(" ".join(cmd), shell=True, check=True)
        duration = datetime.now() - start_time

        info_file.write_text("")  # Create empty info file
        log_script_command(info_file)
        log_message(info_file, f"Initialized VEP database: {args.name}")
        log_message(info_file, f"Input file: {args.i}")
        log_message(info_file, f"Input file MD5: {input_md5}")
        log_message(info_file, f"Command executed: {' '.join(cmd)}")
        log_message(info_file, f"Processing completed in {duration.total_seconds():.2f} seconds")
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error during bcftools operation: {e}")

def replace_db(db_bcf: Path, temp_output_bcf: Path) -> None:
    # Add validation for both files
    is_valid, error = validate_bcf_header(db_bcf)
    if not is_valid:
        sys.exit(f"Invalid input BCF: {error}")
    is_valid, error = validate_bcf_header(temp_output_bcf)
    if not is_valid:
        sys.exit(f"Invalid database BCF: {error}")

    # Get statistics before replacement
    info_file = Path(db_bcf).parent / "vep_db.bcf.info"
    stats = get_bcf_stats(temp_output_bcf)
    new_md5 = compute_md5(temp_output_bcf)

    # Replace files
    os.replace(temp_output_bcf, db_bcf)
    os.replace(str(temp_output_bcf) + '.csi', str(db_bcf) + '.csi')

    # Log statistics
    log_message(info_file, "Database statistics after update:")
    log_message(info_file, f"Database MD5: {new_md5}")
    if "number of SNPs:" in stats:
        log_message(info_file, f"Number of SNPs: {stats['number of SNPs']}")
    if "number of indels:" in stats:
        log_message(info_file, f"Number of indels: {stats['number of indels']}")
    if "number of records:" in stats:
        log_message(info_file, f"Total variants: {stats['number of records']}")
    if "number of samples:" in stats:
        log_message(info_file, f"Number of samples: {stats['number of samples']}")


def create_unique_annotation_dir(db_dir: Path, workflow_hash: str) -> Path:
    """Create a unique directory for this annotation run using timestamp and workflow hash"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    annotations_dir = Path(db_dir) / "annotations"
    run_dir = annotations_dir / f"{timestamp}_{workflow_hash[:8]}"  # Use first 8 chars of hash
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def get_workflow_hash(workflow_dir: Path) -> str:
    """Get combined hash of workflow files"""
    workflow_files = ['main.nf', 'nextflow.config']
    combined_content = b''

    for file in workflow_files:
        path = Path(workflow_dir) / file
        if not path.exists():
            raise FileNotFoundError(f"Required workflow file not found: {path}")
        combined_content += path.read_bytes()

    import hashlib
    return hashlib.md5(combined_content).hexdigest()

def store_workflow_files(workflow_dir: Path, target_dir: Path) -> Dict[str, str]:
    """Store workflow files with their hash for version control"""
    workflow_files = {
        'main.nf': workflow_dir / 'main.nf',
        'nextflow.config': workflow_dir / 'nextflow.config'
    }

    stored_files = {}
    for name, path in workflow_files.items():
        if not path.exists():
            raise FileNotFoundError(f"Required workflow file not found: {path}")

        # Copy file
        target_path = target_dir / name
        with open(path, 'rb') as src, open(target_path, 'wb') as dst:
            dst.write(src.read())

        # Store hash
        stored_files[name] = compute_md5(target_path)

    return stored_files

def add_mode(args: argparse.Namespace) -> None:
    db_bcf = Path(args.db) / "blueprint" / "variants.bcf"
    info_file = Path(args.db) / "blueprint" / "variants.info"
    new_vcf = Path(args.vcf)

    ensure_indexed(db_bcf)
    ensure_indexed(new_vcf)

    threads = str(max(int(args.t), 1))

    input_md5 = compute_md5(new_vcf)
    if check_duplicate_md5(info_file, input_md5):
        sys.exit("Error: This file was already added to the database (MD5 match).")

    temp_merged = Path(str(db_bcf) + ".tmp.bcf")
    try:
        start_time = datetime.now()

        subprocess.run([
            "bcftools", "concat",
            "--allow-overlaps",
            "--rm-dup", "all",
            "-Ob",
            "--write-index",
            "-o", str(temp_merged), "--threads", threads,
            str(db_bcf), str(new_vcf)
        ], check=True)

        duration = datetime.now() - start_time

        stats = get_bcf_stats(temp_merged)

        os.replace(temp_merged, db_bcf)
        os.replace(str(temp_merged) + ".csi", str(db_bcf) + ".csi")

        log_script_command(info_file)
        log_message(info_file, f"Added new file: {new_vcf}")
        log_message(info_file, f"Input file MD5: {input_md5}")
        log_message(info_file, "Database statistics after update:")
        for key, value in stats.items():
            log_message(info_file, f"{key}: {value}")
        log_message(info_file, f"Processing completed in {duration.total_seconds():.2f} seconds")

    except subprocess.CalledProcessError as e:
        temp_merged.unlink(missing_ok=True)
        Path(str(temp_merged) + ".csi").unlink(missing_ok=True)
        raise RuntimeError(f"Failed to add variants: {e}")

def annotate_mode(args: argparse.Namespace) -> None:
    """Run annotation workflow on database"""
    db_bcf = Path(args.db) / "blueprint" / "variants.bcf"
    db_info = Path(args.db) / "blueprint" / "variants.info"

    if not db_bcf.exists():
        sys.exit("Error: Database BCF file does not exist.")

    # Create unique run directory
    workflow_hash = get_workflow_hash(args.workflow)
    run_dir = create_unique_annotation_dir(args.db, workflow_hash)
    run_info = run_dir / "annotation.info"

    # Store blueprint snapshot and workflow files
    shutil.copy2(db_info, run_dir / "blueprint.snapshot")
    workflow_files = store_workflow_files(Path(args.workflow), run_dir)

    try:
        # Run nextflow workflow
        cmd = [
            "nextflow", "run", str(run_dir / "main.nf"),
            "-with-trace",
            "--input", str(db_bcf),
            "--output", str(run_dir),
            "--db_mode", "true"
        ]
        if args.params:
            cmd.extend(["--params-file", args.params])

        # Add any additional Nextflow arguments
        if args.nextflow_args:
            # Remove the '--' separator if present
            if args.nextflow_args[0] == '--':
                cmd.extend(args.nextflow_args[1:])
            else:
                cmd.extend(args.nextflow_args)

        start_time = datetime.now()
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)

        # Store workflow DAG and log completion
        store_workflow_dag(run_dir, cmd)
        duration = datetime.now() - start_time

        # Log annotation details
        log_message(run_info, f"Annotated database: {db_bcf}")
        log_message(run_info, f"Workflow hash: {workflow_hash}")
        log_message(run_info, f"Annotation directory: {run_dir}")
        log_message(run_info, f"Command executed: {' '.join(cmd)}")
        for name, hash_value in workflow_files.items():
            log_message(run_info, f"Workflow {name} MD5: {hash_value}")
        log_message(run_info, f"Processing completed in {duration.total_seconds():.2f} seconds")

        # Log tool versions if available
        tool_versions = Path(run_dir) / "tool_version.log"
        if tool_versions.exists():
            with open(tool_versions, 'r') as f:
                for line in f:
                    log_message(run_info, f"Tool version: {line.strip()}")

        # Create archive and clean up
        archive_path = create_annotation_archive(run_dir)
        log_message(db_info, f"Created annotation archive: {archive_path.name}")

    except subprocess.CalledProcessError as e:
        shutil.rmtree(run_dir)
        sys.exit(f"Error during workflow execution: {e}")

def main() -> None:
    check_bcftools_installed()

    parser = argparse.ArgumentParser(
        description="Manage VEP database with BCFtools."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # init command
    init_parser = subparsers.add_parser("init", help="Initialize VEP database")
    init_parser.add_argument("-n", "--name", dest="name", required=True, help="Name of the new database directory")
    init_parser.add_argument("-i", "--vcf", dest="i", help="CSI-indexed BCF file")
    init_parser.add_argument("-t", "--threads", dest="t", help="Use multithreading with <int> worker threads [4]", default=4)
    init_parser.add_argument("-f", "--fasta", dest="fasta", required=True, help="Reference FASTA file for variant normalization")
    init_parser.add_argument("-o", "--output", dest="output", default=".", help="Output directory for the database (default: current directory)")

    # add command
    add_parser = subparsers.add_parser("add", help="Add new VCF to the database")
    add_parser.add_argument("-d", "--db", required=True, help="Path to the existing database directory")
    add_parser.add_argument("-i", "--vcf", help="Path to the VCF file to be added")
    add_parser.add_argument("-t", "--threads", dest="t", help="Use multithreading with <int> worker threads [4]", default=4)
    add_parser.add_argument("-f", "--fasta", dest="fasta", required=True, help="Reference FASTA file for variant normalization")

    # annotate command

    # Annotate command
    annotate_parser = subparsers.add_parser("annotate", help="Run annotation workflow on database")
    annotate_parser.add_argument("-d", "--db", required=True, help="Path to the existing database directory")
    annotate_parser.add_argument("-w", "--workflow", required=True, help="Directory containing workflow files (main.nf, nextflow.config)")
    annotate_parser.add_argument("-p", "--params", help="Optional parameters file for nextflow workflow")
    # Add support for passing through Nextflow args
    annotate_parser.add_argument('nextflow_args', nargs=argparse.REMAINDER, help='Additional arguments for Nextflow')

    args = parser.parse_args()

    if args.command == "init":
        # Create info file first so we can log the command
        output_dir = Path(args.name)
        output_dir.mkdir(parents=True, exist_ok=True)
        info_file = output_dir / "vep_db.bcf.info"
        info_file.write_text("")
        log_script_command(info_file)
        init_mode(args)
    elif args.command == "add":
        # Log command to existing info file
        info_file = Path(args.db) / "vep_db.bcf.info"
        log_script_command(info_file)
        add_mode(args)
    elif args.command == "annotate":
        info_file = Path(args.db) / "vep_db.bcf.info"
        log_script_command(info_file)
        annotate_mode(args)

if __name__ == "__main__":
    main()
