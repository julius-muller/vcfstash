import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Dict
from .base import VEPDatabase, NextflowWorkflow
from ..utils.validation import compute_md5
from pathlib import Path
import subprocess
from datetime import datetime
from typing import Optional
from multiprocessing import Pool
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from .base import VEPDatabase, NextflowWorkflow
from ..utils.validation import validate_vcf_format, get_bcf_stats


class DatabaseAnnotator(VEPDatabase):
    """Handles database annotation workflows"""
    def __init__(self, db_path: Path, workflow_dir: Path, params_file: Optional[Path], nextflow_args: List[str]):
        super().__init__(db_path)
        self._nf_workflow = NextflowWorkflow(workflow_dir, params_file=params_file)
        # raise ValueError('here!!')
        self.workflow_dir = self._nf_workflow.workflow_dir
        self.workflow_path = self._nf_workflow.workflow_path
        self.workflow_config_path = self._nf_workflow.workflow_config_path
        self.params_file = self._nf_workflow.params_file
        self.nextflow_args = nextflow_args
        self.run_dir = self._create_run_dir(db_dir=self.db_path , workflow_hash=self._nf_workflow.workflow_hash)
        self.run_info = self.run_dir / "annotation.info"

    def annotate(self) -> None:
        """Run annotation workflow on database"""
        if not self.variants_bcf.exists():
            raise FileNotFoundError("Database BCF file does not exist.")

        # Create workflow directory in database if it doesn't exist
        db_workflow_dir = self.db_path / "workflow"
        db_workflow_dir.mkdir(exist_ok=True)

        # Copy workflow files to database
        shutil.copy2(self.workflow_path, db_workflow_dir / "main.nf")
        shutil.copy2(self.workflow_config_path, db_workflow_dir / "nextflow.config")
        if self.params_file:
            shutil.copy2(self.params_file, db_workflow_dir / self.params_file.name)

        # Store blueprint snapshot and workflow files
        shutil.copy2(self.info_file, self.run_dir / "blueprint.snapshot")
        workflow_files = self._store_workflow_files(self.workflow_dir, self.run_dir)

        try:
            start_time = datetime.now()
            result = self._nf_workflow.run(
                input_file=self.variants_bcf,
                output_dir=self.run_dir,
                db_mode=True,
                nextflow_args=self.nextflow_args
            )
            duration = datetime.now() - start_time

            # Log annotation details
            self.log_message(f"Annotated database: {self.variants_bcf}")
            self.log_message(f"Workflow hash: {self._nf_workflow.workflow_hash}")
            self.log_message(f"Annotation directory: {self.run_dir}")
            for name, hash_value in workflow_files.items():
                self.log_message(f"Workflow {name} MD5: {hash_value}")
            self.log_message(f"Processing completed in {duration.total_seconds():.2f} seconds")

            # Log tool versions if available
            tool_versions = Path(self.run_dir) / "tool_version.log"
            if tool_versions.exists():
                with open(tool_versions, 'r') as f:
                    for line in f:
                        self.log_message(f"Tool version: {line.strip()}")

            # Create archive and clean up
            archive_path = self._create_archive(self.run_dir)
            self.log_message(f"Created annotation archive: {archive_path.name}")

        except subprocess.CalledProcessError as e:
            print(f"NextFlow Error:\n{e.stderr}", file=sys.stderr)
            shutil.rmtree(self.run_dir)
            sys.exit(1)

    def _create_archive(self, run_dir: Path) -> Path:
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

    def _read_archive(self, archive_path: Path, extract_path: Optional[Path] = None) -> Dict:
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

    def _create_run_dir(self, db_dir: Path, workflow_hash: str) -> Path:
        """Create a unique directory for this annotation run using timestamp and workflow hash"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        annotations_dir = Path(db_dir) / "annotations"
        run_dir = annotations_dir / f"{timestamp}_{workflow_hash[:8]}"  # Use first 8 chars of hash
        run_dir.mkdir(parents=True, exist_ok=True)
        return run_dir

    def _store_workflow_files(self, workflow_dir: Path, target_dir: Path) -> Dict[str, str]:
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


class VCFAnnotator(VEPDatabase):
    """Handles annotation of user VCF files using the database"""

    def __init__(self, db_path: Path, input_vcf: Path, output_dir: Path, workflow_dir: Optional[Path] = None,
                 threads: int = 4):
        super().__init__(db_path)
        self.input_vcf = Path(input_vcf)
        self.output_dir = Path(output_dir)
        self.workflow_dir = workflow_dir or self.db_path / "workflow"
        self.threads = max(int(threads), 1)
        self._nf_workflow = NextflowWorkflow(self.workflow_dir)
        self.run_dir = self._create_run_dir(workflow_hash=self._nf_workflow.workflow_hash)

    @staticmethod
    def _process_region(args: tuple) -> pd.DataFrame:
        """Process a single genomic region from BCF file.

        Args:
            args: Tuple of (bcf_path, region)

        Returns:
            DataFrame with variants from the region
        """
        bcf_path, region = args
        vcf = pysam.VariantFile(bcf_path)
        records = []

        try:
            for record in vcf.fetch(region):
                # Extract basic variant information
                variant = {
                    'CHROM': record.chrom,
                    'POS': record.pos,
                    'ID': record.id or '.',
                    'REF': record.ref,
                    'ALT': ','.join(str(a) for a in record.alts) if record.alts else '.',
                    'QUAL': record.qual if record.qual is not None else '.',
                    'FILTER': ','.join(record.filter.keys()) if record.filter else 'PASS'
                }

                # Add FORMAT fields
                for sample in record.samples:
                    for field in record.format.keys():
                        value = record.samples[sample][field]
                        if isinstance(value, tuple):
                            value = ','.join(str(v) for v in value)
                        variant[f"{sample}_{field}"] = value

                # Add INFO fields
                for field, value in record.info.items():
                    if isinstance(value, tuple):
                        value = ','.join(str(v) for v in value)
                    variant[field] = value

                records.append(variant)

        except Exception as e:
            print(f"Error processing region {region}: {e}", file=sys.stderr)
            raise

        return pd.DataFrame(records)

    @staticmethod
    def validate_vcf_format(vcf_path: Path) -> tuple[bool, str] | tuple[bool, None]:
        """Validate VCF format fields.

        Args:
            vcf_path: Path to the VCF file

        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            vcf = pysam.VariantFile(str(vcf_path))
            required_formats = {'AD', 'DP', 'GT'}
            available_formats = set(vcf.header.formats.keys())

            missing_formats = required_formats - available_formats
            if missing_formats:
                return False, f"Missing required FORMAT fields: {', '.join(missing_formats)}"

            return True, None
        except Exception as e:
            return False, f"Error reading VCF file: {e}"

    def annotate(self, convert_parquet: bool = True) -> Path:
        """Run annotation workflow on input VCF file.

        Args:
            convert_parquet: Whether to convert output to Parquet format

        Returns:
            Path to output file (BCF or Parquet)
        """
        start_time = time.time()

        try:
            self._validate_inputs()
            annotated_bcf = self._run_workflow()
            vep_time = time.time()
            print(f"[{datetime.now()}] VEP workflow completed in {vep_time - start_time:.1f}s")

            if not convert_parquet:
                return annotated_bcf

            print(f"[{datetime.now()}] Converting to Parquet format...")
            output_file = self._convert_to_parquet(annotated_bcf)
            end_time = time.time()

            print(f"[{datetime.now()}] Parquet conversion completed in {end_time - vep_time:.1f}s")
            print(f"[{datetime.now()}] Total execution time: {end_time - start_time:.1f}s")

            return output_file

        except Exception as e:
            print(f"Error during annotation: {e}", file=sys.stderr)
            if self.run_dir.exists():
                shutil.rmtree(self.run_dir)
            raise


    def _validate_inputs(self) -> None:
        """Validate input VCF format and requirements"""
        if not self.input_vcf.exists():
            raise FileNotFoundError(f"Input VCF file not found: {self.input_vcf}")

        is_valid, error = validate_vcf_format(self.input_vcf)
        if not is_valid:
            raise ValueError(f"Invalid VCF file: {error}")

        self.ensure_indexed(self.input_vcf)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def _run_workflow(self) -> Path:
        """Run the Nextflow annotation workflow"""
        sample_name = self.input_vcf.stem.split('.')[0]
        final_bcf = self.run_dir / f"{sample_name}_norm_final.bcf"

        cmd = [
            "nextflow", "run", str(self.workflow_dir / "main.nf"),
            "--input", str(self.input_vcf),
            "--output", str(self.run_dir),
            "--db_bcf", str(self.variants_bcf),
            "--db_mode", "false",
            "--vep_max_forks", str(min(self.threads, 4)),
            "--vep_max_chr_parallel", str(min(self.threads, 4)),
            "-with-trace"
        ]

        try:
            start_time = datetime.now()
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            duration = datetime.now() - start_time

            self._nf_workflow.store_workflow_dag(self.run_dir, cmd)
            stats = get_bcf_stats(final_bcf)

            # Log annotation details
            self.log_message(f"Annotated VCF file: {self.input_vcf}")
            self.log_message(f"Using database: {self.variants_bcf}")
            self.log_message(f"Command executed: {' '.join(cmd)}")
            self.log_message(f"Processing completed in {duration.total_seconds():.2f} seconds")
            self.log_message("Output statistics:")
            for key, value in stats.items():
                self.log_message(f"{key}: {value}")

            return final_bcf

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Annotation workflow failed: {e}")

    def _convert_to_parquet(self, bcf_path: Path) -> Path:
        """Convert annotated BCF to optimized Parquet format"""
        vcf = pysam.VariantFile(str(bcf_path))
        regions = list(vcf.header.contigs.keys())
        args_list = [(str(bcf_path), region) for region in regions]

        with Pool(self.threads) as pool:
            dataframes = pool.map(self._process_region, args_list)

        # Filter and combine dataframes
        dataframes = [df for df in dataframes if not df.empty]
        if not dataframes:
            raise ValueError("No valid variants found in annotated file")

        combined_df = pd.concat(dataframes, ignore_index=True)
        output_file = self.output_dir / f"{self.input_vcf.stem}.parquet"

        # Write optimized parquet
        table = pa.Table.from_pandas(combined_df)
        pq.write_table(
            table,
            output_file,
            compression="snappy",
            use_dictionary=True,
            row_group_size=100000,
            data_page_size=65536,
            write_statistics=True
        )

        return output_file

    def _create_run_dir(self, workflow_hash: str) -> Path:
        """Create a unique directory for this annotation run.

        Args:
            workflow_hash: Hash of the workflow files

        Returns:
            Path to the created run directory
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        annotations_dir = self.output_dir / "annotations"
        run_dir = annotations_dir / f"{timestamp}_{workflow_hash[:8]}"  # Use first 8 chars of hash
        run_dir.mkdir(parents=True, exist_ok=True)
        return run_dir
