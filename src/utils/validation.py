import re
import sys
from typing import Tuple, Optional, Dict
from pathlib import Path
import subprocess
import pysam


def check_duplicate_md5(db_info: dict, new_md5: str) -> bool:
    """Check if a file with the same MD5 was already added"""
    try:
        return any(f["md5"] == new_md5 for f in db_info.get("input_files", []))
    except KeyError:
        return False


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



def check_bcftools_installed() -> None:
    try:
        subprocess.run(["bcftools", "--version"], check=True, capture_output=True)
    except FileNotFoundError:
        sys.exit("Error: bcftools is not installed or not in PATH.")




def compute_md5(file_path: Path) -> str:
    """

    :param file_path:
    :return:
    file_path = Path('~/projects/vcfstash/tests/data/nodata/dbsnp_test.bcf')
    """
    try:
        print(f'Computing MD5 for {file_path} ...')
        result = subprocess.run(
            ["md5sum", file_path],
            check=True,
            capture_output=True,
            text=True
        )
        return result.stdout.split()[0]  # The MD5 hash is the first word in the output
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error computing MD5 checksum: {e}")

def validate_vcf_format(vcf_path: Path) -> tuple[bool, str | None]:
    """Validate VCF format fields.

    Args:
        vcf_path: Path to the VCF file

    Returns:
        Tuple of (is_valid, error_message)
    """
    try:
        vcf = pysam.VariantFile(str(vcf_path))

        # Ensure the file can be read
        try:
            next(vcf.fetch())
        except StopIteration:
            return False, "VCF file is empty"
        except Exception as e:
            return False, f"Error reading VCF file: {e}"

        # Check for minimal required FORMAT fields
        required_formats = {'GT', 'AD'}  # Removed DP requirement
        available_formats = set(vcf.header.formats.keys())

        missing_formats = required_formats - available_formats
        if missing_formats:
            return False, f"Missing required FORMAT fields: {', '.join(missing_formats)}"

        return True, None

    except Exception as e:
        return False, f"Error reading VCF file: {e}"


def generate_test_command(
        vcfstash_path="~/projects/vcfstash/vcfstash.py",
        vcf_path="~/projects/vcfstash/tests/data/nodata/crayz_db.bcf",
        output_dir="~/tmp/vcfstash/test_stash",
        config_path="~/projects/vcfstash/tests/config/env_test.config",
        yaml_path="~/projects/vcfstash/tests/config/user_params.yaml",
        annotation_config = "~/projects/vcfstash/tests/config/annotation.config",
        add_vcf_path="~/projects/vcfstash/tests/data/nodata/crayz_db2.bcf",
        input_vcf_path="~/projects/vcfstash/tests/data/nodata/sample4.bcf",
        annotate_name="testor",
        annotation_db="~/tmp/vcfstash/test_stash/stash/testor",
        annotation_output="~/tmp/vcfstash/aout",
        force=True
):
    """
    Generate a nicely formatted test command string for vcfstash operations.

    Returns:
        str: A copy-pastable command string with proper formatting
    """

    cmd_init = (
        f"{vcfstash_path} stash-init "
        f"--vcf {vcf_path} "
        f"--output {output_dir} "
        f"-y {yaml_path} "
        f"{'-f' if force else ''} "
    ).strip()

    cmd_add = (
        f"{vcfstash_path} stash-add "
        f"--db {output_dir} "
        f"-i {add_vcf_path} "
    ).strip()

    cmd_annotate = (
        f"{vcfstash_path} stash-annotate "
        f"--name {annotate_name} "
        f"-a {annotation_config} "
        f"--db {output_dir} "
        f"{'-f' if force else ''} "
    ).strip()

    cmd_annotatevcf = (
        f"{vcfstash_path} annotate "
        f"-a {annotation_db} "
        f"--vcf {input_vcf_path} "
        f"{'-f' if force else ''} "
        f"--output {annotation_output} "
    ).strip()

    # Combine commands
    full_cmd = f"{cmd_init} ; {cmd_add} ; {cmd_annotate} ; {cmd_annotatevcf}"

    # Also create a nicely formatted display version for easier reading
    formatted_cmds = f"""
# INITIALIZE
{cmd_init}

# ADD
{cmd_add}

# STASH ANNOTATE
{cmd_annotate}

# ANNOTATE VCF
{cmd_annotatevcf}

# COMBINED COMMAND (for copy-paste)
alias stx="{full_cmd}"
"""

    print(formatted_cmds)
    return full_cmd

# generate_test_command()


# Example usage into test dir in repo:
# ~/projects/vcfstash/vcfstash.py stash-init --name nftest --vcf ~/projects/vcfstash/tests/data/nodata/crayz_db.bcf --output /home/j380r/tmp/test/test_out -f --test -vv
# ~/projects/vcfstash/vcfstash.py stash-add --db ~/projects/vcfstash/tests/data/test_out/nftest -i ~/projects/vcfstash/tests/data/nodata/crayz_db2.bcf --test -vv
# ~/projects/vcfstash/vcfstash.py stash-annotate --name testor --db ~/projects/vcfstash/tests/data/test_out/nftest --test -vv -f
# ... or locally
# ~/projects/vcfstash/vcfstash.py stash-init --vcf ~/projects/vcfstash/tests/data/nodata/crayz_db.bcf --output ~/tmp/vcfstash/test_stash -c ~/projects/vcfstash/tests/config/env_test.config -f
# ~/projects/vcfstash/vcfstash.py stash-add --db /home/j380r/tmp/test/test_out -i ~/projects/vcfstash/tests/data/nodata/crayz_db2.bcf
# ~/projects/vcfstash/vcfstash.py stash-annotate --name testor --db test_out/nftest --test -vv -f
# ~/projects/vcfstash/vcfstash.py annotate --a ~/tmp/test/test_out/nftest/stash/testor --vcf ~/projects/vcfstash/tests/data/nodata/sample4.bcf --output ~/tmp/test/aout --test -f

# as one:
cmd = """alias stx="
~/projects/vcfstash/vcfstash.py stash-init --vcf ~/projects/vcfstash/tests/data/nodata/crayz_db.bcf --output ~/tmp/vcfstash/test_stash -y ~/projects/vcfstash/tests/config/user_params.yaml -f;
~/projects/vcfstash/vcfstash.py stash-add --db ~/tmp/vcfstash/test_stash/ -i ~/projects/vcfstash/tests/data/nodata/crayz_db2.bcf ; 
~/projects/vcfstash/vcfstash.py stash-annotate --name testor -a ~/projects/vcfstash/tests/config/annotation.config --db ~/tmp/vcfstash/test_stash -f;
~/projects/vcfstash/vcfstash.py annotate -a ~/tmp/vcfstash/test_stash/stash/testor --vcf ~/projects/vcfstash/tests/data/nodata/sample4.bcf --output ~/tmp/vcfstash/aout -f
"""



# on gvpre
cmd2 = """
~/projects/vcfstash/vcfstash.py stash-init --vcf /mnt/data/resources/gnomad/vcf_gnomad_v4_hg19_exomes/gnomad.exomes.v4.1.sites.grch37.trimmed_liftover_norm_1e-1.bcf --output gnomad_1e-1  -c ~/projects/vcfstash/tests/config/nextflow_gnomadhg19.config;
~/projects/vcfstash/vcfstash.py stash-add --db gnomad_1e-1 -i /mnt/data/resources/gnomad/vcf_gnomad_v4_hg19_genomes/gnomad.genomes.v4.1.sites.grch37.trimmed_liftover_norm_1e-1.bcf;
~/projects/vcfstash/vcfstash.py stash-annotate --name gen_ex -a ~/projects/vcfstash/tests/config/annotation.config --db gnomad_1e-1;
~/projects/vcfstash/vcfstash.py annotate -a gnomad_1e-1/stash/gen_ex --vcf /mnt/data/samples/test_mgm/mgm_WGS_32.gatkWGS_norm.bcf --output mgm_WGS_32 -p;
"""
