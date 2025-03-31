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
    file_path = Path('~/projects/vepstash/tests/data/nodata/dbsnp_test.bcf')
    """
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



def get_vep_version_from_cmd(vep_cmd):
    """
    Executes the VEP command and extracts the version number from its output.

    Args:
        vep_cmd (str): The VEP command to execute.

    Returns:
        str: The VEP version number (e.g., "113.0").

    Raises:
        ValueError: If version cannot be extracted from output or command fails.

        vep_cmd ="docker run --user \\$(id -u):\\$(id -g) -i -v ${vep_cache}:${vep_cache} -v ${refdir}:${refdir} --rm ensemblorg/ensembl-vep:release_113.0 vep"
    """
    try:
        shell_cmd = vep_cmd.replace("\\$(", "$(").replace("\\${", "${")

        # Add --help flag to trigger version display without requiring inputs
        cmd = f"{shell_cmd} --help"
        # print(cmd)
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            error_msg = f"VEP command failed with exit code {result.returncode}: {result.stderr.strip()}"
            raise ValueError(error_msg)

        # Parse the output to find the version
        output = result.stdout

        # Look for the ensembl-vep version line
        version_match = re.search(r"ensembl-vep\s*:\s*(\d+\.\d+(?:\.\d+)?)", output)
        if not version_match:
            # Alternative pattern for older versions or different formats
            version_match = re.search(r"VEP version (\d+\.\d+(?:\.\d+)?)", output)

        if version_match:
            version = version_match.group(1)
            return version
        else:
            error_msg = "Could not determine VEP version from command output."
            raise ValueError(error_msg)

    except Exception as e:
        error_msg = f"Error determining VEP version: {str(e)}"
        raise ValueError(error_msg)


def get_echtvar_version_from_cmd(echtvar_cmd):
    """
    Executes the echtvar command and extracts the version number from its output.

    Args:
        echtvar_cmd (str): The echtvar command to execute.

    Returns:
        str: The echtvar version number (e.g., "0.2.1").

    Raises:
        ValueError: If version cannot be extracted from output or command fails.
    """
    try:
        # Use -V flag to get version information
        cmd = f"{echtvar_cmd} -V"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            # Try alternative --version flag if -V fails
            cmd = f"{echtvar_cmd} --version"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.returncode != 0:
                error_msg = f"Echtvar command failed with exit code {result.returncode}: {result.stderr.strip()}"
                raise ValueError(error_msg)

        # Parse the output to find the version
        output = result.stdout

        # Try pattern for "echtvar 0.2.1" format
        version_match = re.search(r"echtvar\s+(\d+\.\d+\.\d+)", output)
        if not version_match:
            # Alternative pattern just looking for version number
            version_match = re.search(r"(\d+\.\d+\.\d+)", output)

        if version_match:
            version = version_match.group(1)
            return version
        else:
            error_msg = "Could not determine echtvar version from command output."
            raise ValueError(error_msg)

    except Exception as e:
        error_msg = f"Error determining echtvar version: {str(e)}"
        raise ValueError(error_msg)


