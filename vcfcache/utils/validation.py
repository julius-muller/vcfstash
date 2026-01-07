# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""Validation utilities for the vcfcache package.

This module provides functions for validating VCF/BCF files, checking dependencies,
computing MD5 checksums, and other validation-related tasks.
"""

import pysam
import logging
import os
import subprocess
import sys
import shutil
import re
from pathlib import Path
from typing import Dict, Optional, Tuple
import yaml

# Minimum required bcftools version
MIN_BCFTOOLS_VERSION = "1.20"

def check_duplicate_md5(db_info: dict, new_md5: str) -> bool:
    """Check if a file with the same MD5 was already added."""
    try:
        return any(f["md5"] == new_md5 for f in db_info.get("input_files", []))
    except KeyError:
        return False


def get_bcf_stats(bcf_path: Path, bcftools_path: Path = None) -> Dict[str, str]:
    """Get statistics from BCF file using bcftools stats.

    Args:
        bcf_path: Path to the BCF file
        bcftools_path: Path to the bcftools binary (required)
    """
    if bcftools_path is None:
        raise ValueError("bcftools_path must be provided. A specific bcftools path is required.")

    try:
        result = subprocess.run(
            [str(bcftools_path), "stats", bcf_path], capture_output=True, text=True, check=True
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


def validate_bcf_header(
    bcf_path: Path, norm: bool = True, bcftools_path: Path = None
) -> Tuple[bool, Optional[str]]:
    """Validate BCF header for required normalization command and contig format.

    Args:
        bcf_path: Path to the BCF file
        norm: Whether to check for normalization command
        bcftools_path: Path to the bcftools binary (required)

    Returns:
        tuple: (is_valid, error_message)
    """
    if bcftools_path is None:
        raise ValueError("bcftools_path must be provided. A specific bcftools path is required.")

    try:
        header = subprocess.run(
            [str(bcftools_path), "view", "-h", bcf_path],
            check=True,
            capture_output=True,
            text=True,
        ).stdout

        if norm:
            # Check normalization command
            norm_lines = [
                line
                for line in header.splitlines()
                if line.startswith("##bcftools_normCommand")
            ]

            if not norm_lines:
                return False, "Missing bcftools_normCommand in header"

            norm_cmd = norm_lines[0]
            required_options = ["norm", "-c x", "-m-"]
            missing_options = [opt for opt in required_options if opt not in norm_cmd]

            if missing_options:
                return (
                    False,
                    f"Missing required normalization options: {', '.join(missing_options)}",
                )

        # Check contig presence (allow both with or without chr prefix)
        contig_lines = [
            line for line in header.splitlines() if line.startswith("##contig=")
        ]

        if not contig_lines:
            return False, "No contig lines found in header"

        return True, None

    except subprocess.CalledProcessError as e:
        return False, f"Error reading BCF header: {e}"


def parse_bcftools_version(version_string: str) -> tuple[int, int, int]:
    """Parse bcftools version string into tuple of integers.

    Args:
        version_string: Version string like "1.20" or "1.22+htslib-1.22"

    Returns:
        Tuple of (major, minor, patch) integers

    Examples:
        >>> parse_bcftools_version("1.20")
        (1, 20, 0)
        >>> parse_bcftools_version("1.22+htslib-1.22")
        (1, 22, 0)
    """
    # Extract version numbers from string like "1.20" or "1.22+htslib-1.22"
    match = re.match(r'(\d+)\.(\d+)(?:\.(\d+))?', version_string)
    if not match:
        raise ValueError(f"Could not parse bcftools version: {version_string}")

    major = int(match.group(1))
    minor = int(match.group(2))
    patch = int(match.group(3)) if match.group(3) else 0

    return (major, minor, patch)


def compare_versions(version1: str, version2: str) -> int:
    """Compare two version strings.

    Args:
        version1: First version string
        version2: Second version string

    Returns:
        -1 if version1 < version2
        0 if version1 == version2
        1 if version1 > version2
    """
    v1 = parse_bcftools_version(version1)
    v2 = parse_bcftools_version(version2)

    if v1 < v2:
        return -1
    elif v1 > v2:
        return 1
    else:
        return 0


def find_bcftools() -> Optional[str]:
    """Find bcftools binary, checking VCFCACHE_BCFTOOLS env var first, then system PATH.

    Environment variable VCFCACHE_BCFTOOLS can be used to override the system bcftools.
    This is useful when the system bcftools is too old or not available.

    Returns:
        Path to bcftools binary if found, None otherwise

    Example:
        >>> # Override system bcftools
        >>> os.environ['VCFCACHE_BCFTOOLS'] = '/opt/bcftools-1.22/bin/bcftools'
        >>> bcftools_path = find_bcftools()
    """
    # Check environment variable first
    env_path = os.environ.get("VCFCACHE_BCFTOOLS")
    if env_path:
        # Expand user paths and resolve
        env_path_resolved = Path(env_path).expanduser().resolve()
        # Verify it exists and is executable
        if env_path_resolved.exists() and os.access(env_path_resolved, os.X_OK):
            return str(env_path_resolved)
        # If env var is set but invalid, warn but continue to PATH search
        logger = logging.getLogger("vcfcache")
        logger.warning(
            f"VCFCACHE_BCFTOOLS points to invalid or non-executable path: {env_path_resolved}, "
            "falling back to PATH search"
        )

    # Fall back to system PATH
    return shutil.which("bcftools")


def check_bcftools_version(bcftools_path: str) -> str:
    """Check bcftools version and return it.

    Args:
        bcftools_path: Path to bcftools binary

    Returns:
        Version string

    Raises:
        RuntimeError: If version check fails
    """
    try:
        result = subprocess.run(
            [bcftools_path, "--version-only"],
            capture_output=True,
            text=True,
            check=True,
            timeout=5
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Failed to get bcftools version: {e}")
    except subprocess.TimeoutExpired:
        raise RuntimeError("bcftools version check timed out")


def check_bcftools_installed(min_version: str = MIN_BCFTOOLS_VERSION) -> str:
    """Check if bcftools is installed and meets minimum version requirement.

    This function:
    1. Checks VCFCACHE_BCFTOOLS environment variable first
    2. Falls back to system PATH using shutil.which()
    3. Checks the version is >= min_version
    4. Returns the path to the binary

    Args:
        min_version: Minimum required version (default: 1.20)

    Returns:
        Path to bcftools binary

    Raises:
        FileNotFoundError: If bcftools is not found in PATH
        RuntimeError: If version is too old or version check fails

    Example:
        >>> bcftools_path = check_bcftools_installed()
        >>> print(f"Using bcftools at: {bcftools_path}")
        >>> # Or override system bcftools
        >>> os.environ['VCFCACHE_BCFTOOLS'] = '/opt/bcftools-1.22/bin/bcftools'
        >>> bcftools_path = check_bcftools_installed()
    """
    logger = logging.getLogger("vcfcache")

    # Find bcftools in PATH
    bcftools_path = find_bcftools()

    if not bcftools_path:
        error_msg = (
            "bcftools not found in PATH. "
            f"Please install bcftools >= {min_version}.\n\n"
            "Installation instructions:\n"
            "  - Ubuntu/Debian: sudo apt-get install bcftools\n"
            "  - macOS: brew install bcftools\n"
            "  - Conda: conda install -c bioconda bcftools\n"
            "  - From source: http://www.htslib.org/download/\n"
            "  - OR use the docker image with bundled bcftools\n\n"
            "Alternatively, set VCFCACHE_BCFTOOLS to point to a specific bcftools binary:\n"
            "  export VCFCACHE_BCFTOOLS=/path/to/bcftools"
        )
        logger.error(error_msg)
        raise FileNotFoundError(error_msg)

    # Check version
    try:
        version = check_bcftools_version(bcftools_path)
        logger.debug(f"Found bcftools {version} at {bcftools_path}")

        # Compare versions
        if compare_versions(version, min_version) < 0:
            error_msg = (
                f"bcftools version {version} is too old. "
                f"Minimum required version is {min_version}.\n"
                f"Found at: {bcftools_path}\n\n"
                "Please upgrade bcftools:\n"
                "  - Ubuntu/Debian: sudo apt-get update && sudo apt-get install --only-upgrade bcftools\n"
                "  - macOS: brew upgrade bcftools\n"
                "  - Conda: conda update bcftools\n"
                "  - From source: http://www.htslib.org/download/\n\n"
                "Alternatively, set VCFCACHE_BCFTOOLS to override system bcftools:\n"
                "  export VCFCACHE_BCFTOOLS=/path/to/newer/bcftools"
            )
            logger.error(error_msg)
            raise RuntimeError(error_msg)

        logger.info(f"Using bcftools {version} at {bcftools_path}")
        return bcftools_path

    except RuntimeError:
        raise
    except Exception as e:
        error_msg = f"Error checking bcftools installation: {e}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)



def compute_md5(file_path: Path) -> str:
    """Compute MD5 checksum for a file.

    Args:
        file_path: Path to the file to compute MD5 for

    Returns:
        MD5 checksum as a string

    Example:
        >>> compute_md5(Path('~/projects/vcfcache/tests/data/nodata/dbsnp_test.bcf'))
    """
    import hashlib

    if not file_path.exists():
        raise FileNotFoundError(str(file_path))

    print(f"Computing MD5 for {file_path} ...")
    md5_hash = hashlib.md5()
    with file_path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


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
        required_formats = {"GT", "AD"}  # Removed DP requirement
        available_formats = set(vcf.header.formats.keys())

        missing_formats = required_formats - available_formats
        if missing_formats:
            return (
                False,
                f"Missing required FORMAT fields: {', '.join(missing_formats)}",
            )

        return True, None

    except Exception as e:
        return False, f"Error reading VCF file: {e}"


def generate_test_command(
    vcfcache_path="${VCFCACHE_ROOT}/vcfcache.py",
    vcf_path="${VCFCACHE_ROOT}/tests/data/nodata/crayz_db.bcf",
    output_dir="/tmp/vcfcache/test_cache",
    yaml_path="${VCFCACHE_ROOT}/tests/config/example_params.yaml",
    annotation_config="${VCFCACHE_ROOT}/tests/config/example_annotation.config",
    add_vcf_path="${VCFCACHE_ROOT}/tests/data/nodata/crayz_db2.bcf",
    input_vcf_path="${VCFCACHE_ROOT}/tests/data/nodata/sample4.bcf",
    annotate_name="testor",
    annotation_db="/tmp/vcfcache/test_cache/cache/testor",
    annotation_output="/tmp/vcfcache/aout.bcf",
    annotation_stats_dir="/tmp/vcfcache/aout_stats",
    force=True,
):
    """Generate a nicely formatted test command string for vcfcache operations.

    Returns:
        str: A copy-pastable command string with proper formatting
    """
    cmd_init = (
        f"{vcfcache_path} blueprint-init "
        f"--vcf {vcf_path} "
        f"--output {output_dir} "
        f"-y {yaml_path} "
        f"{'-f' if force else ''} "
    ).strip()

    cmd_add = (
        f"{vcfcache_path} blueprint-extend " f"--db {output_dir} " f"-i {add_vcf_path} "
    ).strip()

    cmd_annotate = (
        f"{vcfcache_path} cache-build "
        f"--name {annotate_name} "
        f"-a {annotation_config} "
        f"--db {output_dir} "
        f"{'-f' if force else ''} "
    ).strip()

    cmd_annotatevcf = (
        f"{vcfcache_path} annotate "
        f"-a {annotation_db} "
        f"--vcf {input_vcf_path} "
        f"{'-f' if force else ''} "
        f"--output {annotation_output} "
        f"--stats-dir {annotation_stats_dir} "
    ).strip()

    # Combine commands
    full_cmd = f"{cmd_init} ; {cmd_add} ; {cmd_annotate} ; {cmd_annotatevcf}"

    # Also create a nicely formatted display version for easier reading
    formatted_cmds = f"""
# INITIALIZE
{cmd_init}

# ADD
{cmd_add}

# CACHE BUILD
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
# ~/projects/vcfcache/vcfcache.py blueprint-init --name nftest --vcf ~/projects/vcfcache/tests/data/nodata/crayz_db.bcf --output /home/j380r/tmp/test/test_out -f --test -vv
# ~/projects/vcfcache/vcfcache.py blueprint-extend --db ~/projects/vcfcache/tests/data/test_out/nftest -i ~/projects/vcfcache/tests/data/nodata/crayz_db2.bcf --test -vv
# ~/projects/vcfcache/vcfcache.py cache-build --name testor --db ~/projects/vcfcache/tests/data/test_out/nftest --test -vv -f
# ... or locally
# ~/projects/vcfcache/vcfcache.py blueprint-init --vcf ~/projects/vcfcache/tests/data/nodata/crayz_db.bcf --output ~/tmp/vcfcache/test_cache -c ~/projects/vcfcache/tests/config/env_test.config -f
# ~/projects/vcfcache/vcfcache.py blueprint-extend --db /home/j380r/tmp/test/test_out -i ~/projects/vcfcache/tests/data/nodata/crayz_db2.bcf
# ~/projects/vcfcache/vcfcache.py cache-build --name testor --db test_out/nftest --test -vv -f
# ~/projects/vcfcache/vcfcache.py annotate --a ~/tmp/test/test_out/nftest/cache/testor --vcf ~/projects/vcfcache/tests/data/nodata/sample4.bcf --output ~/tmp/test/aout.bcf --stats-dir /tmp/test/aout_stats --test -f

# as one:
cmd = """alias stx="
~/projects/vcfcache/vcfcache.py blueprint-init --vcf ~/projects/vcfcache/tests/data/nodata/crayz_db.bcf --output ~/tmp/vcfcache/test_cache -y ~/projects/vcfcache/tests/config/example_params.yaml -f;
~/projects/vcfcache/vcfcache.py blueprint-extend --db ~/tmp/vcfcache/test_cache/ -i ~/projects/vcfcache/tests/data/nodata/crayz_db2.bcf ; 
~/projects/vcfcache/vcfcache.py cache-build --name testor -a ~/projects/vcfcache/tests/config/example_annotation.config --db ~/tmp/vcfcache/test_cache -f;
~/projects/vcfcache/vcfcache.py annotate -a ~/tmp/vcfcache/test_cache/cache/testor --vcf ~/projects/vcfcache/tests/data/nodata/sample4.bcf --output ~/tmp/vcfcache/aout.bcf --stats-dir /tmp/vcfcache/aout_stats -f
"""


# on gvpre
cmd2 = """
~/projects/vcfcache/vcfcache.py blueprint-init --vcf /mnt/data/resources/gnomad/vcf_gnomad_v4_hg19_exomes/gnomad.exomes.v4.1.sites.grch37.trimmed_liftover_norm_1e-1.bcf --output gnomad_1e-1  -c ~/projects/vcfcache/tests/config/nextflow_gnomadhg19.config;
~/projects/vcfcache/vcfcache.py blueprint-extend --db gnomad_1e-1 -i /mnt/data/resources/gnomad/vcf_gnomad_v4_hg19_genomes/gnomad.genomes.v4.1.sites.grch37.trimmed_liftover_norm_1e-1.bcf;
~/projects/vcfcache/vcfcache.py cache-build --name gen_ex -a ~/projects/vcfcache/tests/config/example_annotation.config --db gnomad_1e-1;
~/projects/vcfcache/vcfcache.py annotate -a gnomad_1e-1/cache/gen_ex --vcf /mnt/data/samples/test_mgm/mgm_WGS_32.gatkWGS_norm.bcf --output mgm_WGS_32.bcf --stats-dir mgm_WGS_32_stats -p;
"""
