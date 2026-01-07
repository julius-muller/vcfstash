#!/usr/bin/env python3
"""
update_reference_data.py - Utility script for generating reference (golden) datasets for vcfcache tests.

This script was originally intended to create a golden dataset for regression testing, but was later deprecated for automated test validation due to issues with versioning and reproducibility. However, it remains useful for generating output from the current test files, which can be used to manually compare results across different vcfcache versions or environments.

The script runs the full vcfcache pipeline (blueprint-init, blueprint-extend, cache-build, annotate) using test data, normalizes timestamps for reproducibility, and writes results to a designated output directory. This allows developers to inspect and compare outputs for debugging or validation purposes.
"""

import re
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
from vcfcache.utils.paths import get_vcfcache_root
from vcfcache.utils.validation import check_bcftools_installed

# Use Path for better path handling
TEST_ROOT = Path(os.path.dirname(os.path.abspath(__file__)))
VCFCACHE_CMD = "vcfcache"
TEST_DATA_DIR = TEST_ROOT / "data" / "nodata"
TEST_PARAMS = TEST_ROOT / "config" / "example_params.yaml"
TEST_VCF = TEST_DATA_DIR / "crayz_db.bcf"
EXPECTED_OUTPUT_DIR = TEST_ROOT / "data" / "expected_output"
TEST_ANNO_CONFIG = TEST_ROOT / "config" / "test_annotation.yaml"


def normalize_bcf_timestamps(bcf_file):
    """Normalize timestamps in BCF file to make tests more stable."""
    if not os.path.exists(bcf_file):
        print(f"Warning: BCF file not found at {bcf_file}")
        return

    print(f"Normalizing timestamps in {bcf_file}")

    # Create a temporary VCF file
    temp_vcf = bcf_file + ".temp.vcf"

    # Get bcftools path (respects VCFCACHE_BCFTOOLS, otherwise uses PATH)
    bcftools_path = str(check_bcftools_installed())

    # Convert BCF to VCF
    result = subprocess.run(
        [bcftools_path, "view", bcf_file, "-o", temp_vcf],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    if result.returncode != 0:
        print(f"Error converting BCF to VCF: {result.stderr}")
        if os.path.exists(temp_vcf):
            os.remove(temp_vcf)
        return

    # Read and modify the VCF content
    modified_lines = []
    with open(temp_vcf, 'r') as f:
        for line in f:
            # Replace timestamp patterns in header lines
            if line.startswith('##'):
                line = re.sub(r'\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+\d+\s+\d+:\d+:\d+\s+\d+\b',
                              'Jan 1 00:00:00 2023', line)
                # Also handle fileDate format
                if "##fileDate=" in line:
                    line = "##fileDate=20230101\n"
            modified_lines.append(line)

    # Write modified content back to the temporary file
    with open(temp_vcf, 'w') as f:
        f.writelines(modified_lines)

    # Convert back to BCF
    result = subprocess.run(
        [bcftools_path, "view", "-O", "b", "-o", bcf_file, temp_vcf],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    if result.returncode != 0:
        print(f"Error converting normalized VCF back to BCF: {result.stderr}")

    # Clean up
    if os.path.exists(temp_vcf):
        os.remove(temp_vcf)


def normalize_text_file_timestamps(file_path):
    """Normalize timestamps in text files to make tests more stable."""
    if not os.path.exists(file_path):
        print(f"Warning: File not found at {file_path}")
        return

    print(f"Normalizing timestamps in {file_path}")

    # Read the file content
    with open(file_path, 'r') as f:
        content = f.read()

    # Replace timestamp patterns
    # ISO format: 2023-04-02T18:14:50
    normalized_content = re.sub(
        r'\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}',
        '2023-01-01T00:00:00',
        content
    )

    # Date formats like "Apr 2 18:14:50 2025"
    normalized_content = re.sub(
        r'\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+\d+\s+\d+:\d+:\d+\s+\d+\b',
        'Jan 1 00:00:00 2023',
        normalized_content
    )

    # Write normalized content back
    with open(file_path, 'w') as f:
        f.write(normalized_content)


def update_golden_reference_dataset(force=True):
    """Update the golden reference dataset using test data.

    This function runs all the commands (blueprint-init, blueprint-extend, cache-build, annotate)
    in sequence and uses two output directories for the data. It uses relative paths to
    make it work in any environment.

    Args:
        force: If True, overwrite existing reference data. Defaults to True.

    Returns:
        bool: True if the update was successful, False otherwise.
    """
    print("=== Updating golden reference dataset ===")

    # Use subdirectories in the expected output directory
    cache_dir = os.path.join(EXPECTED_OUTPUT_DIR, "cache_result")
    annotate_file = os.path.join(EXPECTED_OUTPUT_DIR, "annotate_result.bcf")
    annotate_stats_dir = os.path.join(EXPECTED_OUTPUT_DIR, "annotate_stats")

    # Ensure the directories don't exist
    for dir_path in [cache_dir, annotate_stats_dir]:
        if os.path.exists(dir_path):
            if force:
                print(f"Removing existing directory: {dir_path}")
                shutil.rmtree(dir_path)
            else:
                print(f"Directory {dir_path} already exists. Use --force to overwrite.")
                return False
    if os.path.exists(annotate_file):
        if force:
            print(f"Removing existing file: {annotate_file}")
            os.unlink(annotate_file)
        else:
            print(f"File {annotate_file} already exists. Use --force to overwrite.")
            return False

    # Create a temporary params file with the correct paths
    temp_params_file = None
    try:
        # Get the VCFCACHE_ROOT directory
        vcfcache_root = str(get_vcfcache_root())

        # Create a temporary params file with the correct paths
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as temp_file:
            temp_params_file = temp_file.name

            # Read the original params file
            with open(TEST_PARAMS, 'r') as f:
                params_content = f.read()

            # Replace ${VCFCACHE_ROOT} with the actual value
            params_content = params_content.replace('${VCFCACHE_ROOT}', vcfcache_root)

            # Write the modified content to the temporary file
            temp_file.write(params_content)

        # Define the test files
        test_vcf = str(Path(TEST_DATA_DIR) / "crayz_db.bcf")
        test_vcf2 = str(Path(TEST_DATA_DIR) / "crayz_db2.bcf")
        test_sample = str(Path(TEST_DATA_DIR) / "sample4.bcf")

        # Define the annotation name
        annotate_name = "testor"

        # 1. Run blueprint-init
        print("Running blueprint-init...")
        init_cmd = [
            VCFCACHE_CMD,
            "blueprint-init",
            "--vcf", test_vcf,
            "--output", cache_dir,
            "-y", temp_params_file,
            "-f"
        ]

        init_result = subprocess.run(
            init_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if init_result.returncode != 0:
            print(f"blueprint-init failed: {init_result.stderr}")
            return False

        # 2. Run blueprint-extend
        print("Running blueprint-extend...")
        add_cmd = [
            VCFCACHE_CMD,
            "blueprint-extend",
            "--db", cache_dir,
            "-i", test_vcf2
        ]

        add_result = subprocess.run(
            add_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if add_result.returncode != 0:
            print(f"blueprint-extend failed: {add_result.stderr}")
            return False

        # 3. Run cache-build
        print("Running cache-build...")
        annotate_cmd = [
            VCFCACHE_CMD,
            "cache-build",
            "--name", annotate_name,
            "-a", TEST_ANNO_CONFIG,
            "--db", cache_dir,
            "-y", temp_params_file,
            "-f"
        ]

        annotate_result = subprocess.run(
            annotate_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if annotate_result.returncode != 0:
            print(f"cache-build failed: {annotate_result.stderr}")
            return False

        # 4. Run annotate
        print("Running annotate...")
        # Use the annotation directory path
        annotation_db = os.path.join(cache_dir, "cache", annotate_name)

        annotate_vcf_cmd = [
            VCFCACHE_CMD,
            "annotate",
            "-a", annotation_db,
            "--vcf", test_sample,
            "--output", annotate_file,
            "--stats-dir", annotate_stats_dir,
            "-y", temp_params_file,
            "-f"
        ]

        annotate_vcf_result = subprocess.run(
            annotate_vcf_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if annotate_vcf_result.returncode != 0:
            print(f"annotate failed: {annotate_vcf_result.stderr}")
            return False

        # Print the commands that were run (similar to the ones in the issue description)
        print("\nCommands that were run:")
        print(f"{VCFCACHE_CMD} blueprint-init --vcf {test_vcf} --output {cache_dir} -y {temp_params_file} -f")
        print(f"{VCFCACHE_CMD} blueprint-extend --db {cache_dir} -i {test_vcf2}")
        print(f"{VCFCACHE_CMD} cache-build --name {annotate_name} -a {TEST_ANNO_CONFIG} --db {cache_dir} -y {temp_params_file} -f")
        print(f"{VCFCACHE_CMD} annotate -a {annotation_db} --vcf {test_sample} --output {annotate_file} --stats-dir {annotate_stats_dir} -y {temp_params_file} -f")

        print("\nOutput directories:")
        print(f"Cache directory: {cache_dir}")
        print(f"Annotate stats directory: {annotate_stats_dir}")

        return True

    except Exception as e:
        print(f"Error during golden reference dataset update: {str(e)}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        # Don't clean up the temporary directories, as they are the output of the function
        # But do clean up the temporary params file
        if temp_params_file and os.path.exists(temp_params_file):
            os.unlink(temp_params_file)


# Update the main part of the script to include the new function
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Update reference data for vcfcache tests")
    parser.add_argument('--force', action='store_true', help='Force overwrite of existing reference data')
    parser.add_argument('--golden', action='store_true', help='Update golden reference dataset')

    args = parser.parse_args()


    # Update golden reference dataset
    if args.golden:
        success = update_golden_reference_dataset(force=args.force)
        if not success:
            print("Failed to update golden reference dataset")
