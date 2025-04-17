#!/usr/bin/env python3
"""
update_reference_data.py - Script to update reference data for tests
"""

import re
import sys
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import uuid

# Use Path for better path handling
TEST_ROOT = Path(os.path.dirname(os.path.abspath(__file__)))
PROJECT_ROOT = TEST_ROOT.parent
VCFSTASH_CMD = str(PROJECT_ROOT / "vcfstash.py")
TEST_DATA_DIR = str(TEST_ROOT / "data" / "nodata")
TEST_CONFIG = str(TEST_ROOT / "config" / "env_test.config")
TEST_PARAMS = str(TEST_ROOT / "config" / "user_params.yaml")
TEST_VCF = str(Path(TEST_DATA_DIR) / "crayz_db.bcf")
EXPECTED_OUTPUT_DIR = str(TEST_ROOT / "data" / "expected_output")
TEST_ANNO_CONFIG = str(TEST_ROOT / "config" / "annotation.config")



def normalize_bcf_timestamps(bcf_file):
    """Normalize timestamps in BCF file to make tests more stable."""
    if not os.path.exists(bcf_file):
        print(f"Warning: BCF file not found at {bcf_file}")
        return

    print(f"Normalizing timestamps in {bcf_file}")

    # Create a temporary VCF file
    temp_vcf = bcf_file + ".temp.vcf"

    # Convert BCF to VCF
    result = subprocess.run(
        ["bcftools", "view", bcf_file, "-o", temp_vcf],
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
        ["bcftools", "view", "-O", "b", "-o", bcf_file, temp_vcf],
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

    This function runs all the commands (stash-init, stash-add, stash-annotate, annotate)
    in sequence and uses two output directories for the data. It uses relative paths to
    make it work in any environment.

    Args:
        force: If True, overwrite existing reference data. Defaults to True.

    Returns:
        bool: True if the update was successful, False otherwise.
    """
    print("=== Updating golden reference dataset ===")

    # Use subdirectories in the expected output directory
    stash_dir = os.path.join(EXPECTED_OUTPUT_DIR, "stash_result")
    annotate_dir = os.path.join(EXPECTED_OUTPUT_DIR, "annotate_result")

    # Ensure the directories don't exist
    for dir_path in [stash_dir, annotate_dir]:
        if os.path.exists(dir_path):
            if force:
                print(f"Removing existing directory: {dir_path}")
                shutil.rmtree(dir_path)
            else:
                print(f"Directory {dir_path} already exists. Use --force to overwrite.")
                return False

    try:
        # Define the test files
        test_vcf = str(Path(TEST_DATA_DIR) / "crayz_db.bcf")
        test_vcf2 = str(Path(TEST_DATA_DIR) / "crayz_db2.bcf")
        test_sample = str(Path(TEST_DATA_DIR) / "sample4.bcf")

        # Define the annotation name
        annotate_name = "testor"

        # 1. Run stash-init
        print("Running stash-init...")
        init_cmd = [
            sys.executable,
            VCFSTASH_CMD,
            "stash-init",
            "--vcf", test_vcf,
            "--output", stash_dir,
            "-y", TEST_PARAMS,
            "-f"
        ]

        init_result = subprocess.run(
            init_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if init_result.returncode != 0:
            print(f"stash-init failed: {init_result.stderr}")
            return False

        # 2. Run stash-add
        print("Running stash-add...")
        add_cmd = [
            sys.executable,
            VCFSTASH_CMD,
            "stash-add",
            "--db", stash_dir,
            "-i", test_vcf2
        ]

        add_result = subprocess.run(
            add_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if add_result.returncode != 0:
            print(f"stash-add failed: {add_result.stderr}")
            return False

        # 3. Run stash-annotate
        print("Running stash-annotate...")
        annotate_cmd = [
            sys.executable,
            VCFSTASH_CMD,
            "stash-annotate",
            "--name", annotate_name,
            "-a", TEST_ANNO_CONFIG,
            "--db", stash_dir,
            "-y", TEST_PARAMS,
            "-f"
        ]

        annotate_result = subprocess.run(
            annotate_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        if annotate_result.returncode != 0:
            print(f"stash-annotate failed: {annotate_result.stderr}")
            return False

        # 4. Run annotate
        print("Running annotate...")
        # Use the annotation directory path
        annotation_db = os.path.join(stash_dir, "stash", annotate_name)

        annotate_vcf_cmd = [
            sys.executable,
            VCFSTASH_CMD,
            "annotate",
            "-a", annotation_db,
            "--vcf", test_sample,
            "--output", annotate_dir,
            "-y", TEST_PARAMS,
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
        print(f"{VCFSTASH_CMD} stash-init --vcf {test_vcf} --output {stash_dir} -y {TEST_PARAMS} -f")
        print(f"{VCFSTASH_CMD} stash-add --db {stash_dir} -i {test_vcf2}")
        print(f"{VCFSTASH_CMD} stash-annotate --name {annotate_name} -a {TEST_ANNO_CONFIG} --db {stash_dir} -y {TEST_PARAMS} -f")
        print(f"{VCFSTASH_CMD} annotate -a {annotation_db} --vcf {test_sample} --output {annotate_dir} -y {TEST_PARAMS} -f")

        print("\nOutput directories:")
        print(f"Stash directory: {stash_dir}")
        print(f"Annotate directory: {annotate_dir}")

        return True

    except Exception as e:
        print(f"Error during golden reference dataset update: {str(e)}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        # Don't clean up the temporary directories, as they are the output of the function
        pass


# Update the main part of the script to include the new function
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Update reference data for vcfstash tests")
    parser.add_argument('--force', action='store_true', help='Force overwrite of existing reference data')
    parser.add_argument('--golden', action='store_true', help='Update golden reference dataset')

    args = parser.parse_args()


    # Update golden reference dataset
    if args.golden:
        success = update_golden_reference_dataset(force=args.force)
        if not success:
            print("Failed to update golden reference dataset")
