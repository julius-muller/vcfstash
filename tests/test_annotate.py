"""Test annotate functionality of VCFstash."""

import os
import sys
import json
import pytest
from pathlib import Path
import subprocess
import shutil
import tempfile
from vcfstash.utils.paths import get_vcfstash_root, get_resource_path

# Constants
TEST_ROOT = get_vcfstash_root() / "tests"
TEST_DATA_DIR = TEST_ROOT / "data" / "nodata"
TEST_VCF = TEST_DATA_DIR / "crayz_db.bcf"
TEST_VCF2 = TEST_DATA_DIR / "crayz_db2.bcf"
TEST_SAMPLE = TEST_DATA_DIR / "sample4.bcf"
TEST_PARAMS = TEST_ROOT / "config" / "test_params.yaml"
TEST_ANNO_CONFIG = TEST_ROOT / "config" / "test_annotation.config"
VCFSTASH_CMD = "vcfstash"


# Fixture for temporary test output directory
@pytest.fixture
def test_output_dir():
    """Provides a path for a directory that doesn't exist yet."""
    # Create a path for a temporary directory, but don't create it
    temp_dir = tempfile.mkdtemp(prefix="vcfstash_test_")

    # Remove the directory immediately - vcfstash will create it
    os.rmdir(temp_dir)

    # Return the path to the test function
    yield temp_dir

    # Clean up after the test is done
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir, ignore_errors=True)


def create_temp_params_file():
    """Create a temporary params file with the correct paths and return its path."""
    # Get the VCFSTASH_ROOT directory
    vcfstash_root = str(get_vcfstash_root())

    # Create a temporary params file with the correct paths
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
    temp_params_file = temp_file.name

    # Read the original params file
    with open(TEST_PARAMS, 'r') as f:
        params_content = f.read()

    # Replace ${VCFSTASH_ROOT} with the actual value
    params_content = params_content.replace('${VCFSTASH_ROOT}', vcfstash_root)

    # Write the modified content to the temporary file
    temp_file.write(params_content)
    temp_file.close()

    return temp_params_file


def run_stash_init(input_vcf, output_dir, force=False):
    """Run the stash-init command and return the process result."""
    # Make sure the directory doesn't exist (clean start)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    # Create a temporary params file
    temp_params_file = create_temp_params_file()

    try:
        cmd = [
            str(VCFSTASH_CMD),
            "stash-init",
            "--vcf", str(input_vcf),
            "--output", str(output_dir),
            "-y", temp_params_file
        ]

        if force:
            cmd.append("-f")

        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return result
    finally:
        # Clean up the temporary file
        if os.path.exists(temp_params_file):
            os.unlink(temp_params_file)


def run_stash_add(db_dir, input_vcf):
    """Run the stash-add command and return the process result."""
    cmd = [
        str(VCFSTASH_CMD),
        "stash-add",
        "--db", str(db_dir),
        "-i", str(input_vcf)
    ]

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    return result


def run_stash_annotate(db_dir, name, force=False):
    """Run the stash-annotate command and return the process result."""
    # Create a temporary params file
    temp_params_file = create_temp_params_file()

    try:
        cmd = [
            str(VCFSTASH_CMD),
            "stash-annotate",
            "--name", name,
            "-a", str(TEST_ANNO_CONFIG),
            "--db", str(db_dir),
            "-y", temp_params_file
        ]

        if force:
            cmd.append("-f")
        #print(f"RRRRRRRRRRRRRRRunning stash-annotate with commands: {cmd}")
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return result
    except Exception as e:
        print(f"Error running stash-annotate: {e}\nRuuning commands: {cmd}")
        raise e
    finally:
        # Clean up the temporary file
        if os.path.exists(temp_params_file):
            os.unlink(temp_params_file)


def run_annotate(annotation_db, input_vcf, output_dir, force=False):
    """Run the annotate command and return the process result."""
    # Create a temporary params file
    temp_params_file = create_temp_params_file()

    try:
        cmd = [
            str(VCFSTASH_CMD),
            "annotate",
            "-a", str(annotation_db),
            "--vcf", str(input_vcf),
            "--output", str(output_dir),
            "-y", temp_params_file
        ]

        if force:
            cmd.append("-f")

        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return result
    finally:
        # Clean up the temporary file
        if os.path.exists(temp_params_file):
            os.unlink(temp_params_file)


def test_sample_file_validity(test_output_dir):
    """Test that the sample BCF file is valid."""
    print("\n=== Testing sample file validity ===")

    # Get bcftools path for verification
    bcftools_path = get_resource_path('tools/bcftools')
    if not bcftools_path.exists():
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Check if the sample BCF file exists
    assert TEST_SAMPLE.exists(), f"Sample BCF file not found: {TEST_SAMPLE}"
    print(f"Sample file exists: {TEST_SAMPLE}")

    # Check if the sample BCF file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", str(TEST_SAMPLE)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Sample BCF file is not valid: {view_result.stderr}"
    print("Sample file has valid header")

    # Check if the sample BCF file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", str(TEST_SAMPLE)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert stats_result.returncode == 0, f"Failed to get stats for sample BCF file: {stats_result.stderr}"
    assert "number of records:" in stats_result.stdout, "Sample BCF file has no variants"

    # Extract the number of records
    num_records = 0
    for line in stats_result.stdout.splitlines():
        if "number of records:" in line:
            num_records = int(line.split(":")[-1].strip())
            break

    print(f"Sample file has {num_records} variants")
    print("Successfully verified sample file validity")



def test_full_annotation_workflow(test_output_dir):
    """Test the full annotation workflow from stash-init to annotate."""
    print("\n=== Testing full annotation workflow ===")

    # Step 1: Run stash-init
    print("Running stash-init...")
    init_result = run_stash_init(TEST_VCF, test_output_dir, force=True)
    assert init_result.returncode == 0, f"stash-init failed: {init_result.stderr}"

    # Step 2: Run stash-add
    print("Running stash-add...")
    add_result = run_stash_add(test_output_dir, TEST_VCF2)
    assert add_result.returncode == 0, f"stash-add failed: {add_result.stderr}"

    # Print information about the workflow directory and files
    workflow_dir = os.path.join(test_output_dir, "workflow")
    print(f"Workflow directory exists: {os.path.exists(workflow_dir)}")
    if os.path.exists(workflow_dir):
        print(f"Workflow directory contents: {os.listdir(workflow_dir)}")

    # Step 3: Run stash-annotate
    print("Running stash-annotate...")
    annotate_name = "test_annotation"
    annotate_result = run_stash_annotate(test_output_dir, annotate_name, force=True)
    if annotate_result.returncode != 0:
        print(f"Command output: {annotate_result.stdout}")
        print(f"Command error: {annotate_result.stderr}")
        print(f"Working directory contents: {os.listdir(test_output_dir)}")
        print(f"Workflow directory contents: {os.listdir(workflow_dir)}")
    assert annotate_result.returncode == 0, f"stash-annotate failed: {annotate_result.stderr}"

    # Get bcftools path for verification
    bcftools_path = get_resource_path('tools/bcftools')
    if not bcftools_path.exists():
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Step 4: Verify the annotation directory was created
    stash_dir = os.path.join(test_output_dir, "stash")
    annotation_dir = os.path.join(stash_dir, annotate_name)
    assert os.path.exists(annotation_dir), f"Annotation directory not found: {annotation_dir}"
    print(f"Annotation directory created: {annotation_dir}")

    # Step 5: Create output directory for annotate
    output_dir = os.path.join(test_output_dir, "full_workflow_output")

    # Step 6: Run annotate
    print("Running annotate...")
    annotate_result = run_annotate(annotation_dir, TEST_SAMPLE, output_dir, force=True)
    if annotate_result.returncode != 0:
        print(f"Command output: {annotate_result.stdout}")
        print(f"Command error: {annotate_result.stderr}")
        print(f"Working directory contents: {os.listdir(test_output_dir)}")
        print(f"Workflow directory contents: {os.listdir(workflow_dir)}")
    assert annotate_result.returncode == 0, f"annotate failed: {annotate_result.stderr}"

    # Step 7: Verify the output directory exists
    assert os.path.exists(output_dir), f"Output directory not found: {output_dir}"
    print(f"Output directory created: {output_dir}")

    # Step 8: Verify the output file exists
    output_file = os.path.join(output_dir, os.path.splitext(os.path.basename(str(TEST_SAMPLE)))[0] + "_vst.bcf")
    if not os.path.exists(output_file):
        print(
            f"Files in output_dir ({output_dir}): {os.listdir(output_dir) if os.path.exists(output_dir) else 'Directory does not exist'}")
        raise FileNotFoundError(f"Output file not found: {output_file}")
    print(f"Output file created: {output_file}")

    # Step 9: Verify the output file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Output file is not valid: {view_result.stderr}"
    print("Output file has valid header")

    # Step 10: Verify the MOCK_ANNO tag is present in the header
    assert "MOCK_ANNO" in view_result.stdout, "MOCK_ANNO tag not found in the header"
    print("MOCK_ANNO tag found in header")

    # Step 11: Verify the output file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert stats_result.returncode == 0, f"Failed to get stats for output file: {stats_result.stderr}"
    assert "number of records:" in stats_result.stdout, "Output file has no variants"

    # Extract the number of records
    num_records = 0
    for line in stats_result.stdout.splitlines():
        if "number of records:" in line:
            num_records = int(line.split(":")[-1].strip())
            break

    print(f"Output file has {num_records} variants")

    # Step 12: Verify the MOCK_ANNO tag is present in the variants
    variants_result = subprocess.run(
        [str(bcftools_path), "view", output_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert variants_result.returncode == 0, f"Failed to view output file: {variants_result.stderr}"
    assert "MOCK_ANNO=" in variants_result.stdout, "MOCK_ANNO tag not found in the variants"
    print("MOCK_ANNO tag found in variants")

    print("Successfully tested full annotation workflow")


def test_cached_vs_uncached_annotation(test_output_dir, params_file):
    """Test that cached and uncached annotations produce identical results."""
    print("\n=== Testing cached vs uncached annotation ===")

    # Step 1: Create a database
    print("Creating database...")
    init_cmd = [VCFSTASH_CMD, "stash-init", "-i", str(TEST_VCF),
                "-o", test_output_dir, "-y", params_file, "-f"]
    init_result = subprocess.run(init_cmd, capture_output=True, text=True)
    assert init_result.returncode == 0, f"stash-init failed: {init_result.stderr}"

    # Step 2: Run stash-annotate to create the annotation cache
    print("Creating annotation cache...")
    annotate_name = "test_annotation"
    stash_annotate_cmd = [VCFSTASH_CMD, "stash-annotate", "--name", annotate_name,
                          "--db", test_output_dir, "-a", str(TEST_ANNO_CONFIG),
                          "-y", params_file, "-f"]
    stash_annotate_result = subprocess.run(stash_annotate_cmd, capture_output=True, text=True)
    assert stash_annotate_result.returncode == 0, f"stash-annotate failed: {stash_annotate_result.stderr}"

    # Step 3: Run annotation with caching
    print("Running cached annotation...")
    cached_output = f"{test_output_dir}/cached_output"
    cached_cmd = [VCFSTASH_CMD, "annotate", "-a", f"{test_output_dir}/stash/{annotate_name}",
                  "-i", str(TEST_SAMPLE), "-o", cached_output,
                  "-y", params_file, "-f"]
    cached_result = subprocess.run(cached_cmd, capture_output=True, text=True)
    assert cached_result.returncode == 0, f"Cached annotation failed: {cached_result.stderr}"

    # Step 4: Run annotation without caching
    print("Running uncached annotation...")
    uncached_output = f"{test_output_dir}/uncached_output"
    uncached_cmd = [VCFSTASH_CMD, "annotate", "-a", f"{test_output_dir}/stash/{annotate_name}",
                    "-i", str(TEST_SAMPLE), "-o", uncached_output,
                    "-y", params_file, "--uncached", "-f"]
    uncached_result = subprocess.run(uncached_cmd, capture_output=True, text=True)
    assert uncached_result.returncode == 0, f"Uncached annotation failed: {uncached_result.stderr}"

    # Step 5: Compare the outputs
    print("Comparing outputs...")
    
    # Get bcftools path for comparison
    bcftools_path = get_resource_path('tools/bcftools')
    if not bcftools_path.exists():
        bcftools_path = 'bcftools'

    # Compare headers
    cached_header = subprocess.run(
        [str(bcftools_path), "view", "-h", f"{cached_output}/annotated.bcf"],
        capture_output=True, text=True
    )
    uncached_header = subprocess.run(
        [str(bcftools_path), "view", "-h", f"{uncached_output}/annotated.bcf"],
        capture_output=True, text=True
    )
    assert cached_header.stdout == uncached_header.stdout, "Headers differ between cached and uncached outputs"

    # Compare variants
    cached_variants = subprocess.run(
        [str(bcftools_path), "view", f"{cached_output}/annotated.bcf"],
        capture_output=True, text=True
    )
    uncached_variants = subprocess.run(
        [str(bcftools_path), "view", f"{uncached_output}/annotated.bcf"],
        capture_output=True, text=True
    )

    # Sort and compare variant lines
    cached_lines = sorted(cached_variants.stdout.splitlines())
    uncached_lines = sorted(uncached_variants.stdout.splitlines())
    
    # Compare line by line
    for i, (cached_line, uncached_line) in enumerate(zip(cached_lines, uncached_lines)):
        assert cached_line == uncached_line, f"Variant mismatch at line {i+1}:\nCached:   {cached_line}\nUncached: {uncached_line}"

    # Verify we have the same number of variants
    assert len(cached_lines) == len(uncached_lines), f"Different number of variants: cached={len(cached_lines)}, uncached={len(uncached_lines)}"

    print("Successfully verified that cached and uncached annotations produce identical results")
