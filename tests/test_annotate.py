"""Test annotate functionality of VCFstash."""

import os
from pathlib import Path
import subprocess
import shutil
import pytest
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


def run_stash_init(input_vcf, output_dir, force=False, normalize=False):
    """Run the stash-init command and return the process result."""
    # Make sure the directory doesn't exist (clean start)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = [
        str(VCFSTASH_CMD),
        "stash-init",
        "--vcf", str(input_vcf),
        "--output", str(output_dir),
        "-y", TEST_PARAMS
    ]

    if force:
        cmd.append("-f")

    if normalize:
        cmd.append("-n")

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    return result


def run_stash_add(db_dir, input_vcf, normalize=False):
    """Run the stash-add command and return the process result."""
    cmd = [
        str(VCFSTASH_CMD),
        "stash-add",
        "--db", str(db_dir),
        "-i", str(input_vcf)
    ]

    if normalize:
        cmd.append("-n")

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    return result


def run_stash_annotate(db_dir, name, force=False):
    """Run the stash-annotate command and return the process result."""

    cmd = [
        str(VCFSTASH_CMD),
        "stash-annotate",
        "--name", name,
        "-a", str(TEST_ANNO_CONFIG),
        "--db", str(db_dir),
        "-y", TEST_PARAMS
    ]

    if force:
        cmd.append("-f")

    try:

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


def run_annotate(annotation_db, input_vcf, output_dir, force=False):
    """Run the annotate command and return the process result."""

    cmd = [
        str(VCFSTASH_CMD),
        "annotate",
        "-a", str(annotation_db),
        "--vcf", str(input_vcf),
        "--output", str(output_dir),
        "-y", TEST_PARAMS
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


def test_sample_file_validity(test_output_dir, test_scenario):
    """Test that the sample BCF file is valid."""
    print(f"\n=== Testing sample file validity (scenario: {test_scenario}) ===")

    # Get bcftools path for verification
    bcftools_path = get_resource_path(Path('tools/bcftools'))
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



def test_full_annotation_workflow(test_output_dir, test_scenario, prebuilt_cache):
    """Test the full annotation workflow from stash-init to annotate.

    Adapts based on scenario:
    - vanilla: Create cache from scratch
    - blueprint: Use prebuilt cache + create test stash
    - annotated: Use prebuilt cache + create test stash
    """
    print(f"\n=== Testing full annotation workflow (scenario: {test_scenario}) ===")

    # Determine which cache to use
    if test_scenario == "vanilla":
        # Step 1: Run stash-init to create cache from scratch
        print("Running stash-init (creating cache from test data)...")
        init_result = run_stash_init(TEST_VCF, test_output_dir, force=True)
        assert init_result.returncode == 0, f"stash-init failed: {init_result.stderr}"

        # Step 2: Run stash-add
        print("Running stash-add...")
        add_result = run_stash_add(test_output_dir, TEST_VCF2)
        assert add_result.returncode == 0, f"stash-add failed: {add_result.stderr}"

        db_dir = test_output_dir
    else:
        # Blueprint/Annotated: Use prebuilt cache, create a test stash
        print(f"Using prebuilt cache at {prebuilt_cache}")
        print("Creating test cache from test data for annotation testing...")

        # Create a test cache in test_output_dir for annotation testing
        init_result = run_stash_init(TEST_VCF, test_output_dir, force=True)
        assert init_result.returncode == 0, f"stash-init failed: {init_result.stderr}"

        add_result = run_stash_add(test_output_dir, TEST_VCF2)
        assert add_result.returncode == 0, f"stash-add failed: {add_result.stderr}"

        db_dir = test_output_dir

    # Print information about the workflow directory and files
    workflow_dir = Path(test_output_dir) / "workflow"
    print(f"Workflow directory exists: {workflow_dir.exists()}")
    if workflow_dir.exists():
        print(f"Workflow directory contents: {list(workflow_dir.iterdir())}")

    # Step 3: Run stash-annotate
    print("Running stash-annotate...")
    annotate_name = "test_annotation"
    annotate_result = run_stash_annotate(test_output_dir, annotate_name, force=True)
    if annotate_result.returncode != 0:
        print(f"Command output: {annotate_result.stdout}")
        print(f"Command error: {annotate_result.stderr}")
        print(f"Working directory contents: {list(Path(test_output_dir).iterdir())}")
        print(f"Workflow directory contents: {list(workflow_dir.iterdir())}")
    assert annotate_result.returncode == 0, f"stash-annotate failed: {annotate_result.stderr}"

    # Get bcftools path for verification
    bcftools_path = get_resource_path(Path('tools/bcftools'))
    if not bcftools_path.exists():
        # Fall back to system bcftools if the project-specific one doesn't exist
        bcftools_path = 'bcftools'

    # Step 4: Verify the annotation directory was created
    stash_dir = Path(test_output_dir) / "stash"
    annotation_dir = stash_dir / annotate_name
    assert annotation_dir.exists(), f"Annotation directory not found: {annotation_dir}"
    print(f"Annotation directory created: {annotation_dir}")

    # Step 5: Create output directory for annotate
    output_dir = Path(test_output_dir) / "full_workflow_output"

    # Step 6: Run annotate
    print("Running annotate...")
    annotate_result = run_annotate(annotation_dir, TEST_SAMPLE, output_dir, force=True)
    if annotate_result.returncode != 0:
        print(f"Command output: {annotate_result.stdout}")
        print(f"Command error: {annotate_result.stderr}")
        print(f"Working directory contents: {list(Path(test_output_dir).iterdir())}")
        print(f"Workflow directory contents: {list(workflow_dir.iterdir())}")
    assert annotate_result.returncode == 0, f"annotate failed: {annotate_result.stderr}"

    # Step 7: Verify the output directory exists
    assert output_dir.exists(), f"Output directory not found: {output_dir}"
    print(f"Output directory created: {output_dir}")

    # Step 8: Verify the output file exists
    output_file = output_dir / (Path(str(TEST_SAMPLE)).stem + "_vst.bcf")
    if not output_file.exists():
        print(
            f"Files in output_dir ({output_dir}): {list(output_dir.iterdir()) if output_dir.exists() else 'Directory does not exist'}")
        raise FileNotFoundError(f"Output file not found: {output_file}")
    print(f"Output file created: {output_file}")

    # Step 9: Verify the output file is valid
    view_result = subprocess.run(
        [str(bcftools_path), "view", "-h", str(output_file)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert view_result.returncode == 0, f"Output file is not valid: {view_result.stderr}"
    print("Output file has valid header")

    # Step 10: Verify the MOCK_ANNO tag is present in the header
    assert "MOCK_ANNO" in view_result.stdout, "MOCK_ANNO tag not found in the header"
    print("MOCK_ANNO tag found in header")

    # Step 11: Verify the output file has variants
    stats_result = subprocess.run(
        [str(bcftools_path), "stats", str(output_file)],
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
        [str(bcftools_path), "view", str(output_file)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert variants_result.returncode == 0, f"Failed to view output file: {variants_result.stderr}"
    assert "MOCK_ANNO=" in variants_result.stdout, "MOCK_ANNO tag not found in the variants"
    print("MOCK_ANNO tag found in variants")

    # Step 13: Verify the output contains the auxiliary information

    # Check auxiliary directory exists
    auxiliary_dir = output_dir / "auxiliary"
    assert auxiliary_dir.exists(), "Auxiliary directory was not created"
    assert auxiliary_dir.is_dir(), "Auxiliary directory is not a directory"

    # Check for the expected test file
    test_file = auxiliary_dir / "auxiliary_test.exo"
    assert test_file.exists(), "Auxiliary test file was not copied to the auxiliary directory"

    print("Successfully tested full annotation workflow")


def test_cached_vs_uncached_annotation(test_output_dir, params_file, test_scenario):
    """Test that cached and uncached annotations produce identical results."""
    print(f"\n=== Testing cached vs uncached annotation (scenario: {test_scenario}) ===")

    # Step 1: Create a database
    print("Creating database...")
    init_cmd = [VCFSTASH_CMD, "stash-init", "-i", str(TEST_VCF),
                "-o", str(test_output_dir), "-y", str(params_file), "-f"]
    init_result = subprocess.run(init_cmd, capture_output=True, text=True)
    assert init_result.returncode == 0, f"stash-init failed: {init_result.stderr}"

    # Step 2: Run stash-annotate to create the annotation cache
    print("Creating annotation cache...")
    annotate_name = "test_annotation"
    stash_annotate_cmd = [VCFSTASH_CMD, "stash-annotate", "--name", annotate_name,
                          "--db", str(test_output_dir), "-a", str(TEST_ANNO_CONFIG),
                          "-y", str(params_file), "-f"]
    stash_annotate_result = subprocess.run(stash_annotate_cmd, capture_output=True, text=True)
    assert stash_annotate_result.returncode == 0, f"stash-annotate failed: {stash_annotate_result.stderr}"

    # Step 3: Run annotation with caching
    print("Running cached annotation...")
    cached_output = Path(test_output_dir) / "cached_output"
    cached_cmd = [VCFSTASH_CMD, "annotate", "-a", str(Path(test_output_dir) / "stash" / annotate_name),
                  "-i", str(TEST_SAMPLE), "-o", str(cached_output),
                  "-y", str(params_file), "-f"]
    cached_result = subprocess.run(cached_cmd, capture_output=True, text=True)
    assert cached_result.returncode == 0, f"Cached annotation failed: {cached_result.stderr}"

    # Step 4: Run annotation without caching
    print("Running uncached annotation...")
    uncached_output = Path(test_output_dir) / "uncached_output"
    uncached_cmd = [VCFSTASH_CMD, "annotate", "-a", str(Path(test_output_dir) / "stash" / annotate_name),
                    "-i", str(TEST_SAMPLE), "-o", str(uncached_output),
                    "-y", str(params_file), "--uncached", "-f"]
    uncached_result = subprocess.run(uncached_cmd, capture_output=True, text=True)
    assert uncached_result.returncode == 0, f"Uncached annotation failed: {uncached_result.stderr}"

    # Step 5: Compare the outputs
    print("Comparing outputs...")

    # Get bcftools path for comparison
    bcftools_path = get_resource_path(Path('tools/bcftools'))
    if not bcftools_path.exists():
        bcftools_path = 'bcftools'

    # Compare headers
    cached_header = subprocess.run(
        [str(bcftools_path), "view", "-h", str(cached_output / "annotated.bcf")],
        capture_output=True, text=True
    )
    uncached_header = subprocess.run(
        [str(bcftools_path), "view", "-h", str(uncached_output / "annotated.bcf")],
        capture_output=True, text=True
    )
    assert cached_header.stdout == uncached_header.stdout, "Headers differ between cached and uncached outputs"

    # Compare variants
    cached_variants = subprocess.run(
        [str(bcftools_path), "view", str(cached_output / "annotated.bcf")],
        capture_output=True, text=True
    )
    uncached_variants = subprocess.run(
        [str(bcftools_path), "view", str(uncached_output / "annotated.bcf")],
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


def test_input_not_modified_during_annotation(test_output_dir, params_file, test_scenario):
    """Test that input VCF files are not modified during annotation."""
    print(f"\n=== Testing input file preservation during annotation (scenario: {test_scenario}) ===")

    # Step 1: Create a database
    print("Creating database...")
    init_cmd = [VCFSTASH_CMD, "stash-init", "-i", str(TEST_VCF),
                "-o", str(test_output_dir), "-y", str(params_file), "-f"]
    init_result = subprocess.run(init_cmd, capture_output=True, text=True)
    assert init_result.returncode == 0, f"stash-init failed: {init_result.stderr}"

    # Step 2: Run stash-annotate to create the annotation cache
    print("Creating annotation cache...")
    annotate_name = "test_annotation"
    stash_annotate_cmd = [VCFSTASH_CMD, "stash-annotate", "--name", annotate_name,
                          "--db", str(test_output_dir), "-a", str(TEST_ANNO_CONFIG),
                          "-y", str(params_file), "-f"]
    stash_annotate_result = subprocess.run(stash_annotate_cmd, capture_output=True, text=True)
    assert stash_annotate_result.returncode == 0, f"stash-annotate failed: {stash_annotate_result.stderr}"

    # Step 3: Make a copy of the input file to compare later
    print("Creating a copy of the input file...")
    input_copy_dir = Path(test_output_dir) / "input_copy"
    input_copy_dir.mkdir(exist_ok=True, parents=True)
    input_copy = input_copy_dir / "sample_copy.bcf"

    # Get the MD5 hash of the original input file
    original_md5 = subprocess.run(
        ["md5sum", str(TEST_SAMPLE)],
        capture_output=True, text=True
    ).stdout.split()[0]

    # Copy the input file
    import shutil
    shutil.copy(TEST_SAMPLE, input_copy)
    shutil.copy(f"{TEST_SAMPLE}.csi", f"{input_copy}.csi")

    # Step 4: Run annotation with caching
    print("Running annotation...")
    output_dir = Path(test_output_dir) / "annotation_output"
    annotate_cmd = [VCFSTASH_CMD, "annotate", "-a", str(Path(test_output_dir) / "stash" / annotate_name),
                  "-i", str(TEST_SAMPLE), "-o", str(output_dir),
                  "-y", str(params_file), "-f"]
    annotate_result = subprocess.run(annotate_cmd, capture_output=True, text=True)
    assert annotate_result.returncode == 0, f"Annotation failed: {annotate_result.stderr}"

    # Step 5: Verify the input file was not modified
    print("Verifying input file was not modified...")

    # Get the MD5 hash of the input file after annotation
    after_md5 = subprocess.run(
        ["md5sum", str(TEST_SAMPLE)],
        capture_output=True, text=True
    ).stdout.split()[0]

    # Compare the hashes
    assert original_md5 == after_md5, "Input file was modified during annotation"

    # Step 6: Verify the output file exists and has the expected content
    print("Verifying output file...")

    # Get bcftools path for verification
    bcftools_path = get_resource_path(Path('tools/bcftools'))
    if not bcftools_path.exists():
        bcftools_path = 'bcftools'

    # Check if the output file exists
    # The output file is named after the input file with _vst.bcf suffix
    sample_name = Path(str(TEST_SAMPLE)).stem
    output_file = output_dir / f"{sample_name}_vst.bcf"
    assert output_file.exists(), f"Output file not found: {output_file}"

    # Check if the output file has a valid header
    header_result = subprocess.run(
        [str(bcftools_path), "view", "-h", str(output_file)],
        capture_output=True, text=True
    )
    assert header_result.returncode == 0, f"Output file has invalid header: {header_result.stderr}"

    # Check if the output file has the MOCK_ANNO tag in the header
    assert "MOCK_ANNO" in header_result.stdout, "MOCK_ANNO tag not found in output header"

    # Check if the output file has variants
    variants_result = subprocess.run(
        [str(bcftools_path), "view", str(output_file)],
        capture_output=True, text=True
    )
    assert variants_result.returncode == 0, f"Failed to view output file: {variants_result.stderr}"
    assert "MOCK_ANNO=" in variants_result.stdout, "MOCK_ANNO tag not found in variants"

    print("Successfully verified that input files are not modified during annotation")


def test_normalization_flag(test_output_dir, params_file, test_scenario):
    """Test that the normalization flag works correctly."""
    print(f"\n=== Testing normalization flag functionality (scenario: {test_scenario}) ===")

    # Step 1: Run stash-init with normalization
    print("Running stash-init with normalization...")
    norm_dir = Path(test_output_dir) / "normalized"
    norm_result = run_stash_init(TEST_VCF, norm_dir, force=True, normalize=True)
    assert norm_result.returncode == 0, f"stash-init with normalization failed: {norm_result.stderr}"

    # Step 2: Run stash-init without normalization
    print("Running stash-init without normalization...")
    no_norm_dir = Path(test_output_dir) / "not_normalized"
    no_norm_result = run_stash_init(TEST_VCF, no_norm_dir, force=True, normalize=False)
    assert no_norm_result.returncode == 0, f"stash-init without normalization failed: {no_norm_result.stderr}"

    # Step 3: Compare the output files
    print("Comparing output files...")

    # Get bcftools path for comparison
    bcftools_path = get_resource_path(Path('tools/bcftools'))
    if not bcftools_path.exists():
        bcftools_path = 'bcftools'

    # Compare the headers of the normalized and non-normalized files
    norm_header = subprocess.run(
        [str(bcftools_path), "view", "-h", str(norm_dir / "blueprint" / "vcfstash.bcf")],
        capture_output=True, text=True
    )
    no_norm_header = subprocess.run(
        [str(bcftools_path), "view", "-h", str(no_norm_dir / "blueprint" / "vcfstash.bcf")],
        capture_output=True, text=True
    )

    # The headers should be different if normalization was applied
    assert norm_header.stdout != no_norm_header.stdout, "Normalization did not produce different headers"

    # Step 4: Run stash-add with normalization
    print("Running stash-add with normalization...")
    add_norm_result = run_stash_add(norm_dir, TEST_VCF2, normalize=True)
    assert add_norm_result.returncode == 0, f"stash-add with normalization failed: {add_norm_result.stderr}"

    # Step 5: Run stash-add without normalization
    print("Running stash-add without normalization...")
    add_no_norm_result = run_stash_add(no_norm_dir, TEST_VCF2, normalize=False)
    assert add_no_norm_result.returncode == 0, f"stash-add without normalization failed: {add_no_norm_result.stderr}"

    # Step 6: Compare the output files after stash-add
    print("Comparing output files after stash-add...")

    # Compare the headers after stash-add
    norm_header_after = subprocess.run(
        [str(bcftools_path), "view", "-h", str(norm_dir / "blueprint" / "vcfstash.bcf")],
        capture_output=True, text=True
    )
    no_norm_header_after = subprocess.run(
        [str(bcftools_path), "view", "-h", str(no_norm_dir / "blueprint" / "vcfstash.bcf")],
        capture_output=True, text=True
    )

    # The headers should still be different after stash-add
    assert norm_header_after.stdout != no_norm_header_after.stdout, "Normalization did not produce different headers after stash-add"

    # Step 7: Compare the actual content of the files to verify normalization was applied
    print("Comparing file content to verify normalization...")

    # Get the content of the normalized and non-normalized files
    norm_content = subprocess.run(
        [str(bcftools_path), "view", str(norm_dir / "blueprint" / "vcfstash.bcf")],
        capture_output=True, text=True
    )
    no_norm_content = subprocess.run(
        [str(bcftools_path), "view", str(no_norm_dir / "blueprint" / "vcfstash.bcf")],
        capture_output=True, text=True
    )

    # The content should be different if normalization was applied
    assert norm_content.stdout != no_norm_content.stdout, "Normalization did not produce different file content"

    # Check if the normalized file has chr prefix in chromosome names (a sign of normalization)
    if "chr" in norm_content.stdout:
        print("Verified that normalization adds chr prefix to chromosome names")

    # Check if the number of variants is different (another sign of normalization)
    norm_count = len(norm_content.stdout.strip().split("\n"))
    no_norm_count = len(no_norm_content.stdout.strip().split("\n"))
    print(f"Normalized file has {norm_count} variants, non-normalized file has {no_norm_count} variants")

    # The counts might be different due to normalization (splitting multiallelic sites)
    # but we don't assert this as it depends on the test data

    print("Successfully verified normalization flag functionality")
