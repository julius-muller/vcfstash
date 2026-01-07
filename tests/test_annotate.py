"""Test annotate functionality of VCFcache."""

import os
from pathlib import Path
import subprocess
import shutil
import pytest
import random
from vcfcache.utils.paths import get_vcfcache_root, get_resource_path
from vcfcache.utils.validation import compute_md5

# Constants
TEST_ROOT = get_vcfcache_root() / "tests"
TEST_DATA_DIR = TEST_ROOT / "data" / "nodata"
TEST_VCF = TEST_DATA_DIR / "crayz_db.bcf"
TEST_VCF2 = TEST_DATA_DIR / "crayz_db2.bcf"
TEST_SAMPLE = TEST_DATA_DIR / "sample4.bcf"
TEST_PARAMS = TEST_ROOT / "config" / "test_params.yaml"
TEST_ANNO_CONFIG = TEST_ROOT / "config" / "test_annotation.yaml"
import sys

# Use python -m vcfcache to ensure we use the installed package
VCFCACHE_CMD = [sys.executable, "-m", "vcfcache"]
VCFCACHE_ROOT = get_vcfcache_root()


def _env():
    env = os.environ.copy()
    env["VCFCACHE_ROOT"] = str(VCFCACHE_ROOT)
    return env


def extract_canary_variants(cache_bcf, output_bcf, num_variants=5, seed=42):
    """Extract random canary variants from cache for validation.

    Args:
        cache_bcf: Path to the cache BCF file
        output_bcf: Path where canary variants should be written
        num_variants: Number of random variants to extract
        seed: Random seed for reproducibility

    Returns:
        List of variant IDs (CHROM:POS:REF:ALT) for the extracted variants
    """
    # Get total number of variants
    stats_result = subprocess.run(
        ["bcftools", "stats", str(cache_bcf)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
    )

    total_variants = 0
    for line in stats_result.stdout.splitlines():
        if "number of records:" in line:
            total_variants = int(line.split(":")[-1].strip())
            break

    if total_variants == 0:
        raise ValueError(f"No variants found in {cache_bcf}")

    # Select random variant indices
    random.seed(seed)
    num_to_extract = min(num_variants, total_variants)
    selected_indices = sorted(random.sample(range(total_variants), num_to_extract))

    # Extract header
    header_result = subprocess.run(
        ["bcftools", "view", "-h", str(cache_bcf)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
    )

    # Get all variants
    variants_result = subprocess.run(
        ["bcftools", "view", "-H", str(cache_bcf)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
    )

    all_variants = variants_result.stdout.strip().split("\n")
    selected_variants = [all_variants[i] for i in selected_indices]

    # Create VCF with selected variants
    vcf_content = header_result.stdout + "\n".join(selected_variants) + "\n"

    # Write to temporary VCF, then convert to BCF
    temp_vcf = str(output_bcf).replace(".bcf", ".vcf")
    with open(temp_vcf, 'w') as f:
        f.write(vcf_content)

    # Convert to BCF and index
    subprocess.run(
        ["bcftools", "view", "-Ob", "-o", str(output_bcf), temp_vcf],
        check=True
    )
    subprocess.run(
        ["bcftools", "index", str(output_bcf)],
        check=True
    )

    # Extract variant IDs for comparison
    variant_ids = []
    for variant_line in selected_variants:
        fields = variant_line.split("\t")
        chrom, pos, _, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        variant_ids.append(f"{chrom}:{pos}:{ref}:{alt}")

    # Clean up temp VCF
    os.remove(temp_vcf)

    return variant_ids


def compare_info_tag_values(bcf1, bcf2, info_tag, variant_ids):
    """Compare INFO tag values between two BCF files for specific variants.

    Args:
        bcf1: Path to first BCF file
        bcf2: Path to second BCF file
        info_tag: Name of the INFO tag to compare (e.g., 'CSQ', 'MOCK_ANNO')
        variant_ids: List of variant IDs to compare (CHROM:POS:REF:ALT format)

    Returns:
        dict: Comparison results with 'matches', 'mismatches', and 'details'
    """
    # Extract INFO tag values from both files
    query_format = f"%CHROM:%POS:%REF:%ALT\\t%INFO/{info_tag}\\n"

    result1 = subprocess.run(
        ["bcftools", "query", "-f", query_format, str(bcf1)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
    )

    result2 = subprocess.run(
        ["bcftools", "query", "-f", query_format, str(bcf2)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
    )

    # Parse results into dictionaries
    def parse_query_output(output):
        values = {}
        for line in output.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) == 2:
                var_id, tag_value = parts
                values[var_id] = tag_value
        return values

    values1 = parse_query_output(result1.stdout)
    values2 = parse_query_output(result2.stdout)

    # Compare values for specified variants
    matches = []
    mismatches = []

    for var_id in variant_ids:
        val1 = values1.get(var_id, "")
        val2 = values2.get(var_id, "")

        if val1 == val2:
            matches.append(var_id)
        else:
            mismatches.append({
                'variant': var_id,
                'cached_value': val1,
                'fresh_value': val2
            })

    return {
        'matches': matches,
        'mismatches': mismatches,
        'total': len(variant_ids),
        'match_count': len(matches),
        'mismatch_count': len(mismatches)
    }


def run_blueprint_init(input_vcf, output_dir, force=False, normalize=False):
    """Run the blueprint-init command and return the process result."""
    # Make sure the directory doesn't exist (clean start)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = VCFCACHE_CMD + [
        "blueprint-init",
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
        text=True,
        env=_env(),
    )
    return result


def run_blueprint_extend(db_dir, input_vcf, normalize=False):
    """Run the blueprint-extend command and return the process result."""
    cmd = VCFCACHE_CMD + [
        "blueprint-extend",
        "--db", str(db_dir),
        "-i", str(input_vcf)
    ]

    if normalize:
        cmd.append("-n")

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=_env(),
    )
    return result


def run_cache_build(db_dir, name, force=False):
    """Run the cache-build command and return the process result."""

    cmd = VCFCACHE_CMD + [
        "cache-build",
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
            text=True,
            env=_env(),
        )
        return result

    except Exception as e:
        print(f"Error running cache-build: {e}\nRuuning commands: {cmd}")
        raise e


def run_annotate(annotation_db, input_vcf, output_dir, force=False):
    """Run the annotate command and return the process result."""

    cmd = VCFCACHE_CMD + [
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
        text=True,
        env=_env(),
    )
    return result


def test_sample_file_validity(test_output_dir, test_scenario):
    """Test that the sample BCF file is valid."""
    print(f"\n=== Testing sample file validity (scenario: {test_scenario}) ===")

    # Use bcftools from PATH (respects setup_test_environment fixture)
    # In annotated images, this will be /opt/bcftools/bin/bcftools (compiled 1.22)
    # In other scenarios, this will be the bundled or system bcftools
    from tests.conftest import get_bcftools_cmd
    bcftools_path = get_bcftools_cmd()

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
    """Test the full annotation workflow from blueprint-init to annotate.

    Adapts based on scenario:
    - vanilla: Create cache from scratch
    - blueprint: Use prebuilt cache + create test cache
    - annotated: Use prebuilt cache + create test cache
    """
    print(f"\n=== Testing full annotation workflow (scenario: {test_scenario}) ===")

    # Determine which cache to use
    if test_scenario == "vanilla":
        # Step 1: Run blueprint-init to create cache from scratch
        print("Running blueprint-init (creating cache from test data)...")
        init_result = run_blueprint_init(TEST_VCF, test_output_dir, force=True)
        assert init_result.returncode == 0, f"blueprint-init failed: {init_result.stderr}"

        # Step 2: Run blueprint-extend
        print("Running blueprint-extend...")
        add_result = run_blueprint_extend(test_output_dir, TEST_VCF2)
        assert add_result.returncode == 0, f"blueprint-extend failed: {add_result.stderr}"

        db_dir = test_output_dir
    else:
        # Blueprint/Annotated: Use prebuilt cache, create a test cache
        print(f"Using prebuilt cache at {prebuilt_cache}")
        print("Creating test cache from test data for annotation testing...")

        # Create a test cache in test_output_dir for annotation testing
        init_result = run_blueprint_init(TEST_VCF, test_output_dir, force=True)
        assert init_result.returncode == 0, f"blueprint-init failed: {init_result.stderr}"

        add_result = run_blueprint_extend(test_output_dir, TEST_VCF2)
        assert add_result.returncode == 0, f"blueprint-extend failed: {add_result.stderr}"

        db_dir = test_output_dir

    # Print information about the workflow directory and files
    workflow_dir = Path(test_output_dir) / "workflow"
    print(f"Workflow directory exists: {workflow_dir.exists()}")
    if workflow_dir.exists():
        print(f"Workflow directory contents: {list(workflow_dir.iterdir())}")

    # Step 3: Run cache-build
    print("Running cache-build...")
    annotate_name = "test_annotation"
    annotate_result = run_cache_build(test_output_dir, annotate_name, force=True)
    if annotate_result.returncode != 0:
        print(f"Command output: {annotate_result.stdout}")
        print(f"Command error: {annotate_result.stderr}")
        print(f"Working directory contents: {list(Path(test_output_dir).iterdir())}")
        print(f"Workflow directory contents: {list(workflow_dir.iterdir())}")
    assert annotate_result.returncode == 0, f"cache-build failed: {annotate_result.stderr}"

    # Use bcftools from PATH (respects setup_test_environment fixture)
    # In annotated images, this will be /opt/bcftools/bin/bcftools (compiled 1.22)
    # In other scenarios, this will be the bundled or system bcftools
    from tests.conftest import get_bcftools_cmd
    bcftools_path = get_bcftools_cmd()

    # Step 4: Verify the annotation directory was created
    cache_dir = Path(test_output_dir) / "cache"
    annotation_dir = cache_dir / annotate_name
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
    output_file = output_dir / (Path(str(TEST_SAMPLE)).stem + "_vc.bcf")
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
    init_cmd = VCFCACHE_CMD + [ "blueprint-init", "-i", str(TEST_VCF),
                "-o", str(test_output_dir), "-y", str(params_file), "-f"]
    init_result = subprocess.run(init_cmd, capture_output=True, text=True)
    assert init_result.returncode == 0, f"blueprint-init failed: {init_result.stderr}"

    # Step 2: Run cache-build to create the annotation cache
    print("Creating annotation cache...")
    annotate_name = "test_annotation"
    cache_build_cmd = VCFCACHE_CMD + [ "cache-build", "--name", annotate_name,
                          "--db", str(test_output_dir), "-a", str(TEST_ANNO_CONFIG),
                          "-y", str(params_file), "-f"]
    cache_build_result = subprocess.run(cache_build_cmd, capture_output=True, text=True)
    assert cache_build_result.returncode == 0, f"cache-build failed: {cache_build_result.stderr}"

    # Step 3: Run annotation with caching
    print("Running cached annotation...")
    cached_output = Path(test_output_dir) / "cached_output"
    cached_cmd = VCFCACHE_CMD + [ "annotate", "-a", str(Path(test_output_dir) / "cache" / annotate_name),
                  "-i", str(TEST_SAMPLE), "-o", str(cached_output),
                  "-y", str(params_file), "-f"]
    cached_result = subprocess.run(cached_cmd, capture_output=True, text=True)
    assert cached_result.returncode == 0, f"Cached annotation failed: {cached_result.stderr}"

    # Step 4: Run annotation without caching
    print("Running uncached annotation...")
    uncached_output = Path(test_output_dir) / "uncached_output"
    uncached_cmd = VCFCACHE_CMD + [ "annotate", "-a", str(Path(test_output_dir) / "cache" / annotate_name),
                    "-i", str(TEST_SAMPLE), "-o", str(uncached_output),
                    "-y", str(params_file), "--uncached", "-f"]
    uncached_result = subprocess.run(uncached_cmd, capture_output=True, text=True)
    assert uncached_result.returncode == 0, f"Uncached annotation failed: {uncached_result.stderr}"

    # Step 5: Compare the outputs
    print("Comparing outputs...")

    # Use bcftools from PATH (respects setup_test_environment fixture)
    from tests.conftest import get_bcftools_cmd
    bcftools_path = get_bcftools_cmd()

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
    init_cmd = VCFCACHE_CMD + [ "blueprint-init", "-i", str(TEST_VCF),
                "-o", str(test_output_dir), "-y", str(params_file), "-f"]
    init_result = subprocess.run(init_cmd, capture_output=True, text=True)
    assert init_result.returncode == 0, f"blueprint-init failed: {init_result.stderr}"

    # Step 2: Run cache-build to create the annotation cache
    print("Creating annotation cache...")
    annotate_name = "test_annotation"
    cache_build_cmd = VCFCACHE_CMD + [ "cache-build", "--name", annotate_name,
                          "--db", str(test_output_dir), "-a", str(TEST_ANNO_CONFIG),
                          "-y", str(params_file), "-f"]
    cache_build_result = subprocess.run(cache_build_cmd, capture_output=True, text=True)
    assert cache_build_result.returncode == 0, f"cache-build failed: {cache_build_result.stderr}"

    # Step 3: Make a copy of the input file to compare later
    print("Creating a copy of the input file...")
    input_copy_dir = Path(test_output_dir) / "input_copy"
    input_copy_dir.mkdir(exist_ok=True, parents=True)
    input_copy = input_copy_dir / "sample_copy.bcf"

    # Get the MD5 hash of the original input file
    original_md5 = compute_md5(TEST_SAMPLE)

    # Copy the input file
    import shutil
    shutil.copy(TEST_SAMPLE, input_copy)
    shutil.copy(f"{TEST_SAMPLE}.csi", f"{input_copy}.csi")

    # Step 4: Run annotation with caching
    print("Running annotation...")
    output_dir = Path(test_output_dir) / "annotation_output"
    annotate_cmd = VCFCACHE_CMD + [ "annotate", "-a", str(Path(test_output_dir) / "cache" / annotate_name),
                  "-i", str(TEST_SAMPLE), "-o", str(output_dir),
                  "-y", str(params_file), "-f"]
    annotate_result = subprocess.run(annotate_cmd, capture_output=True, text=True)
    assert annotate_result.returncode == 0, f"Annotation failed: {annotate_result.stderr}"

    # Step 5: Verify the input file was not modified
    print("Verifying input file was not modified...")

    # Get the MD5 hash of the input file after annotation
    after_md5 = compute_md5(TEST_SAMPLE)

    # Compare the hashes
    assert original_md5 == after_md5, "Input file was modified during annotation"

    # Step 6: Verify the output file exists and has the expected content
    print("Verifying output file...")

    # Use bcftools from PATH (respects setup_test_environment fixture)
    from tests.conftest import get_bcftools_cmd
    bcftools_path = get_bcftools_cmd()

    # Check if the output file exists
    # The output file is named after the input file with _vc.bcf suffix
    sample_name = Path(str(TEST_SAMPLE)).stem
    output_file = output_dir / f"{sample_name}_vc.bcf"
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

    # Step 1: Run blueprint-init with normalization
    print("Running blueprint-init with normalization...")
    norm_dir = Path(test_output_dir) / "normalized"
    norm_result = run_blueprint_init(TEST_VCF, norm_dir, force=True, normalize=True)
    assert norm_result.returncode == 0, f"blueprint-init with normalization failed: {norm_result.stderr}"

    # Step 2: Run blueprint-init without normalization
    print("Running blueprint-init without normalization...")
    no_norm_dir = Path(test_output_dir) / "not_normalized"
    no_norm_result = run_blueprint_init(TEST_VCF, no_norm_dir, force=True, normalize=False)
    assert no_norm_result.returncode == 0, f"blueprint-init without normalization failed: {no_norm_result.stderr}"

    # Step 3: Compare the output files
    print("Comparing output files...")

    # Use bcftools from PATH (respects setup_test_environment fixture)
    from tests.conftest import get_bcftools_cmd
    bcftools_path = get_bcftools_cmd()

    # Compare the headers of the normalized and non-normalized files
    norm_header = subprocess.run(
        [str(bcftools_path), "view", "-h", str(norm_dir / "blueprint" / "vcfcache.bcf")],
        capture_output=True, text=True
    )
    no_norm_header = subprocess.run(
        [str(bcftools_path), "view", "-h", str(no_norm_dir / "blueprint" / "vcfcache.bcf")],
        capture_output=True, text=True
    )

    # The headers should be different if normalization was applied
    assert norm_header.stdout != no_norm_header.stdout, "Normalization did not produce different headers"

    # Step 4: Run blueprint-extend with normalization
    print("Running blueprint-extend with normalization...")
    add_norm_result = run_blueprint_extend(norm_dir, TEST_VCF2, normalize=True)
    assert add_norm_result.returncode == 0, f"blueprint-extend with normalization failed: {add_norm_result.stderr}"

    # Step 5: Run blueprint-extend without normalization
    print("Running blueprint-extend without normalization...")
    add_no_norm_result = run_blueprint_extend(no_norm_dir, TEST_VCF2, normalize=False)
    assert add_no_norm_result.returncode == 0, f"blueprint-extend without normalization failed: {add_no_norm_result.stderr}"

    # Step 6: Compare the output files after blueprint-extend
    print("Comparing output files after blueprint-extend...")

    # Compare the headers after blueprint-extend
    norm_header_after = subprocess.run(
        [str(bcftools_path), "view", "-h", str(norm_dir / "blueprint" / "vcfcache.bcf")],
        capture_output=True, text=True
    )
    no_norm_header_after = subprocess.run(
        [str(bcftools_path), "view", "-h", str(no_norm_dir / "blueprint" / "vcfcache.bcf")],
        capture_output=True, text=True
    )

    # The headers should still be different after blueprint-extend
    assert norm_header_after.stdout != no_norm_header_after.stdout, "Normalization did not produce different headers after blueprint-extend"

    # Step 7: Compare the actual content of the files to verify normalization was applied
    print("Comparing file content to verify normalization...")

    # Get the content of the normalized and non-normalized files
    norm_content = subprocess.run(
        [str(bcftools_path), "view", str(norm_dir / "blueprint" / "vcfcache.bcf")],
        capture_output=True, text=True
    )
    no_norm_content = subprocess.run(
        [str(bcftools_path), "view", str(no_norm_dir / "blueprint" / "vcfcache.bcf")],
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


@pytest.mark.skipif(
    os.environ.get("SKIP_CANARY_TEST", "0") == "1",
    reason="Canary validation test skipped via SKIP_CANARY_TEST=1"
)
def test_canary_validation(test_output_dir, params_file, test_scenario):
    """Test that cached annotations match fresh annotations for canary variants.

    This critical test validates that the annotation tool produces identical results
    when run on the same variants, ensuring:
    - The annotation pipeline is working correctly
    - The cached annotations are valid
    - The annotation tool configuration hasn't changed

    The test:
    1. Creates a cache with annotated variants
    2. Extracts random "canary" variants from the cache
    3. Runs those same variants through the annotation pipeline again
    4. Compares INFO tag values to ensure they're identical

    Can be skipped via: SKIP_CANARY_TEST=1 pytest
    """
    import yaml

    print(f"\n=== Testing canary validation (scenario: {test_scenario}) ===")

    # Step 1: Create a database and cache
    print("Creating database and cache...")
    init_cmd = VCFCACHE_CMD + [ "blueprint-init", "-i", str(TEST_VCF),
                "-o", str(test_output_dir), "-y", str(params_file), "-f"]
    init_result = subprocess.run(init_cmd, capture_output=True, text=True)
    assert init_result.returncode == 0, f"blueprint-init failed: {init_result.stderr}"

    # Step 2: Run cache-build to create the annotation cache
    print("Creating annotation cache...")
    annotate_name = "canary_test_annotation"
    cache_build_cmd = VCFCACHE_CMD + [ "cache-build", "--name", annotate_name,
                       "--db", str(test_output_dir), "-a", str(TEST_ANNO_CONFIG),
                       "-y", str(params_file), "-f"]
    cache_build_result = subprocess.run(cache_build_cmd, capture_output=True, text=True)
    assert cache_build_result.returncode == 0, f"cache-build failed: {cache_build_result.stderr}"

    # Step 3: Load annotation config to get the INFO tag we need to validate
    with open(TEST_ANNO_CONFIG, 'r') as f:
        annotation_config = yaml.safe_load(f)

    info_tag = annotation_config['must_contain_info_tag']
    print(f"Validating INFO tag: {info_tag}")

    # Step 4: Extract canary variants from the annotated cache
    cache_dir = Path(test_output_dir) / "cache" / annotate_name
    annotated_cache_bcf = cache_dir / "vcfcache_annotated.bcf"

    assert annotated_cache_bcf.exists(), f"Annotated cache not found: {annotated_cache_bcf}"

    canary_bcf = Path(test_output_dir) / "canary_variants.bcf"
    print("Extracting canary variants from cache...")
    variant_ids = extract_canary_variants(
        annotated_cache_bcf,
        canary_bcf,
        num_variants=5,  # Test with 5 random variants
        seed=42  # Reproducible selection
    )

    print(f"Extracted {len(variant_ids)} canary variants:")
    for var_id in variant_ids:
        print(f"  - {var_id}")

    # Step 5: Strip annotations from canary variants to create fresh input
    print("Creating unannotated canary variants...")
    unannotated_canary = Path(test_output_dir) / "canary_unannotated.bcf"

    # Remove all INFO tags to simulate fresh input
    subprocess.run(
        ["bcftools", "annotate", "-x", "INFO", "-Ob", "-o", str(unannotated_canary),
         str(canary_bcf)],
        check=True
    )
    subprocess.run(
        ["bcftools", "index", str(unannotated_canary)],
        check=True
    )

    # Step 6: Run the annotation command on the unannotated canaries
    print("Running annotation pipeline on canary variants...")

    # Read the annotation command from the config
    annotation_cmd = annotation_config['annotation_cmd']

    # Load params for variable substitution
    with open(params_file, 'r') as f:
        params_content = f.read()
        # Replace VCFCACHE_ROOT
        params_content = params_content.replace('${VCFCACHE_ROOT}', str(VCFCACHE_ROOT))
    params = yaml.safe_load(params_content)

    # Substitute variables in annotation command
    fresh_annotated = Path(test_output_dir) / "canary_fresh_annotated.bcf"
    auxiliary_dir = Path(test_output_dir) / "auxiliary"
    auxiliary_dir.mkdir(exist_ok=True)

    cmd = annotation_cmd
    cmd = cmd.replace('${INPUT_BCF}', str(unannotated_canary))
    cmd = cmd.replace('${OUTPUT_BCF}', str(fresh_annotated))
    cmd = cmd.replace('${AUXILIARY_DIR}', str(auxiliary_dir))

    # Substitute params variables
    for key, value in params.items():
        if isinstance(value, str):
            cmd = cmd.replace(f'${{params.{key}}}', value)

    # Execute the annotation command
    result = subprocess.run(
        cmd,
        shell=True,
        cwd=test_output_dir,
        capture_output=True,
        text=True,
        env=_env()
    )

    if result.returncode != 0:
        print(f"Annotation command failed!")
        print(f"Command: {cmd}")
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        assert False, f"Fresh annotation failed: {result.stderr}"

    assert fresh_annotated.exists(), f"Fresh annotation output not found: {fresh_annotated}"

    # Step 7: Compare INFO tag values between cached and fresh annotations
    print(f"\nComparing {info_tag} values between cached and fresh annotations...")

    comparison = compare_info_tag_values(
        annotated_cache_bcf,
        fresh_annotated,
        info_tag,
        variant_ids
    )

    print(f"\nComparison results:")
    print(f"  Total variants: {comparison['total']}")
    print(f"  Matches: {comparison['match_count']}")
    print(f"  Mismatches: {comparison['mismatch_count']}")

    if comparison['mismatches']:
        print(f"\nMismatched variants:")
        for mismatch in comparison['mismatches']:
            print(f"  Variant: {mismatch['variant']}")
            print(f"    Cached:  {mismatch['cached_value']}")
            print(f"    Fresh:   {mismatch['fresh_value']}")

    # Assert all variants match
    assert comparison['mismatch_count'] == 0, \
        f"Canary validation failed: {comparison['mismatch_count']} variants have different annotations"

    print(f"\n✓ Canary validation passed: All {comparison['match_count']} variants have identical annotations")
    print("Successfully validated that cached annotations match fresh annotations")
