"""Test annotation functionality of VCFstash."""

import os
import sys
import pytest
from pathlib import Path
import subprocess
import shutil
import tempfile
import time

# Constants
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TEST_VCF = os.path.join(TEST_DATA_DIR, "nodata", "crayz_db.bcf")
TEST_CONFIG = os.path.join(os.path.dirname(__file__), "config", "env_test.config")
TEST_PARAMS = os.path.join(os.path.dirname(__file__), "config", "user_params.yaml")
TEST_ANNO_CONFIG = os.path.join(os.path.dirname(__file__), "config", "annotation.config")
VCFSTASH_CMD = os.path.join(os.path.dirname(os.path.dirname(__file__)), "vcfstash.py")
EXPECTED_OUTPUT_DIR = os.path.join(TEST_DATA_DIR, "expected_output")

def get_variants(bcf_path):
    """Read BCF contents excluding headers."""
    try:
        bcf_text = subprocess.run(
            ["bcftools", "view", str(bcf_path)],
            capture_output=True,
            text=True,
            check=True
        ).stdout
        return [line for line in bcf_text.splitlines()
                if line and not line.startswith('#')]
    except subprocess.CalledProcessError as e:
        pytest.fail(f"Failed to read BCF file {bcf_path}: {e.stderr}")
        return []

def test_cached_vs_uncached_annotation():
    """Test that the annotated BCF file exists and has content."""
    # The update_reference.py script now creates only one BCF file directly in the annotate_result directory
    bcf_file = Path(EXPECTED_OUTPUT_DIR) / "annotate_result" / "sample4_vst.bcf"

    assert bcf_file.exists(), f"Annotated BCF not found at {bcf_file}"

    variants = get_variants(bcf_file)
    assert len(variants) > 0, f"No variants found in BCF {bcf_file}"

    # Print the number of variants found
    print(f"Found {len(variants)} variants in {bcf_file}")

def test_cached_vs_uncached_variants_identical():
    """Test that variants from cached and uncached annotation are identical."""
    # Skip if reference files don't exist
    test_input = Path(TEST_DATA_DIR) / "nodata" / "sample4.bcf"
    test_cache = Path(EXPECTED_OUTPUT_DIR) / "stash_result/stash/testor"

    if not test_input.exists() or not test_cache.exists():
        pytest.skip(f"Required test files not found: {test_input} or {test_cache}")

    # Create temporary output directories for cached and uncached runs
    with tempfile.TemporaryDirectory() as temp_dir:
        cached_output_dir = Path(temp_dir) / "cached"
        uncached_output_dir = Path(temp_dir) / "uncached"

        # Run with cache
        cached_cmd = [
            sys.executable, VCFSTASH_CMD,
            "annotate",
            "-i", str(test_input),
            "-a", str(test_cache),
            "-o", str(cached_output_dir),
            "-y", TEST_PARAMS,
            "-f"
        ]

        # Run without cache
        uncached_cmd = [
            sys.executable, VCFSTASH_CMD,
            "annotate",
            "-i", str(test_input),
            "-a", str(test_cache),
            "-o", str(uncached_output_dir),
            "-y", TEST_PARAMS,
            "--uncached",
            "-f"
        ]

        # Execute both commands
        try:
            print("Running annotation with cache...")
            cached_result = subprocess.run(
                cached_cmd, 
                capture_output=True, 
                text=True, 
                check=True,
                timeout=300
            )

            print("Running annotation without cache...")
            uncached_result = subprocess.run(
                uncached_cmd, 
                capture_output=True, 
                text=True, 
                check=True,
                timeout=300
            )

            # Check that output files exist
            cached_output_file = cached_output_dir / "sample4_vst.bcf"
            uncached_output_file = uncached_output_dir / "sample4_vst.bcf"

            assert cached_output_file.exists(), f"Cached output file not created: {cached_output_file}"
            assert uncached_output_file.exists(), f"Uncached output file not created: {uncached_output_file}"

            # Get variants from both files
            cached_variants = get_variants(cached_output_file)
            uncached_variants = get_variants(uncached_output_file)

            # Check that both files have variants
            assert len(cached_variants) > 0, f"No variants found in cached BCF {cached_output_file}"
            assert len(uncached_variants) > 0, f"No variants found in uncached BCF {uncached_output_file}"

            # Print the number of variants found
            print(f"Found {len(cached_variants)} variants in cached BCF")
            print(f"Found {len(uncached_variants)} variants in uncached BCF")

            # Compare variant counts
            assert len(cached_variants) == len(uncached_variants), (
                f"Number of variants differs: cached={len(cached_variants)}, "
                f"uncached={len(uncached_variants)}"
            )

            # Compare each variant
            for i, (cached, uncached) in enumerate(zip(cached_variants, uncached_variants)):
                assert cached == uncached, (
                    f"Variants differ at line {i + 1}:\n"
                    f"Cached:   {cached}\n"
                    f"Uncached: {uncached}"
                )

            print("All variants are identical between cached and uncached annotation.")

        except subprocess.CalledProcessError as e:
            pytest.fail(f"Command failed with exit code {e.returncode}:\n{e.stderr}")
        except subprocess.TimeoutExpired:
            pytest.fail("Command timed out after 5 minutes")

@pytest.mark.parametrize("use_cache", [True, False])
def test_annotate_command(use_cache):
    """Test that the annotate command works correctly with and without cache."""
    # Skip if reference files don't exist
    test_input = Path(TEST_DATA_DIR) / "nodata" / "sample4.bcf"
    test_cache = Path(EXPECTED_OUTPUT_DIR) / "stash_result/stash/testor"

    if not test_input.exists() or not test_cache.exists():
        pytest.skip(f"Required test files not found: {test_input} or {test_cache}")

    # Create temporary output directory
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir) / "output"
        # Don't create the directory - let vcfstash create it

        # Build command
        cmd = [
            sys.executable, VCFSTASH_CMD,
            "annotate",
            "-i", str(test_input),
            "-a", str(test_cache),
            "-o", str(output_dir),
            "-y", TEST_PARAMS,
            "-f"
        ]

        # Add uncached flag if testing without cache
        if not use_cache:
            cmd.append("--uncached")

        # Run command
        try:
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                check=True,
                timeout=300
            )

            # Check that output file exists
            output_file = output_dir / "sample4_vst.bcf"
            assert output_file.exists(), f"Output file not created: {output_file}"

            # Check that file has content
            file_size = output_file.stat().st_size
            assert file_size > 0, f"Output file is empty: {output_file}"

        except subprocess.CalledProcessError as e:
            pytest.fail(f"Command failed with exit code {e.returncode}:\n{e.stderr}")
        except subprocess.TimeoutExpired:
            pytest.fail("Command timed out after 5 minutes")
