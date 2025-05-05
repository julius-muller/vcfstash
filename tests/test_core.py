"""Test core functionality of VCFstash."""

import os
import json
import pytest
from pathlib import Path
import subprocess
import shutil
import tempfile
from vcfstash.utils.paths import get_vcfstash_root
from vcfstash.utils.validation import compute_md5

# Constants
TEST_ROOT = get_vcfstash_root() / "tests"
TEST_DATA_DIR = TEST_ROOT / "data" / "nodata"
TEST_VCF = TEST_DATA_DIR / "crayz_db.bcf"
TEST_VCF2 = TEST_DATA_DIR / "crayz_db2.bcf"
VCFSTASH_CMD = "vcfstash"

@pytest.fixture
def test_output_dir():
    """Provides a temporary directory for testing."""
    temp_dir = tempfile.mkdtemp(prefix="vcfstash_test_")
    os.rmdir(temp_dir)
    yield temp_dir
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir, ignore_errors=True)

def test_basic_workflow(test_output_dir, params_file):
    """Test the basic stash workflow: init and add."""
    # Test stash-init
    init_cmd = [VCFSTASH_CMD, "stash-init", "-i", str(TEST_VCF),
                "-o", test_output_dir, "-y", params_file, "-f"]
    result = subprocess.run(init_cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"stash-init failed: {result.stderr}"

    # Verify directory structure
    blueprint_dir = Path(test_output_dir) / "blueprint"
    assert blueprint_dir.exists(), "Blueprint directory not created"
    assert (blueprint_dir / "vcfstash.bcf").exists(), "VCF file not created"
    assert (blueprint_dir / "sources.info").exists(), "Sources file not created"

    # Test stash-add
    add_cmd = [VCFSTASH_CMD, "stash-add", "--db", test_output_dir,
               "-i", str(TEST_VCF2)]
    result = subprocess.run(add_cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"stash-add failed: {result.stderr}"

    # Verify updated sources.info
    with open(blueprint_dir / "sources.info") as f:
        sources = json.load(f)
        assert len(sources["input_files"]) == 2, "Second file not added"

def test_error_handling(test_output_dir, params_file):
    """Test error conditions and edge cases."""
    # Test with non-existent input file
    init_cmd = [VCFSTASH_CMD, "stash-init", "-i", "nonexistent.bcf",
                "-o", test_output_dir, "-y", params_file]
    result = subprocess.run(init_cmd, capture_output=True, text=True)
    assert result.returncode != 0, "Should fail with non-existent input"

    # Test with invalid output location
    init_cmd = [VCFSTASH_CMD, "stash-init", "-i", str(TEST_VCF),
                "-o", test_output_dir, "-y", "nonexistent.yaml"]
    result = subprocess.run(init_cmd, capture_output=True, text=True)
    assert result.returncode != 0, "Should fail with invalid yaml"

    # Test add without init
    add_cmd = [VCFSTASH_CMD, "stash-add", "--db", test_output_dir,
               "-i", str(TEST_VCF)]
    result = subprocess.run(add_cmd, capture_output=True, text=True)
    assert result.returncode != 0, "Should fail without initialization"

def test_file_validation(test_output_dir: str):
    """Test file validation and integrity checks."""
    # Create test file
    ref_file = TEST_ROOT / "data/references/reference.fasta"

    # Test MD5 calculation
    md5_hash = compute_md5(ref_file)
    assert md5_hash == "ec59e3976d29e276414191a6283499f7"
