"""Test core functionality of VCFstash."""

from pathlib import Path
import subprocess
import pytest
from vcfstash.utils.paths import get_vcfstash_root
from vcfstash.utils.validation import compute_md5

# Constants
TEST_ROOT = get_vcfstash_root() / "tests"
TEST_DATA_DIR = TEST_ROOT / "data" / "nodata"
TEST_VCF = TEST_DATA_DIR / "crayz_db.bcf"
TEST_VCF2 = TEST_DATA_DIR / "crayz_db2.bcf"
VCFSTASH_CMD = "vcfstash"


def test_error_handling(test_output_dir, params_file, test_scenario):
    """Test error conditions and edge cases."""
    print(f"\n=== Testing error handling (scenario: {test_scenario}) ===")

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

def test_file_validation(test_output_dir: str, test_scenario):
    """Test file validation and integrity checks."""
    print(f"\n=== Testing file validation (scenario: {test_scenario}) ===")

    # Create test file
    ref_file = TEST_ROOT / "data/references/reference.fasta"

    # Test MD5 calculation
    md5_hash = compute_md5(ref_file)
    assert md5_hash == "ec59e3976d29e276414191a6283499f7"

def test_vcf_reference_validation(test_scenario):
    """Test VCF reference validation."""
    print(f"\n=== Testing VCF reference validation (scenario: {test_scenario}) ===")

    from vcfstash.database.base import VCFDatabase
    import logging

    # Set up test files
    vcf_file = TEST_DATA_DIR / "crayz_db.bcf"
    ref_file = TEST_ROOT / "data/references/reference.fasta"

    # Create a VCFDatabase instance
    db = VCFDatabase(Path(TEST_ROOT), 2, True, TEST_ROOT.parent / "tools" / "bcftools")
    db.logger = logging.getLogger("test")

    # Test validation with valid files
    result, error = db.validate_vcf_reference(vcf_file, ref_file)
    assert result, f"Validation should pass but failed with: {error}"

    # Test validation with non-existent VCF file
    result, error = db.validate_vcf_reference(Path("nonexistent.bcf"), ref_file)
    assert not result, "Validation should fail with non-existent VCF file"
    assert "not found" in error

    # Test validation with non-existent reference file
    result, error = db.validate_vcf_reference(vcf_file, Path("nonexistent.fasta"))
    assert not result, "Validation should fail with non-existent reference file"
    assert "not found" in error
