"""Test core functionality of VCFcache."""

from pathlib import Path
import subprocess
import pytest
from vcfcache.utils.paths import get_vcfcache_root
from vcfcache.utils.validation import compute_md5

# Constants
TEST_ROOT = get_vcfcache_root() / "tests"
TEST_DATA_DIR = TEST_ROOT / "data" / "nodata"
TEST_VCF = TEST_DATA_DIR / "crayz_db.bcf"
TEST_VCF2 = TEST_DATA_DIR / "crayz_db2.bcf"
import sys

# Use python -m vcfcache to ensure we use the installed package
VCFCACHE_CMD = [sys.executable, "-m", "vcfcache"]


def test_error_handling(test_output_dir, params_file, test_scenario):
    """Test error conditions and edge cases."""
    print(f"\n=== Testing error handling (scenario: {test_scenario}) ===")

    # Test with non-existent input file
    init_cmd = VCFCACHE_CMD + [
        "blueprint-init",
        "-i",
        "nonexistent.bcf",
        "-o",
        test_output_dir,
        "-y",
        params_file
    ]

    result = subprocess.run(init_cmd, capture_output=True, text=True)
    assert result.returncode != 0, "Should fail with non-existent input"

    # Test with invalid output location
    init_cmd = VCFCACHE_CMD + [
        "blueprint-init",
        "-i",
        str(TEST_VCF),
        "-o",
        test_output_dir,
        "-y",
        "nonexistent.yaml",
    ]
    result = subprocess.run(init_cmd, capture_output=True, text=True)
    assert result.returncode != 0, "Should fail with invalid yaml"

    # Test add without init
    add_cmd = VCFCACHE_CMD + [ "blueprint-extend", "--db", test_output_dir, "-i", str(TEST_VCF)]
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

    from vcfcache.database.base import VCFDatabase
    from vcfcache.utils.validation import check_bcftools_installed
    import logging

    # Set up test files
    vcf_file = TEST_DATA_DIR / "crayz_db.bcf"
    ref_file = TEST_ROOT / "data/references/reference.fasta"

    # Get system bcftools path
    bcftools_path = check_bcftools_installed()

    # Create a VCFDatabase instance
    db = VCFDatabase(Path(TEST_ROOT), 2, True, Path(bcftools_path))
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


def test_show_command_outputs_annotation_tool_cmd(test_output_dir):
    """--show-command should print the frozen annotation tool command."""

    annotation_dir = Path(test_output_dir) / "cache" / "vep_gnomad"
    annotation_dir.mkdir(parents=True, exist_ok=True)

    expected_cmd = "docker run --rm -i vep:latest"
    (annotation_dir / "annotation.yaml").write_text(
        f'annotation_tool_cmd: "{expected_cmd}"\n'
    )

    result = subprocess.run(
        VCFCACHE_CMD + [
            "annotate",
            "--show-command",
            "-a",
            str(annotation_dir),
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert expected_cmd in result.stdout


def test_list_shows_available_annotations(test_output_dir):
    """--list should enumerate annotation directories containing vcfcache_annotated.bcf."""

    cache_dir = Path(test_output_dir) / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    anno1 = cache_dir / "vep_gnomad"
    anno1.mkdir()
    (anno1 / "vcfcache_annotated.bcf").write_text("dummy")

    anno2 = cache_dir / "custom"
    anno2.mkdir()
    (anno2 / "vcfcache_annotated.bcf").write_text("dummy")

    result = subprocess.run(
        VCFCACHE_CMD + [
            "annotate",
            "--list",
            "-a",
            str(cache_dir),
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "vep_gnomad" in result.stdout
    assert "custom" in result.stdout
