"""Tests for vanilla vcfcache package (no cache, no external tools).

These tests should pass in a plain Python environment without any
cache files, Nextflow, or annotation tools.
"""

import pytest
import subprocess
import tempfile
from pathlib import Path


def test_module_imports(test_scenario):
    """Test that all core modules can be imported."""
    # Core modules
    import vcfcache
    import vcfcache.cli
    import vcfcache.database
    import vcfcache.utils

    # Utils
    import vcfcache.utils.paths
    import vcfcache.utils.validation
    import vcfcache.utils.logging

    # Database classes
    from vcfcache.database.base import VCFDatabase
    from vcfcache.database.initializer import DatabaseInitializer
    from vcfcache.database.updater import DatabaseUpdater
    from vcfcache.database.annotator import DatabaseAnnotator, VCFAnnotator

    # Verify constants are set
    assert hasattr(vcfcache, 'EXPECTED_BCFTOOLS_VERSION')
    assert isinstance(vcfcache.EXPECTED_BCFTOOLS_VERSION, str)


def test_cli_help(test_scenario):
    """Test that CLI help command works."""
    result = subprocess.run(
        ["vcfcache", "--help"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0
    # Check for key content in help text (case-insensitive)
    stdout_lower = result.stdout.lower()
    assert "vcf annotation" in stdout_lower or "blueprint-init" in stdout_lower
    assert "blueprint-init" in result.stdout
    assert "blueprint-extend" in result.stdout
    assert "cache-build" in result.stdout
    assert "annotate" in result.stdout


def test_cli_version(test_scenario):
    """Test that CLI version command works."""
    result = subprocess.run(
        ["vcfcache", "--version"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0
    # Version should be in stdout or stderr (check for semantic version format)
    version_output = result.stdout + result.stderr
    # Check for version format like "0.1.0" or "0.2.0"
    import re
    assert re.search(r'\d+\.\d+\.\d+', version_output), f"No version found in: {version_output}"


def test_blueprint_init_help(test_scenario):
    """Test blueprint-init help command."""
    result = subprocess.run(
        ["vcfcache", "blueprint-init", "--help"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0
    assert "blueprint-init" in result.stdout
    assert "--vcf" in result.stdout or "-i" in result.stdout
    assert "--output" in result.stdout or "-o" in result.stdout


def test_error_handling_missing_vcf(test_scenario):
    """Test error handling for missing VCF file."""
    result = subprocess.run(
        ["vcfcache", "blueprint-init",
         "--vcf", "nonexistent.bcf",
         "--output", "/tmp/test",
         "-y", "nonexistent.yaml"],
        capture_output=True,
        text=True
    )
    # Should fail with non-zero exit code
    assert result.returncode != 0


def test_error_handling_missing_params(test_scenario):
    """Test error handling for missing params file."""
    # Create a temporary VCF path (doesn't need to exist for this test)
    result = subprocess.run(
        ["vcfcache", "blueprint-init",
         "--vcf", "test.bcf",
         "--output", "/tmp/test",
         "-y", "definitely_nonexistent_file.yaml"],
        capture_output=True,
        text=True
    )
    # Should fail with non-zero exit code
    assert result.returncode != 0


def test_check_duplicate_md5(test_scenario):
    """Test duplicate MD5 checking logic."""
    from vcfcache.utils.validation import check_duplicate_md5

    # Test with empty info
    db_info = {"input_files": []}
    assert not check_duplicate_md5(db_info, "abc123")

    # Test with matching MD5
    db_info = {"input_files": [{"md5": "abc123", "file": "test.bcf"}]}
    assert check_duplicate_md5(db_info, "abc123")

    # Test with non-matching MD5
    assert not check_duplicate_md5(db_info, "xyz789")

    # Test with missing key
    db_info = {}
    assert not check_duplicate_md5(db_info, "abc123")


def test_md5_computation(test_scenario):
    """Test MD5 computation with a temporary file."""
    # Create a temporary file with known content
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        f.write("test content\n")
        temp_file = Path(f.name)

    try:
        from vcfcache.utils.validation import compute_md5

        # Compute MD5
        md5_hash = compute_md5(temp_file)

        # Verify it's a valid MD5 hash (32 hex characters)
        assert len(md5_hash) == 32
        assert all(c in '0123456789abcdef' for c in md5_hash)

        # Compute again to verify consistency
        md5_hash2 = compute_md5(temp_file)
        assert md5_hash == md5_hash2

    finally:
        temp_file.unlink()


def test_path_resolution(test_scenario):
    """Test VCFCACHE_ROOT path resolution."""
    from vcfcache.utils.paths import get_vcfcache_root, get_resource_path

    # Get root path
    root = get_vcfcache_root()
    assert root.exists()
    assert root.is_dir()

    # Verify it contains expected directories based on environment
    # In development: root/vcfcache, root/resources, root/tools
    # In Docker: root/resources, root/tools (vcfcache is in venv)
    has_dev_structure = (root / "vcfcache").exists()
    has_docker_structure = (root / "resources").exists() and (root / "tools").exists()

    assert has_dev_structure or has_docker_structure, \
        f"Expected either dev structure (vcfcache/) or Docker structure (resources/, tools/) in {root}"

    # Test resource path resolution
    resource_path = get_resource_path(Path("resources/chr_add.txt"))
    # Path should be constructed, whether or not it exists
    assert isinstance(resource_path, Path)


def test_bcftools_expected_version(test_scenario):
    """Test that expected bcftools version is defined."""
    from vcfcache import EXPECTED_BCFTOOLS_VERSION

    assert isinstance(EXPECTED_BCFTOOLS_VERSION, str)
    assert len(EXPECTED_BCFTOOLS_VERSION) > 0
    # Should be in format like "1.20" or "1.20.1"
    assert EXPECTED_BCFTOOLS_VERSION[0].isdigit()
