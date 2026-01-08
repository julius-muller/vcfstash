"""Tests for cache validation (all scenarios).

These tests validate cache structure when available. They adapt based on scenario:
- vanilla: Skip cache tests (no pre-built cache)
- blueprint: Validate /cache structure
- annotated: Validate /cache structure
"""

import pytest
import subprocess
import sys
from pathlib import Path
from tests.conftest import get_bcftools_cmd

# Command helpers
VCFCACHE_CMD = [sys.executable, "-m", "vcfcache"]


def get_bcftools():
    """Get bcftools command."""
    return get_bcftools_cmd()


def _skip_if_vanilla(test_scenario):
    """Helper to skip cache tests in vanilla scenario."""
    if test_scenario == "vanilla":
        pytest.skip("Vanilla scenario has no pre-built cache")


def test_cache_directory_exists(test_scenario):
    """Test that the cache directory exists (blueprint/annotated only)."""
    print(f"\n=== Testing cache directory (scenario: {test_scenario}) ===")

    if test_scenario == "vanilla":
        pytest.skip("Vanilla scenario has no pre-built cache")

    cache_dir = Path("/cache")
    assert cache_dir.exists(), "Cache directory /cache does not exist"
    assert cache_dir.is_dir(), "Cache path /cache is not a directory"



def test_cache_structure(test_scenario):
    """Test that the cache has the expected directory structure."""
    _skip_if_vanilla(test_scenario)

    cache_dir = Path("/cache")

    # Check for db directory
    db_dir = cache_dir / "db"
    assert db_dir.exists(), "DB directory does not exist"
    assert db_dir.is_dir(), "DB path is not a directory"

    # Check for blueprint directory (inside db)
    blueprint_dir = db_dir / "blueprint"
    assert blueprint_dir.exists(), "Blueprint directory does not exist"
    assert blueprint_dir.is_dir(), "Blueprint path is not a directory"



def test_blueprint_cache_file(test_scenario):
    """Test that the blueprint cache BCF file exists and is valid."""
    _skip_if_vanilla(test_scenario)

    cache_bcf = Path("/cache/db/blueprint/vcfcache.bcf")
    assert cache_bcf.exists(), f"Blueprint cache file not found: {cache_bcf}"

    # Check index exists
    cache_csi = Path("/cache/db/blueprint/vcfcache.bcf.csi")
    assert cache_csi.exists(), f"Blueprint cache index not found: {cache_csi}"

    # Verify it's a valid BCF file using bcftools
    result = subprocess.run(
        [get_bcftools(), "view", "-h", str(cache_bcf)],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, f"Blueprint cache is not a valid BCF file: {result.stderr}"
    assert "##fileformat=VCF" in result.stdout, "Blueprint cache missing VCF header"



def test_blueprint_has_variants(test_scenario):
    """Test that the blueprint cache contains variants."""
    _skip_if_vanilla(test_scenario)

    cache_bcf = Path("/cache/db/blueprint/vcfcache.bcf")

    # Get variant count
    result = subprocess.run(
        [get_bcftools(), "view", "-H", str(cache_bcf)],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, f"Failed to read blueprint cache: {result.stderr}"

    # Count variant lines (non-header, non-empty)
    variant_lines = [line for line in result.stdout.strip().split("\n") if line and not line.startswith("#")]
    variant_count = len(variant_lines)

    print(f"Blueprint cache contains {variant_count} variants")
    assert variant_count > 0, "Blueprint cache is empty (no variants)"



def test_cache_metadata(test_scenario):
    """Test that cache blueprints files exist."""
    _skip_if_vanilla(test_scenario)

    # Check for sources.info
    sources_info = Path("/cache/db/blueprint/sources.info")
    assert sources_info.exists(), "sources.info file not found"

    # Verify it contains expected blueprints
    with open(sources_info, 'r') as f:
        content = f.read()
        # Should contain some blueprints about the source BCF files
        assert len(content) > 0, "sources.info is empty"



def test_workflow_directory(test_scenario):
    """Test that workflow directory exists."""
    _skip_if_vanilla(test_scenario)

    workflow_dir = Path("/cache/db/workflow")
    assert workflow_dir.exists(), "Workflow directory does not exist"

    # The workflow directory exists for storing config snapshots

    # Check for init.yaml (configuration snapshot)
    init_yaml = workflow_dir / "init.yaml"
    assert init_yaml.exists(), "init.yaml not found in workflow directory"



def test_bcftools_version(test_scenario):
    """Test that bcftools is available and reports version."""
    result = subprocess.run(
        [get_bcftools(), "--version"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, "bcftools not available or failed"
    assert "bcftools" in result.stdout.lower(), "bcftools version output unexpected"

    # Extract version
    version_line = result.stdout.split("\n")[0]
    print(f"bcftools version: {version_line}")



def test_vcfcache_cli_available(test_scenario):
    """Test that vcfcache CLI is available."""
    result = subprocess.run(
        VCFCACHE_CMD + [ "--version"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, "vcfcache CLI not available"

    # Verify version output (check for semantic version format)
    version_output = result.stdout + result.stderr
    import re
    assert re.search(r'\d+\.\d+\.\d+', version_output), f"No version found in: {version_output}"



def test_cache_query_performance(test_scenario):
    """Test that we can query the cache efficiently."""
    _skip_if_vanilla(test_scenario)

    cache_bcf = Path("/cache/db/blueprint/vcfcache.bcf")

    # Query a specific region (should be fast)
    result = subprocess.run(
        [get_bcftools(), "view", "-H", "-r", "chr1:1-100000", str(cache_bcf)],
        capture_output=True,
        text=True,
        timeout=10  # Should complete in under 10 seconds
    )
    assert result.returncode == 0, f"Failed to query cache: {result.stderr}"

    print(f"Query completed successfully")



def test_cache_contig_format(test_scenario):
    """Test that cache uses correct chromosome naming (chr prefix)."""
    _skip_if_vanilla(test_scenario)

    cache_bcf = Path("/cache/db/blueprint/vcfcache.bcf")

    # Get header
    result = subprocess.run(
        [get_bcftools(), "view", "-h", str(cache_bcf)],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0

    # Check for chr-prefixed contigs
    contig_lines = [line for line in result.stdout.split("\n") if line.startswith("##contig=")]
    assert len(contig_lines) > 0, "No contig lines found in cache header"

    # Verify at least some contigs have chr prefix
    chr_contigs = [line for line in contig_lines if "ID=chr" in line]
    assert len(chr_contigs) > 0, "Cache does not use chr-prefixed chromosome names"

    print(f"Found {len(chr_contigs)} chromosomes with chr prefix")



def test_python_environment(test_scenario):
    """Test that Python environment is set up correctly in Docker image."""
    # Test that PYTHONPATH includes venv packages
    import sys

    # Should have access to vcfcache modules
    import vcfcache
    import vcfcache.cli
    import vcfcache.database

    # Verify installation
    from vcfcache.utils.validation import MIN_BCFTOOLS_VERSION
    assert MIN_BCFTOOLS_VERSION is not None

    print(f"VCFcache minimum required bcftools version: {MIN_BCFTOOLS_VERSION}")
    print(f"Python path: {sys.path}")
