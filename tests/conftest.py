"""Shared pytest fixtures for VCFstash tests."""

import os
import tempfile
import pytest
import subprocess
from pathlib import Path
from vcfstash.utils.paths import get_vcfstash_root
import shutil

# Constants
TEST_ROOT = get_vcfstash_root() / "tests"
TEST_DATA_DIR = TEST_ROOT / "data" / "nodata"
TEST_PARAMS = TEST_ROOT / "config" / "test_params.yaml"
TEST_ANNO_CONFIG = TEST_ROOT / "config" / "test_annotation.config"


# ============================================================================
# Scenario Detection Fixtures
# ============================================================================

@pytest.fixture(scope="session")
def test_scenario():
    """Detect which scenario we're running in.

    Returns:
        str: One of 'vanilla', 'blueprint', or 'annotated'

    Detection logic:
    - 'annotated': /cache exists AND vep is available
    - 'blueprint': /cache exists but vep is NOT available
    - 'vanilla': /cache does not exist
    """
    cache_exists = Path("/cache").exists()
    vep_available = shutil.which("vep") is not None

    if cache_exists and vep_available:
        return "annotated"
    elif cache_exists:
        return "blueprint"
    else:
        return "vanilla"


@pytest.fixture(scope="session")
def prebuilt_cache(test_scenario):
    """Provide path to pre-built cache if available.

    Returns:
        Path or None: Path to /cache if in blueprint/annotated scenario, None otherwise
    """
    if test_scenario in ["blueprint", "annotated"]:
        return Path("/cache")
    return None


@pytest.fixture
def annotation_config(test_scenario):
    """Provide appropriate annotation config based on scenario.

    Returns:
        Path: Path to annotation config file
        - 'annotated': Could use VEP config (for now using mock)
        - 'blueprint': Mock annotation with bcftools
        - 'vanilla': Mock annotation with bcftools
    """
    # For now, all scenarios use mock annotation
    # TODO: Add VEP annotation config for annotated scenario
    return TEST_ANNO_CONFIG


@pytest.fixture
def use_prebuilt_cache(test_scenario):
    """Boolean flag indicating whether to use pre-built cache.

    Returns:
        bool: True if in blueprint/annotated scenario, False for vanilla
    """
    return test_scenario in ["blueprint", "annotated"]


@pytest.fixture
def mini_cache_dir(test_output_dir, test_scenario, prebuilt_cache):
    """Create a mini test cache for blueprint scenario.

    In blueprint scenario, we don't want to modify /cache, so we create
    a small test cache from the top 10 variants of the blueprint.

    Returns:
        Path: Directory containing mini test cache
    """
    if test_scenario != "blueprint" or not prebuilt_cache:
        # Not needed for vanilla/annotated
        return None

    mini_cache_path = Path(test_output_dir) / "mini_cache"
    mini_cache_path.mkdir(parents=True, exist_ok=True)

    # Extract top 10 variants from blueprint cache
    blueprint_bcf = prebuilt_cache / "blueprint" / "vcfstash.bcf"
    mini_bcf = mini_cache_path / "mini_test.bcf"

    # Use bcftools to extract first 10 variants (header + 10 variant lines)
    # First get header, then get first 10 variants
    subprocess.run(
        f"bcftools view -h {blueprint_bcf} > {mini_cache_path}/header.vcf && "
        f"bcftools view -H {blueprint_bcf} | head -10 >> {mini_cache_path}/header.vcf && "
        f"bcftools view -Ob -o {mini_bcf} {mini_cache_path}/header.vcf",
        shell=True,
        check=True,
        timeout=30
    )

    # Index the mini BCF
    subprocess.run(
        ["bcftools", "index", str(mini_bcf)],
        check=True,
        timeout=30
    )

    return mini_cache_path


# ============================================================================
# Standard Fixtures
# ============================================================================

@pytest.fixture
def params_file():
    """Creates a temporary params file with correct paths."""
    vcfstash_root = str(get_vcfstash_root())
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)

    with open(TEST_PARAMS, 'r') as f:
        content = f.read().replace('${VCFSTASH_ROOT}', vcfstash_root)
        temp_file.write(content)
    temp_file.close()

    yield temp_file.name
    os.unlink(temp_file.name)


# This hook is needed to properly track test results
@pytest.hookimpl(tryfirst=True, hookwrapper=True)
def pytest_runtest_makereport(item, call):
    # Execute all other hooks to obtain the report object
    outcome = yield
    rep = outcome.get_result()

    # Set a report attribute for each phase of a call
    setattr(item, f"rep_{rep.when}", rep)


# Create or update the test_output_dir fixture
@pytest.fixture
def test_output_dir(request):
    """Provides a path for a directory that doesn't exist yet.
    If the test fails, the directory is NOT removed and a big warning is printed.
    """
    temp_dir = Path(tempfile.mkdtemp(prefix="vcfstash_test_"))

    # Remove the directory immediately - vcfstash will create it
    temp_dir.rmdir()

    yield str(temp_dir)

    # After the test, check if it passed
    # We must check for 'rep_call' and check its 'passed' attribute
    # Tests have 3 phases: setup, call, teardown
    # We only care about the 'call' phase result
    if hasattr(request.node, "rep_call") and request.node.rep_call.passed:
        # Test passed, clean up the directory
        if temp_dir.exists():
            shutil.rmtree(temp_dir, ignore_errors=True)
    else:
        # Test failed or no rep_call attribute, keep the directory and print a warning
        if temp_dir.exists():
            banner = "\n" + "=" * 80
            print(
                f"{banner}\n"
                f"TEST FAILED! Temporary directory NOT removed for forensic analysis:\n"
                f"    {temp_dir}\n"
                f"Please clean up manually after investigation.\n"
                f"{banner}\n"
            )