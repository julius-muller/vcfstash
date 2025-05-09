"""Shared pytest fixtures for VCFstash tests."""

import os
import tempfile
import pytest
from pathlib import Path
from vcfstash.utils.paths import get_vcfstash_root
import shutil

# Constants
TEST_ROOT = get_vcfstash_root() / "tests"
TEST_PARAMS = TEST_ROOT / "config" / "test_params.yaml"

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