"""Shared pytest fixtures for VCFstash tests."""

import os
import tempfile
import pytest
from pathlib import Path
from vcfstash.utils.paths import get_vcfstash_root

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