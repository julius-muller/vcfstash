# tests/conftest.py
import os
import pytest
from pathlib import Path


@pytest.fixture(scope="session", autouse=True)
def set_vcfstash_root():
    """Set the VCFSTASH_ROOT environment variable for all tests."""
    if 'VCFSTASH_ROOT' not in os.environ:
        package_root = Path(__file__).parent.parent.absolute()
        os.environ['VCFSTASH_ROOT'] = str(package_root)
        print(f"Set VCFSTASH_ROOT to {package_root}")
    return os.environ['VCFSTASH_ROOT']