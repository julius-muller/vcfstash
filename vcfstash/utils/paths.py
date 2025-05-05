"""
Utility module for handling project paths and environment variables.

This module provides functions to get the project root directory and resource paths
regardless of how the package is installed.
"""

import os
import sys
from pathlib import Path
from importlib import resources


def get_project_root():
    """Get the project root directory regardless of how the package is installed."""
    try:
        # When installed as a package
        path = resources.files("vcfstash") / "__init__.py"
        return path.parent.parent
    except (ImportError, ModuleNotFoundError):
        # When running from source
        return Path(__file__).parent.parent.parent


# Set VCFSTASH_ROOT if not already set
if "VCFSTASH_ROOT" not in os.environ:
    os.environ["VCFSTASH_ROOT"] = str(get_project_root())


def get_vcfstash_root() -> Path:
    """Get the VCFSTASH_ROOT directory."""
    if "VCFSTASH_ROOT" not in os.environ:
        os.environ["VCFSTASH_ROOT"] = str(get_project_root())
    return Path(os.environ["VCFSTASH_ROOT"])


def get_resource_path(relative_path: Path) -> Path:
    """Get the absolute path to a resource file."""
    return get_vcfstash_root() / relative_path
