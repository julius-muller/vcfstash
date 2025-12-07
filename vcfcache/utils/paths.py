"""Utility module for handling project paths and environment variables.

This module provides functions to get the project root directory and resource paths
regardless of how the package is installed.
"""

import os
from importlib import resources
from pathlib import Path


def get_project_root():
    """Get the project root directory regardless of how the package is installed."""
    # Always use VCFCACHE_ROOT if set (for Docker and custom setups)
    if "VCFCACHE_ROOT" in os.environ:
        return Path(os.environ["VCFCACHE_ROOT"])

    # Otherwise, use development directory structure (project root)
    return Path(__file__).parent.parent.parent


# Set VCFCACHE_ROOT if not already set
if "VCFCACHE_ROOT" not in os.environ:
    os.environ["VCFCACHE_ROOT"] = str(get_project_root())


def get_vcfcache_root() -> Path:
    """Get the VCFCACHE_ROOT directory."""
    if "VCFCACHE_ROOT" not in os.environ:
        os.environ["VCFCACHE_ROOT"] = str(get_project_root())
    return Path(os.environ["VCFCACHE_ROOT"])


def get_resource_path(relative_path: Path) -> Path:
    """Get the absolute path to a resource file."""
    return get_vcfcache_root() / relative_path
