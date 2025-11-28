"""Utility module for handling project paths and environment variables.

This module provides functions to get the project root directory and resource paths
regardless of how the package is installed.
"""

import os
from importlib import resources
from pathlib import Path


def get_project_root():
    """Get the project root directory regardless of how the package is installed."""
    # Always use VCFSTASH_ROOT if set (for Docker and custom setups)
    if "VCFSTASH_ROOT" in os.environ:
        return Path(os.environ["VCFSTASH_ROOT"])

    # Otherwise, use development directory structure (project root)
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
