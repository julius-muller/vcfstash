# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

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

    here = Path(__file__).resolve()
    dev_root = here.parent.parent.parent  # repo root in editable/dev installs
    pkg_root = here.parent.parent         # installed package root (site-packages/vcfcache)

    # Prefer repo root if we see pyproject.toml (editable/development)
    if (dev_root / "pyproject.toml").exists():
        return dev_root

    # Fallback to package root when installed from wheel/sdist
    return pkg_root


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
