"""Completion flag utilities for tracking successful vcfcache runs.

This module provides utilities for writing and reading completion blueprints
that signals successful completion of vcfcache operations (cache-build,
blueprint-init, blueprint-extend, annotate).
"""

import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import yaml


def get_git_commit_hash() -> Optional[str]:
    """Get the current git commit hash if running from a git repository.

    Returns:
        Git commit hash or None if not in a git repository
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=2,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
        return None


def write_completion_flag(
    output_dir: Path,
    command: str,
    mode: Optional[str] = None,
    version: Optional[str] = None,
    output_file: Optional[str] = None,
) -> None:
    """Write a completion flag indicating successful run completion.

    Args:
        output_dir: Output directory where the flag should be written
        command: Command that was run (cache-build, blueprint-init, blueprint-extend, annotate)
        mode: Optional mode (cached, uncached, annotate-nocache)
        version: Optional vcfcache version string
    """
    if version is None:
        try:
            from vcfcache import __version__
            version = __version__
        except ImportError:
            version = "unknown"

    completion_data = {
        "completed": True,
        "command": command,
        "mode": mode,
        "version": version,
        "git_commit": (
            get_git_commit_hash()
            or os.environ.get("VCFCACHE_GIT_COMMIT")
            or os.environ.get("GIT_COMMIT")
            or "unknown"
        ),
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    if output_file:
        completion_data["output_file"] = output_file

    flag_file = output_dir / ".vcfcache_complete"
    with open(flag_file, "w") as f:
        yaml.dump(completion_data, f, default_flow_style=False)


def read_completion_flag(output_dir: Path) -> Optional[dict]:
    """Read completion flag from an output directory.

    Args:
        output_dir: Output directory to check for completion flag

    Returns:
        Dictionary with completion blueprints, or None if flag doesn't exist
    """
    flag_file = output_dir / ".vcfcache_complete"
    if not flag_file.exists():
        return None

    try:
        with open(flag_file, "r") as f:
            data = yaml.safe_load(f)
        return data
    except Exception:
        return None


def is_run_complete(output_dir: Path) -> bool:
    """Check if a run completed successfully.

    Args:
        output_dir: Output directory to check

    Returns:
        True if run completed successfully, False otherwise
    """
    data = read_completion_flag(output_dir)
    return data is not None and data.get("completed", False)


def validate_compatibility(dir1: Path, dir2: Path) -> tuple[bool, Optional[str]]:
    """Validate that two runs are compatible for comparison.

    Args:
        dir1: First output directory
        dir2: Second output directory

    Returns:
        Tuple of (is_compatible, error_message)
    """
    data1 = read_completion_flag(dir1)
    data2 = read_completion_flag(dir2)

    if data1 is None:
        return False, f"No completion flag found in {dir1}"

    if data2 is None:
        return False, f"No completion flag found in {dir2}"

    if data1["command"] != data2["command"]:
        return False, f"Different commands: {data1['command']} vs {data2['command']}"

    # Warn about version mismatches but allow comparison
    if data1["version"] != data2["version"]:
        warning = (
            f"Warning: Different versions detected "
            f"({data1['version']} vs {data2['version']}). "
            f"Results may not be directly comparable."
        )
        return True, warning

    return True, None
