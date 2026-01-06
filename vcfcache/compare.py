"""Compare two vcfcache annotate runs and report timing differences.

This module provides functionality to compare two successful vcfcache annotate
runs (typically cached vs uncached) and display performance metrics.
"""

import hashlib
import subprocess
from pathlib import Path
from typing import Optional

from vcfcache.utils.completion import read_completion_flag, validate_compatibility


def get_md5(file_path: Path) -> str:
    """Calculate MD5 checksum of a file.

    Args:
        file_path: Path to the file

    Returns:
        MD5 checksum as hex string
    """
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def extract_timing(output_dir: Path) -> Optional[float]:
    """Extract timing information from a vcfcache annotate run.

    Args:
        output_dir: Output directory from vcfcache annotate

    Returns:
        Duration in seconds, or None if timing not found
    """
    # Check workflow.log for timing information
    workflow_log = output_dir / "workflow.log"
    if not workflow_log.exists():
        return None

    try:
        with open(workflow_log, "r") as f:
            for line in f:
                # Look for timing lines like "Command completed in 123.45s"
                if "completed in" in line and "s" in line:
                    # Extract the number before 's'
                    parts = line.split("completed in")
                    if len(parts) == 2:
                        time_str = parts[1].strip().rstrip("s")
                        try:
                            return float(time_str)
                        except ValueError:
                            continue
    except Exception:
        pass

    # Try timing.txt file in work directory
    timing_file = output_dir / "work" / "timing.txt"
    if timing_file.exists():
        try:
            content = timing_file.read_text().strip()
            return float(content)
        except (ValueError, IOError):
            pass

    return None


def find_output_bcf(output_dir: Path) -> Optional[Path]:
    """Find the output BCF file in the annotate output directory.

    Args:
        output_dir: Output directory from vcfcache annotate

    Returns:
        Path to output BCF file, or None if not found
    """
    # Look for *.bcf files in the output directory (not in work subdirectory)
    bcf_files = [f for f in output_dir.glob("*.bcf") if not f.name.startswith("work")]
    if bcf_files:
        return bcf_files[0]

    return None


def compare_runs(dir1: Path, dir2: Path) -> None:
    """Compare two vcfcache annotate runs and display results.

    Args:
        dir1: First annotate output directory
        dir2: Second annotate output directory
    """
    # Validate both directories exist
    if not dir1.exists():
        raise FileNotFoundError(f"Directory not found: {dir1}")
    if not dir2.exists():
        raise FileNotFoundError(f"Directory not found: {dir2}")

    # Read completion flags
    data1 = read_completion_flag(dir1)
    data2 = read_completion_flag(dir2)

    # Validate compatibility
    is_compatible, message = validate_compatibility(dir1, dir2)
    if not is_compatible:
        raise ValueError(f"Runs are not compatible: {message}")

    if message:  # Warning message
        print(f"\n{message}\n")

    # Determine which is cached and which is uncached
    mode1 = data1.get("mode", "unknown")
    mode2 = data2.get("mode", "unknown")

    if mode1 == "uncached" and mode2 == "cached":
        uncached_dir, cached_dir = dir1, dir2
        uncached_data, cached_data = data1, data2
    elif mode1 == "cached" and mode2 == "uncached":
        uncached_dir, cached_dir = dir2, dir1
        uncached_data, cached_data = data2, data1
    else:
        # Both are same mode or unknown, just compare them as run1 vs run2
        print(
            f"\nNote: Both runs are in {mode1} mode. Comparison may not be meaningful.\n"
        )
        uncached_dir, cached_dir = dir1, dir2
        uncached_data, cached_data = data1, data2

    # Extract timing
    uncached_time = extract_timing(uncached_dir)
    cached_time = extract_timing(cached_dir)

    if uncached_time is None:
        raise ValueError(f"Could not extract timing from {uncached_dir}")
    if cached_time is None:
        raise ValueError(f"Could not extract timing from {cached_dir}")

    # Find output BCF files
    uncached_bcf = find_output_bcf(uncached_dir)
    cached_bcf = find_output_bcf(cached_dir)

    # Calculate MD5s if files exist
    uncached_md5 = get_md5(uncached_bcf) if uncached_bcf and uncached_bcf.exists() else "N/A"
    cached_md5 = get_md5(cached_bcf) if cached_bcf and cached_bcf.exists() else "N/A"

    # Display results
    print("\n" + "=" * 70)
    print("  VCFcache Annotation Comparison")
    print("=" * 70)
    print()
    print(f"Run 1 ({uncached_data.get('mode', 'unknown')}): {uncached_dir}")
    print(f"Run 2 ({cached_data.get('mode', 'unknown')}): {cached_dir}")
    print()
    print("-" * 70)
    print("  Timing Results")
    print("-" * 70)
    print()
    print(f"{'Uncached annotation time:':<30} {uncached_time:>10.2f}s")
    print(f"{'Cached annotation time:':<30} {cached_time:>10.2f}s")
    print()
    print(f"{'Speed-up:':<30} {uncached_time / cached_time:>10.2f}x")
    print(f"{'Time saved:':<30} {uncached_time - cached_time:>10.2f}s ({(uncached_time - cached_time) / 60:.1f} minutes)")
    print()
    print("-" * 70)
    print("  Output Verification")
    print("-" * 70)
    print()
    print(f"Uncached output MD5: {uncached_md5}")
    print(f"Cached output MD5:   {cached_md5}")
    print()
    if uncached_md5 != "N/A" and cached_md5 != "N/A":
        if uncached_md5 == cached_md5:
            print("✓ Output files are identical")
        else:
            print("✗ WARNING: Output files differ!")
    print()
    print("=" * 70)
    print()
