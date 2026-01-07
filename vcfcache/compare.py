"""Compare two vcfcache annotate runs and report timing differences.

This module provides functionality to compare two successful vcfcache annotate
runs (typically cached vs uncached) and display performance metrics.
"""

import hashlib
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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


def parse_workflow_log(output_dir: Path) -> Tuple[Optional[float], List[Dict[str, str]]]:
    """Parse workflow.log to extract total time and detailed step timings.

    Args:
        output_dir: Output directory from vcfcache annotate

    Returns:
        Tuple of (total_time, step_timings)
        - total_time: Total workflow duration in seconds, or None if not found
        - step_timings: List of dicts with 'step', 'command', 'duration' keys
    """
    workflow_log = output_dir / "workflow.log"
    if not workflow_log.exists():
        return None, []

    total_time = None
    step_timings = []

    try:
        with open(workflow_log, "r") as f:
            for line in f:
                # Look for overall workflow completion: "Workflow completed successfully in 4592.2s"
                if "Workflow completed successfully in" in line:
                    match = re.search(r"completed successfully in ([\d.]+)s", line)
                    if match:
                        total_time = float(match.group(1))

                # Look for individual command timings: "Command completed in 32.733s: bcftools norm"
                elif "Command completed in" in line:
                    match = re.search(r"Command completed in ([\d.]+)s: (.+)$", line)
                    if match:
                        duration = float(match.group(1))
                        command = match.group(2).strip()
                        step_timings.append({
                            "duration": duration,
                            "command": command,
                        })

                # Also capture step descriptions for context
                elif "Step " in line and "/4:" in line:
                    # Extract step descriptions like "Step 1/4: Adding cache annotations"
                    match = re.search(r"Step (\d+/\d+): (.+)$", line)
                    if match:
                        step_num = match.group(1)
                        description = match.group(2).strip()
                        # Store this for the next timing entry
                        if step_timings and "step" not in step_timings[-1]:
                            step_timings[-1]["step"] = f"Step {step_num}: {description}"

    except Exception:
        pass

    return total_time, step_timings


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


def format_time(seconds: float) -> str:
    """Format seconds into a human-readable string.

    Args:
        seconds: Duration in seconds

    Returns:
        Formatted string like "1h 23m 45.6s" or "12m 34.5s" or "45.6s"
    """
    if seconds >= 3600:
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        secs = seconds % 60
        return f"{hours}h {minutes}m {secs:.1f}s"
    elif seconds >= 60:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.1f}s"
    else:
        return f"{seconds:.1f}s"


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
            f"\nNote: Both runs are in '{mode1}' mode. Comparison may not be meaningful.\n"
        )
        uncached_dir, cached_dir = dir1, dir2
        uncached_data, cached_data = data1, data2

    # Parse workflow logs
    uncached_time, uncached_steps = parse_workflow_log(uncached_dir)
    cached_time, cached_steps = parse_workflow_log(cached_dir)

    if uncached_time is None:
        raise ValueError(
            f"Could not extract timing from {uncached_dir}. "
            f"Check that workflow.log exists and contains 'Workflow completed successfully in' line."
        )
    if cached_time is None:
        raise ValueError(
            f"Could not extract timing from {cached_dir}. "
            f"Check that workflow.log exists and contains 'Workflow completed successfully in' line."
        )

    # Find output BCF files
    uncached_bcf = find_output_bcf(uncached_dir)
    cached_bcf = find_output_bcf(cached_dir)

    # Calculate MD5s if files exist
    uncached_md5 = get_md5(uncached_bcf) if uncached_bcf and uncached_bcf.exists() else "N/A"
    cached_md5 = get_md5(cached_bcf) if cached_bcf and cached_bcf.exists() else "N/A"

    # Display results
    print("\n" + "=" * 80)
    print("  VCFcache Annotation Comparison")
    print("=" * 80)
    print()
    print(f"Uncached run: {uncached_dir}")
    print(f"  Mode: {uncached_data.get('mode', 'unknown')}")
    print(f"  Version: {uncached_data.get('version', 'unknown')}")
    print(f"  Timestamp: {uncached_data.get('timestamp', 'unknown')}")
    print()
    print(f"Cached run: {cached_dir}")
    print(f"  Mode: {cached_data.get('mode', 'unknown')}")
    print(f"  Version: {cached_data.get('version', 'unknown')}")
    print(f"  Timestamp: {cached_data.get('timestamp', 'unknown')}")
    print()

    # Display detailed step timings
    if cached_steps:
        print("-" * 80)
        print("  Detailed Step Timings (Cached Run)")
        print("-" * 80)
        print()
        for i, step in enumerate(cached_steps, 1):
            duration = step.get("duration", 0)
            command = step.get("command", "unknown")
            print(f"  Step {i}: {command}")
            print(f"    Duration: {format_time(duration)} ({duration:.2f}s)")
        print()

    if uncached_steps:
        print("-" * 80)
        print("  Detailed Step Timings (Uncached Run)")
        print("-" * 80)
        print()
        for i, step in enumerate(uncached_steps, 1):
            duration = step.get("duration", 0)
            command = step.get("command", "unknown")
            print(f"  Step {i}: {command}")
            print(f"    Duration: {format_time(duration)} ({duration:.2f}s)")
        print()

    # Summary comparison
    print("=" * 80)
    print("  Summary")
    print("=" * 80)
    print()
    print(f"  {'Total Time (Uncached):':<35} {format_time(uncached_time):>20} ({uncached_time:,.2f}s)")
    print(f"  {'Total Time (Cached):':<35} {format_time(cached_time):>20} ({cached_time:,.2f}s)")
    print()
    speedup = uncached_time / cached_time if cached_time > 0 else 0
    time_saved = uncached_time - cached_time
    print(f"  {'Speed-up:':<35} {speedup:>20.2f}x")
    print(f"  {'Time Saved:':<35} {format_time(time_saved):>20} ({time_saved:,.2f}s)")
    if time_saved >= 60:
        print(f"  {'':35} ({time_saved / 60:>19.1f} minutes)")
    if time_saved >= 3600:
        print(f"  {'':35} ({time_saved / 3600:>19.2f} hours)")
    print()

    # Output verification
    print("-" * 80)
    print("  Output Verification (MD5)")
    print("-" * 80)
    print()
    print(f"  Uncached: {uncached_md5}")
    print(f"  Cached:   {cached_md5}")
    print()
    if uncached_md5 != "N/A" and cached_md5 != "N/A":
        if uncached_md5 == cached_md5:
            print("  ✓ Output files are IDENTICAL")
        else:
            print("  ✗ WARNING: Output files DIFFER!")
    else:
        print("  ⚠ Could not verify outputs (BCF files not found)")
    print()
    print("=" * 80)
    print()
