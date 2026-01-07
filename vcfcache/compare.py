"""Compare two vcfcache annotate runs and report timing differences.

This module provides functionality to compare two successful vcfcache annotate
runs (typically cached vs uncached) and display performance metrics.
"""

import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml

from vcfcache.utils.completion import read_completion_flag


def read_compare_stats(stats_dir: Path) -> Dict[str, any]:
    """Read compare_stats.yaml from stats directory."""
    stats_file = stats_dir / "compare_stats.yaml"
    if not stats_file.exists():
        return {}
    try:
        with open(stats_file, "r") as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}


def parse_workflow_log(output_dir: Path) -> Tuple[Optional[float], List[Dict[str, str]]]:
    """Parse workflow.log to extract total time and detailed step timings from stats dir.

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
    step_timings: List[Dict[str, str]] = []
    current_step = None

    try:
        with open(workflow_log, "r") as f:
            for line in f:
                # Look for overall workflow completion: "Workflow completed successfully in 4592.2s"
                if "Workflow completed successfully in" in line:
                    match = re.search(r"completed successfully in ([\d.]+)s", line)
                    if match:
                        total_time = float(match.group(1))

                # Capture step descriptions like "Step 1/4: Adding cache annotations"
                elif "Step " in line and "/4:" in line:
                    match = re.search(r"Step (\d+/\d+): (.+)$", line)
                    if match:
                        step_num = match.group(1)
                        description = match.group(2).strip()
                        current_step = f"Step {step_num}: {description}"

                # Look for individual command timings: "Command completed in 32.733s: bcftools norm"
                elif "Command completed in" in line:
                    match = re.search(r"Command completed in ([\d.]+)s: (.+)$", line)
                    if match:
                        duration = float(match.group(1))
                        command = match.group(2).strip()
                        step_timings.append({
                            "duration": duration,
                            "command": command,
                            "step": current_step,
                        })

    except Exception:
        pass

    return total_time, step_timings


def find_output_bcf(output_dir: Path) -> Optional[Path]:
    """Find the output BCF file from the annotate stats directory.

    Args:
        output_dir: Output directory from vcfcache annotate

    Returns:
        Path to output BCF file, or None if not found
    """
    completion = read_completion_flag(output_dir)
    if completion:
        output_file = completion.get("output_file")
        if output_file and output_file not in {"stdout", "-"}:
            output_path = Path(output_file).expanduser()
            if output_path.exists():
                return output_path

    # Fallback: look for *.bcf files in the stats directory
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
        dir1: First annotate stats directory
        dir2: Second annotate stats directory
    """
    # Validate both directories exist
    if not dir1.exists():
        raise FileNotFoundError(f"Directory not found: {dir1}")
    if not dir2.exists():
        raise FileNotFoundError(f"Directory not found: {dir2}")

    stats1 = read_compare_stats(dir1)
    stats2 = read_compare_stats(dir2)
    if not stats1 or not stats2:
        raise ValueError("Missing compare_stats.yaml in one or both stats directories.")

    # Compatibility checks
    input1 = stats1.get("input_name")
    input2 = stats2.get("input_name")
    if input1 and input2 and input1 != input2:
        raise ValueError(f"Input filename mismatch: {input1} vs {input2}")

    anno1 = stats1.get("annotation_yaml_md5")
    anno2 = stats2.get("annotation_yaml_md5")
    if anno1 and anno2 and anno1 != anno2:
        raise ValueError("annotation.yaml mismatch: the two runs used different annotation recipes.")

    warnings = []
    if stats1.get("vcfcache_version") != stats2.get("vcfcache_version"):
        warnings.append(
            f"WARNING: Different vcfcache versions ({stats1.get('vcfcache_version')} vs {stats2.get('vcfcache_version')})."
        )
    if stats1.get("genome_build_params") != stats2.get("genome_build_params"):
        warnings.append(
            f"WARNING: Different genome_build in params.yaml ({stats1.get('genome_build_params')} vs {stats2.get('genome_build_params')})."
        )
    if stats1.get("genome_build_annotation") != stats2.get("genome_build_annotation"):
        warnings.append(
            f"WARNING: Different genome_build in annotation.yaml ({stats1.get('genome_build_annotation')} vs {stats2.get('genome_build_annotation')})."
        )

    # Parse workflow logs
    time1, steps1 = parse_workflow_log(dir1)
    time2, steps2 = parse_workflow_log(dir2)

    # Determine comparator order (A = slower / longer runtime)
    if time1 is not None and time2 is not None and time1 < time2:
        dir_a, dir_b = dir2, dir1
        stats_a, stats_b = stats2, stats1
        time_a, time_b = time2, time1
        steps_a, steps_b = steps2, steps1
    else:
        dir_a, dir_b = dir1, dir2
        stats_a, stats_b = stats1, stats2
        time_a, time_b = time1, time2
        steps_a, steps_b = steps1, steps2

    def _counts(stats: Dict[str, any]) -> tuple[Optional[int], Optional[int], Optional[int], Optional[int]]:
        counts = stats.get("variant_counts", {}) or {}
        total = counts.get("total_output")
        annotated = counts.get("annotated_output")
        tool_annotated = counts.get("tool_annotated")
        dropped = counts.get("dropped_variants")
        mode = stats.get("mode")
        if annotated is None and mode == "uncached":
            annotated = total
        if tool_annotated is None:
            tool_annotated = annotated
        return total, annotated, tool_annotated, dropped

    total_a, annotated_a, tool_annotated_a, dropped_a = _counts(stats_a)
    total_b, annotated_b, tool_annotated_b, dropped_b = _counts(stats_b)

    def _rate(annotated: Optional[int], duration: Optional[float]) -> Optional[float]:
        if annotated is None or duration is None or duration <= 0:
            return None
        return annotated / duration

    rate_a = _rate(tool_annotated_a or annotated_a, time_a)
    rate_b = _rate(tool_annotated_b or annotated_b, time_b)

    md5_a = stats_a.get("variant_md5", {}) or {}
    md5_b = stats_b.get("variant_md5", {}) or {}

    md5_all_a = md5_a.get("all")
    md5_all_b = md5_b.get("all")

    print("\n" + "=" * 80)
    print("  VCFcache Run Comparison")
    print("=" * 80)
    print()
    print(f"Input file: {stats_a.get('input_name', 'unknown')}")
    print()
    if warnings:
        print("WARNINGS:")
        for w in warnings:
            print(f"  {w}")
        print()

    def _print_comparator(label: str, stats: Dict[str, any], total: Optional[int], annotated: Optional[int], tool_annotated: Optional[int], dropped: Optional[int], rate: Optional[float], total_time: Optional[float]) -> None:
        print(f"{label}: {stats.get('stats_dir', 'unknown')}")
        print(f"  Mode: {stats.get('mode', 'unknown')}")
        print(f"  Input file: {stats.get('input_name', 'unknown')}")
        print(f"  Cache name: {stats.get('cache_name', 'unknown')}")
        if stats.get("cache_path"):
            print(f"  Cache path: {stats.get('cache_path')}")
        print(f"  Version: {stats.get('vcfcache_version', 'unknown')}")
        print(f"  Run timestamp: {stats.get('run_timestamp', 'unknown')}")
        print(f"  Threads: {stats.get('threads', 'unknown')}")
        print(f"  Genome build (params.yaml): {stats.get('genome_build_params', 'N/A')}")
        print(f"  Genome build (annotation.yaml): {stats.get('genome_build_annotation', 'N/A')}")
        if total is not None:
            print(f"  Output variants: {total:,}")
        if tool_annotated is not None:
            print(f"  Annotated variants (tool): {tool_annotated:,}")
        if annotated is not None:
            print(f"  Annotated variants (output): {annotated:,}")
        if dropped is not None:
            print(f"  Dropped variants: {dropped:,}")
        if rate is not None:
            print(f"  Annotated variants/sec: {rate:,.2f}")
        if total_time is not None:
            print(f"  Total time: {format_time(total_time)} ({total_time:,.2f}s)")
        print()

    stats_a["stats_dir"] = str(dir_a)
    stats_b["stats_dir"] = str(dir_b)
    print("Comparator A (slower)")
    _print_comparator("  Stats dir", stats_a, total_a, annotated_a, tool_annotated_a, dropped_a, rate_a, time_a)
    print("Comparator B (faster)")
    _print_comparator("  Stats dir", stats_b, total_b, annotated_b, tool_annotated_b, dropped_b, rate_b, time_b)

    if time_a is not None and time_b is not None:
        speedup = time_a / time_b if time_b > 0 else 0
        time_saved = time_a - time_b
        print("Summary:")
        print(f"  Speed-up (A/B): {speedup:.2f}x")
        print(f"  Time saved: {format_time(time_saved)} ({time_saved:,.2f}s)")
        print()

    print("Output verification:")
    print(f"  Top10 MD5 A: {md5_a.get('top10') or 'N/A'}")
    print(f"  Top10 MD5 B: {md5_b.get('top10') or 'N/A'}")
    print(f"  Bottom10 MD5 A: {md5_a.get('bottom10') or 'N/A'}")
    print(f"  Bottom10 MD5 B: {md5_b.get('bottom10') or 'N/A'}")
    print(f"  Total MD5 A (all variants): {md5_all_a or 'N/A'}")
    print(f"  Total MD5 B (all variants): {md5_all_b or 'N/A'}")
    print()
    print("Detailed Step Timings:")
    def _print_steps(label: str, steps: List[Dict[str, str]]) -> None:
        print(f"  {label}:")
        if not steps:
            print("    (no detailed timings found)")
            return
        for entry in steps:
            step = entry.get("step")
            if step:
                print(f"    {step}")
            if entry.get("command") and entry.get("duration") is not None:
                print(f"      {format_time(entry['duration'])}  {entry['command']}")

    _print_steps("Comparator A", steps_a)
    _print_steps("Comparator B", steps_b)
    print()
    print("=" * 80)
    print()
