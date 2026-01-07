"""Compare two vcfcache annotate runs and report timing differences.

This module provides functionality to compare two successful vcfcache annotate
runs (typically cached vs uncached) and display performance metrics.
"""

import os
import re
import subprocess
import sys
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


def parse_missing_variants(output_dir: Path) -> Optional[int]:
    """Parse workflow.log to extract missing variants count for cached runs."""
    workflow_log = output_dir / "workflow.log"
    if not workflow_log.exists():
        return None
    try:
        with open(workflow_log, "r") as f:
            for line in f:
                match = re.search(r"Annotating (\d+) missing variants", line)
                if match:
                    return int(match.group(1))
    except Exception:
        return None
    return None


def _load_params_snapshot(stats_dir: Path) -> Dict[str, any]:
    params_path = stats_dir / "workflow" / "params.snapshot.yaml"
    if not params_path.exists():
        return {}
    try:
        with open(params_path, "r") as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}


def _load_annotation_snapshot(stats_dir: Path) -> str:
    anno_path = stats_dir / "workflow" / "annotation.snapshot.yaml"
    if not anno_path.exists():
        return ""
    try:
        return anno_path.read_text(encoding="utf-8")
    except Exception:
        return ""


def _required_param_keys(annotation_text: str) -> set[str]:
    filtered_lines = []
    for line in annotation_text.splitlines():
        if line.lstrip().startswith("#"):
            continue
        filtered_lines.append(line)
    filtered_text = "\n".join(filtered_lines)
    keys = set(re.findall(r"\$\{params\.([A-Za-z0-9_\\.]+)\}", filtered_text))
    top_level = set()
    for key in keys:
        top_level.add(key.split(".")[0])
    base_required = {
        "genome_build",
        "bcftools_cmd",
        "annotation_tool_cmd",
        "tool_version_command",
        "temp_dir",
        "threads",
    }
    return base_required | top_level


def _extra_param_keys(params: Dict[str, any], annotation_text: str) -> list[str]:
    if not isinstance(params, dict):
        return []
    required = _required_param_keys(annotation_text)
    extra = [k for k in params.keys() if k not in required]
    return sorted(extra)


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
    missing1 = parse_missing_variants(dir1)
    missing2 = parse_missing_variants(dir2)

    # Fill threads/timestamps from snapshots/completion flag
    params1 = _load_params_snapshot(dir1)
    params2 = _load_params_snapshot(dir2)
    anno_text1 = _load_annotation_snapshot(dir1)
    anno_text2 = _load_annotation_snapshot(dir2)
    extra_params1 = _extra_param_keys(params1, anno_text1)
    extra_params2 = _extra_param_keys(params2, anno_text2)

    completion1 = read_completion_flag(dir1) or {}
    completion2 = read_completion_flag(dir2) or {}

    if not stats1.get("threads"):
        stats1["threads"] = params1.get("threads")
    if not stats2.get("threads"):
        stats2["threads"] = params2.get("threads")
    if not stats1.get("run_timestamp"):
        stats1["run_timestamp"] = completion1.get("timestamp")
    if not stats2.get("run_timestamp"):
        stats2["run_timestamp"] = completion2.get("timestamp")

    if stats1.get("threads") is None:
        raise ValueError(f"Missing threads in stats or params snapshot for {dir1}")
    if stats2.get("threads") is None:
        raise ValueError(f"Missing threads in stats or params snapshot for {dir2}")
    if not stats1.get("run_timestamp"):
        raise ValueError(f"Missing run timestamp in stats or completion flag for {dir1}")
    if not stats2.get("run_timestamp"):
        raise ValueError(f"Missing run timestamp in stats or completion flag for {dir2}")

    # Determine comparator order (A = slower / longer runtime)
    if time1 is not None and time2 is not None and time1 < time2:
        dir_a, dir_b = dir2, dir1
        stats_a, stats_b = stats2, stats1
        time_a, time_b = time2, time1
        steps_a, steps_b = steps2, steps1
        missing_a, missing_b = missing2, missing1
    else:
        dir_a, dir_b = dir1, dir2
        stats_a, stats_b = stats1, stats2
        time_a, time_b = time1, time2
        steps_a, steps_b = steps1, steps2
        missing_a, missing_b = missing1, missing2

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

    def _fallback_tool_annotated(
        stats: Dict[str, any],
        tool_annotated: Optional[int],
        missing: Optional[int],
        total: Optional[int],
    ) -> Optional[int]:
        mode = stats.get("mode")
        if tool_annotated is None and mode == "cached":
            return missing
        if mode == "cached" and missing is not None and tool_annotated == total:
            return missing
        return tool_annotated

    tool_annotated_a = _fallback_tool_annotated(stats_a, tool_annotated_a, missing_a, total_a)
    tool_annotated_b = _fallback_tool_annotated(stats_b, tool_annotated_b, missing_b, total_b)

    def _rate(count: Optional[int], duration: Optional[float]) -> Optional[float]:
        if count is None or duration is None or duration <= 0:
            return None
        return count / duration

    def _tool_step_time(steps: List[Dict[str, str]], mode: str) -> Optional[float]:
        if not steps:
            return None
        if mode == "cached":
            annotate_steps = [s for s in steps if s.get("step") and "Annotating" in s.get("step")]
            if annotate_steps:
                return sum(s.get("duration", 0.0) for s in annotate_steps if s.get("duration") is not None)
        view_steps = [s for s in steps if s.get("command", "").startswith("bcftools view")]
        if view_steps:
            return max(s.get("duration", 0.0) for s in view_steps if s.get("duration") is not None)
        return None

    tool_time_a = _tool_step_time(steps_a, stats_a.get("mode", ""))
    tool_time_b = _tool_step_time(steps_b, stats_b.get("mode", ""))

    rate_a = _rate(total_a, time_a)
    rate_b = _rate(total_b, time_b)
    tool_rate_a = _rate(tool_annotated_a or annotated_a, tool_time_a)
    tool_rate_b = _rate(tool_annotated_b or annotated_b, tool_time_b)

    md5_a = stats_a.get("variant_md5", {}) or {}
    md5_b = stats_b.get("variant_md5", {}) or {}

    md5_all_a = md5_a.get("all")
    md5_all_b = md5_b.get("all")

    use_color = sys.stdout.isatty() and os.environ.get("NO_COLOR") is None
    C_RESET = "\033[0m" if use_color else ""
    C_BOLD = "\033[1m" if use_color else ""
    C_CYAN = "\033[36m" if use_color else ""
    C_GREEN = "\033[32m" if use_color else ""
    C_YELLOW = "\033[33m" if use_color else ""

    def _hdr(text: str) -> str:
        return f"{C_BOLD}{C_CYAN}{text}{C_RESET}"

    def _ok(text: str) -> str:
        return f"{C_GREEN}{text}{C_RESET}"

    def _warn(text: str) -> str:
        return f"{C_YELLOW}{text}{C_RESET}"

    print("\n" + "=" * 80)
    print(_hdr("  VCFcache Run Comparison"))
    print("=" * 80)
    print()
    print(_hdr("Samples"))
    print(f"  Input file: {stats_a.get('input_name', 'unknown')}")
    print(f"  Comparator A stats dir: {dir_a}")
    print(f"  Comparator B stats dir: {dir_b}")
    print(f"  Cache A: {stats_a.get('cache_name', 'unknown')}")
    print(f"  Cache B: {stats_b.get('cache_name', 'unknown')}")
    print()
    if warnings:
        print(_hdr("Warnings"))
        for w in warnings:
            print(f"  {_warn(w)}")
        print()

    def _print_comparator(
        label: str,
        stats: Dict[str, any],
        total: Optional[int],
        annotated: Optional[int],
        tool_annotated: Optional[int],
        dropped: Optional[int],
        total_rate: Optional[float],
        total_time: Optional[float],
        tool_rate: Optional[float],
        tool_time: Optional[float],
        extra_params: list[str],
    ) -> None:
        print(label)
        items = [
            ("Mode", stats.get("mode", "unknown")),
            ("Cache name", stats.get("cache_name", "unknown")),
            ("Cache path", stats.get("cache_path")),
            ("Version", stats.get("vcfcache_version", "unknown")),
            ("Run timestamp", stats.get("run_timestamp", "unknown")),
            ("Threads", stats.get("threads", "unknown")),
            ("Genome build (params.yaml)", stats.get("genome_build_params", "N/A")),
            ("Genome build (annotation.yaml)", stats.get("genome_build_annotation", "N/A")),
            ("Extra params", ", ".join(extra_params) if extra_params else "(none)"),
            ("Output variants (total)", f"{total:,}" if total is not None else None),
            ("Annotated variants in output", f"{annotated:,}" if annotated is not None else None),
            ("Annotated variants (tool)", f"{tool_annotated:,}" if tool_annotated is not None else None),
            ("Dropped variants", f"{dropped:,}" if dropped is not None else None),
            ("Output variants/sec (end-to-end)", f"{total_rate:,.2f}" if total_rate is not None else None),
            ("Annotated variants/sec (tool step)", f"{tool_rate:,.2f}" if tool_rate is not None else None),
            ("Total time", f"{format_time(total_time)} ({total_time:,.2f}s)" if total_time is not None else None),
            ("Tool time", f"{format_time(tool_time)} ({tool_time:,.2f}s)" if tool_time is not None else None),
        ]
        key_width = max(len(k) for k, _ in items)
        for key, value in items:
            if value is None:
                continue
            print(f"  {key.ljust(key_width)} : {value}")
        print()

    print(_hdr("Comparator A (slower)"))
    _print_comparator(
        "Details",
        stats_a,
        total_a,
        annotated_a,
        tool_annotated_a,
        dropped_a,
        rate_a,
        time_a,
        tool_rate_a,
        tool_time_a,
        extra_params1 if stats_a is stats1 else extra_params2,
    )
    print(_hdr("Comparator B (faster)"))
    _print_comparator(
        "Details",
        stats_b,
        total_b,
        annotated_b,
        tool_annotated_b,
        dropped_b,
        rate_b,
        time_b,
        tool_rate_b,
        tool_time_b,
        extra_params2 if stats_b is stats2 else extra_params1,
    )

    print("Note: Output variants (total) can include unannotated records (e.g., --preserve-unannotated).")
    print()
    print(_hdr("Detailed Step Timings"))
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
                duration = format_time(entry["duration"])
                print(f"      {duration.rjust(8)}  {entry['command']}")

    _print_steps("Comparator A", steps_a)
    _print_steps("Comparator B", steps_b)
    print()

    print(_hdr("Summary"))
    speedup = None
    time_saved = None
    if time_a is not None and time_b is not None and time_b > 0:
        speedup = time_a / time_b
        time_saved = time_a - time_b

    print("+---------------------------+----------------------+----------------------+")
    print("| Metric                    | Comparator A         | Comparator B         |")
    print("+---------------------------+----------------------+----------------------+")
    print(f"| End-to-end time           | {format_time(time_a) if time_a is not None else 'N/A':<20} | {format_time(time_b) if time_b is not None else 'N/A':<20} |")
    print(f"| Tool time                 | {format_time(tool_time_a) if tool_time_a is not None else 'N/A':<20} | {format_time(tool_time_b) if tool_time_b is not None else 'N/A':<20} |")
    print(f"| End-to-end rate (var/s)   | {f'{rate_a:,.2f}' if rate_a is not None else 'N/A':<20} | {f'{rate_b:,.2f}' if rate_b is not None else 'N/A':<20} |")
    print(f"| Tool rate (var/s)         | {f'{tool_rate_a:,.2f}' if tool_rate_a is not None else 'N/A':<20} | {f'{tool_rate_b:,.2f}' if tool_rate_b is not None else 'N/A':<20} |")
    print("+---------------------------+----------------------+----------------------+")
    if speedup is not None and time_saved is not None:
        print(f"  Speed-up (A/B): {speedup:.2f}x  (end-to-end)")
        print(f"  Time saved: {format_time(time_saved)} ({time_saved:,.2f}s)")
    else:
        print(f"  {_warn('[!]')} Speed-up: N/A (end-to-end timing missing)")

    def _md5_status(a: Optional[str], b: Optional[str]) -> str:
        if not a or not b:
            return _warn("[!]") + " N/A"
        if a == b:
            return _ok("[OK]") + " match"
        return _warn("[!]") + " differ"

    print("  Output verification:")
    print(f"    Top10 MD5: {md5_a.get('top10') or 'N/A'} vs {md5_b.get('top10') or 'N/A'}  {_md5_status(md5_a.get('top10'), md5_b.get('top10'))}")
    print(f"    Bottom10 MD5: {md5_a.get('bottom10') or 'N/A'} vs {md5_b.get('bottom10') or 'N/A'}  {_md5_status(md5_a.get('bottom10'), md5_b.get('bottom10'))}")
    print(f"    Total MD5 (all variants): {md5_all_a or 'N/A'} vs {md5_all_b or 'N/A'}  {_md5_status(md5_all_a, md5_all_b)}")
    print()
    print("=" * 80)
    print()
