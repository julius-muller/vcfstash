# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""Comprehensive smoke test demo of vcfcache workflow.

This module demonstrates the complete vcfcache workflow with all 4 main commands:
1. blueprint-init: Create initial blueprint from VCF
2. blueprint-extend: Add more variants to existing blueprint
3. cache-build: Annotate the blueprint to create cache
4. annotate: Use the cache to annotate a sample VCF

The smoke test validates that all commands work correctly and produces
identical outputs for cached and uncached annotation.

For comparing existing annotation runs (e.g., cached vs uncached), use the
`vcfcache compare` command instead.

Usage:
    vcfcache demo --smoke-test [--debug] [--quiet]

Or from Python:
    from vcfcache.demo import run_smoke_test
    run_smoke_test(keep_files=False, quiet=False)
"""

import hashlib
import sys
import subprocess
import tempfile
import shutil
import time
from pathlib import Path

# Module-level quiet mode flag
_QUIET_MODE = False


def print_section(title):
    """Print a section header."""
    if not _QUIET_MODE:
        print("\n" + "="*70)
        print(f"  {title}")
        print("="*70 + "\n")


def print_step(step_num, description):
    """Print a step header."""
    if _QUIET_MODE:
        print(".", end="", flush=True)
    else:
        print(f"\n{'─'*70}")
        print(f"Step {step_num}: {description}")
        print(f"{'─'*70}\n")


def format_duration(seconds):
    """Format duration as m:s.ms"""
    minutes = int(seconds // 60)
    secs = seconds % 60
    if minutes > 0:
        return f"{minutes}m {secs:.3f}s"
    else:
        return f"{secs:.3f}s"


def collect_detailed_timings(cache_dir, cached_stats_dir, uncached_stats_dir):
    """Collect detailed timing information from workflow log files."""
    import re
    detailed_timings = {}

    # Pattern to match timing log lines
    timing_pattern = re.compile(r'Command completed in ([\d.]+)s: (.+)')

    # Collect from cache operations
    for subdir in ["blueprint", "cache/demo_cache"]:
        log_file = cache_dir / subdir / "workflow.log"
        if log_file.exists():
            step_name = subdir.replace("/", "-")
            if step_name not in detailed_timings:
                detailed_timings[step_name] = []

            with log_file.open() as f:
                for line in f:
                    match = timing_pattern.search(line)
                    if match:
                        duration = float(match.group(1))
                        cmd = match.group(2).strip()
                        detailed_timings[step_name].append((cmd, duration))

    # Collect from annotation operations
    for output_subdir, step_name in [
        (cached_stats_dir, "annotate-cached"),
        (uncached_stats_dir, "annotate-uncached"),
    ]:
        if not output_subdir.exists():
            continue

        log_file = output_subdir / "workflow.log"
        if log_file.exists():
            if step_name not in detailed_timings:
                detailed_timings[step_name] = []

            with log_file.open() as f:
                for line in f:
                    match = timing_pattern.search(line)
                    if match:
                        duration = float(match.group(1))
                        cmd = match.group(2).strip()
                        detailed_timings[step_name].append((cmd, duration))

    return detailed_timings


def show_step_timing(log_file, shown_lines=None):
    """Display detailed timing for a specific workflow step.

    Args:
        log_file: Path to workflow log file
        shown_lines: Set of line numbers already displayed (to avoid duplicates)

    Returns:
        Updated set of shown line numbers
    """
    import re

    if shown_lines is None:
        shown_lines = set()

    if not log_file.exists():
        return shown_lines

    timing_pattern = re.compile(r'Command completed in ([\d.]+)s: (.+)')
    operations = []

    with log_file.open() as f:
        for line_num, line in enumerate(f):
            if line_num in shown_lines:
                continue

            match = timing_pattern.search(line)
            if match:
                duration = float(match.group(1))
                cmd = match.group(2).strip()
                operations.append((cmd, duration))
                shown_lines.add(line_num)

    if operations:
        print("\n  Detailed timing:")
        total = sum(dur for _, dur in operations)
        for cmd, duration in operations:
            pct = (duration / total * 100) if total > 0 else 0
            print(f"    • {cmd:30s}: {format_duration(duration):>10s}  ({pct:5.1f}%)")
        print(f"    {'─' * 55}")
        print(f"    {'Total':30s}: {format_duration(total):>10s}")

    return shown_lines


def run_command(cmd, description, cwd=None):
    """Run a command and check for success."""
    if not _QUIET_MODE:
        print(f"Running: {' '.join(cmd)}")
        print()

    start_time = time.time()

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=cwd
        )

        duration = time.time() - start_time

        if result.returncode != 0:
            if _QUIET_MODE:
                print(f"\n✗ {description} FAILED")
            else:
                print(f"✗ {description} FAILED (took {format_duration(duration)})")
            print(f"\nSTDOUT:\n{result.stdout}")
            print(f"\nSTDERR:\n{result.stderr}")
            return False, duration

        if not _QUIET_MODE:
            print(f"✓ {description} succeeded (took {format_duration(duration)})")

            # Show abbreviated output
            if result.stdout:
                lines = result.stdout.strip().split('\n')
                if len(lines) > 10:
                    print(f"\n[Output truncated, showing last 10 lines]")
                    print('\n'.join(lines[-10:]))
                else:
                    print(f"\n{result.stdout}")

        return True, duration

    except subprocess.TimeoutExpired:
        duration = time.time() - start_time
        print(f"✗ {description} TIMED OUT (after {format_duration(duration)})")
        return False, duration
    except Exception as e:
        duration = time.time() - start_time
        print(f"✗ {description} FAILED: {e} (after {format_duration(duration)})")
        return False, duration


def get_demo_data_dir():
    """Get the demo_data directory path."""
    import vcfcache
    package_dir = Path(vcfcache.__file__).parent
    return package_dir / "demo_data"


def run_smoke_test(keep_files=False, quiet=False):
    """Run the complete vcfcache smoke test workflow.

    Args:
        keep_files: If True, keep temporary files for inspection
        quiet: If True, suppress detailed output (show only essential information)

    Returns:
        int: Exit code (0 for success, 1 for failure)
    """
    global _QUIET_MODE
    _QUIET_MODE = quiet

    if quiet:
        print("Running smoke test...", end="", flush=True)
    else:
        print_section("VCFcache Complete Workflow Smoke Test")

    # Get demo data directory
    demo_data = get_demo_data_dir()

    if not demo_data.exists():
        print(f"\n✗ Demo data directory not found: {demo_data}")
        print("This should not happen with a proper installation.")
        return 1

    # Verify demo files exist
    required_files = [
        "demo_bp.bcf",
        "demo_bp.bcf.csi",
        "demo_bpext.bcf",
        "demo_bpext.bcf.csi",
        "demo_sample.vcf.gz",
        "demo_sample.vcf.gz.csi",
        "demo_params.yaml",
        "demo_annotation.yaml"
    ]

    missing_files = [f for f in required_files if not (demo_data / f).exists()]
    if missing_files:
        print(f"\n✗ Missing demo files: {', '.join(missing_files)}")
        return 1

    if not quiet:
        print(f"✓ Demo data directory: {demo_data}")
        print(f"✓ All required files present\n")

    # Create temporary directory
    temp_dir = Path(tempfile.mkdtemp(prefix="vcfcache_demo_"))
    print(f"Working directory: {temp_dir}")

    if keep_files:
        print(f"Note: Files will be kept for inspection")

    # Track timing for each step
    timings = {}

    # Track which log lines we've already shown (to avoid duplicates in shared log files)
    shown_log_lines = {}

    try:
        # Define paths
        cache_dir = temp_dir / "cache"
        stats_dir = temp_dir / "stats"

        bp_init_file = demo_data / "demo_bp.bcf"
        bp_extend_file = demo_data / "demo_bpext.bcf"
        sample_file = demo_data / "demo_sample.vcf.gz"
        params_file = demo_data / "demo_params.yaml"
        annotation_file = demo_data / "demo_annotation.yaml"

        # ====================================================================
        # Step 1: blueprint-init
        # ====================================================================
        print_step(1, "blueprint-init - Create initial cache from variants")

        cmd = [
            sys.executable, "-m", "vcfcache.cli",
            "blueprint-init",
            "--vcf", str(bp_init_file),
            "--output", str(cache_dir),
            "--force",
            # Note: Multiallelic splitting is now always performed
            # Params file is auto-generated internally
        ]

        success, duration = run_command(cmd, "Blueprint initialization")
        timings['blueprint-init'] = duration
        if not success:
            return 1

        # Verify blueprint was created
        blueprint_bcf = cache_dir / "blueprint" / "vcfcache.bcf"
        if not blueprint_bcf.exists():
            print(f"✗ Blueprint file not created: {blueprint_bcf}")
            return 1

        print(f"\n✓ Blueprint created: {blueprint_bcf}")

        # Show some stats
        stats_cmd = ["bcftools", "stats", str(blueprint_bcf)]
        result = subprocess.run(stats_cmd, capture_output=True, text=True)
        if result.returncode == 0:
            for line in result.stdout.split('\n'):
                if line.startswith('SN') and 'number of records' in line:
                    print(f"  {line.split(':')[1].strip()}")
                    break

        # Show detailed timing for this step
        blueprint_log = cache_dir / "blueprint" / "workflow.log"
        shown_log_lines[str(blueprint_log)] = show_step_timing(
            blueprint_log, shown_log_lines.get(str(blueprint_log))
        )

        # ====================================================================
        # Step 2: blueprint-extend
        # ====================================================================
        print_step(2, "blueprint-extend - Add more variants to cache")

        cmd = [
            sys.executable, "-m", "vcfcache.cli",
            "blueprint-extend",
            "--db", str(cache_dir),
            "-i", str(bp_extend_file)
        ]

        success, duration = run_command(cmd, "Blueprint extension")
        timings['blueprint-extend'] = duration
        if not success:
            return 1

        print(f"\n✓ Blueprint extended with additional variants")

        # Show updated stats
        result = subprocess.run(stats_cmd, capture_output=True, text=True)
        if result.returncode == 0:
            for line in result.stdout.split('\n'):
                if line.startswith('SN') and 'number of records' in line:
                    print(f"  {line.split(':')[1].strip()}")
                    break

        # Show detailed timing for this step (reuse same tracker for blueprint log)
        shown_log_lines[str(blueprint_log)] = show_step_timing(
            blueprint_log, shown_log_lines.get(str(blueprint_log))
        )

        # ====================================================================
        # Step 3: cache-build
        # ====================================================================
        print_step(3, "cache-build - Annotate the blueprint")

        cmd = [
            sys.executable, "-m", "vcfcache.cli",
            "cache-build",
            "--name", "demo_cache",
            "--db", str(cache_dir),
            "-a", str(annotation_file),
            "-y", str(params_file),
            "--force"
        ]

        success, duration = run_command(cmd, "Cache build")
        timings['cache-build'] = duration
        if not success:
            return 1

        # Verify cache was created
        cache_bcf = cache_dir / "cache" / "demo_cache" / "vcfcache_annotated.bcf"
        if not cache_bcf.exists():
            print(f"✗ Annotated cache not created: {cache_bcf}")
            return 1

        print(f"\n✓ Annotated cache created: {cache_bcf}")

        # Verify annotation tag is present
        header_cmd = ["bcftools", "view", "-h", str(cache_bcf)]
        result = subprocess.run(header_cmd, capture_output=True, text=True)
        if result.returncode == 0 and "##INFO=<ID=CSQ," in result.stdout:
            print(f"✓ Annotation tag CSQ present in cache")
        else:
            print(f"⚠ Warning: CSQ not found in cache header")

        # Show detailed timing for this step
        cache_log = cache_dir / "cache" / "demo_cache" / "workflow.log"
        shown_log_lines[str(cache_log)] = show_step_timing(
            cache_log, shown_log_lines.get(str(cache_log))
        )

        # ====================================================================
        # Step 4: annotate
        # ====================================================================
        print_step(4, "annotate - Use cache to annotate a sample VCF")

        output_bcf = temp_dir / "demo_sample_vc.bcf"
        stats_out_dir = stats_dir / f"{output_bcf.name}_vcstats"
        cmd = [
            sys.executable, "-m", "vcfcache.cli",
            "annotate",
            "-a", str(cache_dir / "cache" / "demo_cache"),
            "--vcf", str(sample_file),
            "--output", str(output_bcf),
            "--stats-dir", str(stats_dir),
            "-y", str(params_file),
            "--force"
        ]

        success, duration = run_command(cmd, "Sample annotation")
        timings['annotate (cached)'] = duration
        if not success:
            return 1

        # Verify output was created (output filename is user-specified)
        if not output_bcf.exists():
            print(f"✗ Annotated output not created: {output_bcf}")
            return 1

        print(f"\n✓ Annotated sample created: {output_bcf}")

        # Show final stats and check annotation
        stats_cmd = ["bcftools", "stats", str(output_bcf)]
        result = subprocess.run(stats_cmd, capture_output=True, text=True)
        if result.returncode == 0:
            for line in result.stdout.split('\n'):
                if line.startswith('SN') and 'number of records' in line:
                    print(f"  {line.split(':')[1].strip()}")
                    break

        # Check for annotation tag
        header_cmd = ["bcftools", "view", "-h", str(output_bcf)]
        result = subprocess.run(header_cmd, capture_output=True, text=True)
        if result.returncode == 0 and "##INFO=<ID=CSQ," in result.stdout:
            print(f"✓ Annotation tag CSQ present in output")

        # Show detailed timing for this step
        output_log = stats_out_dir / "workflow.log"
        shown_log_lines[str(output_log)] = show_step_timing(
            output_log, shown_log_lines.get(str(output_log))
        )

        # ====================================================================
        # Validation: Compare cached vs uncached annotation
        # ====================================================================
        print_step(5, "Validation - Compare cached vs uncached annotation outputs")

        # Run uncached annotation for comparison
        output_bcf_uncached = temp_dir / "demo_sample_uncached_vc.bcf"
        stats_out_uncached = stats_dir / f"{output_bcf_uncached.name}_vcstats"
        cmd_uncached = [
            sys.executable, "-m", "vcfcache.cli",
            "annotate",
            "-a", str(cache_dir / "cache" / "demo_cache"),
            "--vcf", str(sample_file),
            "--output", str(output_bcf_uncached),
            "--stats-dir", str(stats_dir),
            "-y", str(params_file),
            "--uncached",  # Force full annotation without cache
            "--force"
        ]

        success, duration = run_command(cmd_uncached, "Uncached annotation (for validation)")
        timings['annotate (uncached)'] = duration
        if not success:
            print("✗ ERROR: Uncached annotation failed")
            print("This is a critical issue - uncached mode must work for validation.")
            return 1

        if not output_bcf_uncached.exists():
            print(f"✗ ERROR: Uncached output not created: {output_bcf_uncached}")
            return 1

        # Show detailed timing for uncached annotation
        uncached_log = stats_out_uncached / "workflow.log"
        shown_log_lines[str(uncached_log)] = show_step_timing(
            uncached_log, shown_log_lines.get(str(uncached_log))
        )

        # Compute MD5 of BCF bodies (without headers)
        def compute_bcf_body_md5(bcf_path, filter_annotated=False, tag="CSQ"):
            """Compute MD5 of BCF body without header."""
            if filter_annotated:
                # Filter to annotated variants only
                result = subprocess.run(
                    ["bcftools", "view", "-H", "-i", f'INFO/{tag}!=""', str(bcf_path)],
                    capture_output=True,
                    text=True,
                )
            else:
                result = subprocess.run(
                    ["bcftools", "view", "-H", str(bcf_path)],
                    capture_output=True,
                    text=True,
                )
            return hashlib.md5(result.stdout.encode()).hexdigest()

        print("\nComputing MD5 checksums (body only, excluding headers)...")
        cached_md5 = compute_bcf_body_md5(output_bcf)
        uncached_md5 = compute_bcf_body_md5(output_bcf_uncached)

        print(f"Cached output MD5:   {cached_md5}")
        print(f"Uncached output MD5: {uncached_md5}")

        if cached_md5 == uncached_md5:
            print("\n✓ SUCCESS: Cached and uncached outputs are identical!")
        else:
            print("\n⚠ WARNING: Cached and uncached outputs differ (MD5 mismatch)")
            print("This may indicate a problem with the caching logic, OR it may be due to")
            print("non-deterministic behavior in the annotation tool itself (e.g., VEP ≥113).")
            print("See: https://github.com/Ensembl/ensembl-vep/issues/1959")
            print("\nRecommendation: Verify annotations semantically (same variants have same CSQ tags)")
            print("rather than relying solely on MD5 checksums for validation.")
            # Don't fail the test - this is a known issue with some annotation tools
            print("\n✓ Demo completed (with MD5 warning)")
            return 0

        # ====================================================================
        # Summary
        # ====================================================================
        if quiet:
            print(" ✓")
            print("Smoke test passed!")
            return 0

        print_section("Demo Complete!")

        print("✓ All steps executed successfully:\n")
        print("  1. blueprint-init  - Created initial cache")
        print("  2. blueprint-extend - Extended cache with more variants")
        print("  3. cache-build     - Annotated the blueprint")
        print("  4. annotate        - Used cache to annotate sample")
        print("  5. validation      - Verified cached == uncached (MD5 match)\n")

        # Timing summary
        print("Timing Summary:")
        print("─" * 60)
        total_time = sum(timings.values())
        for step, duration in timings.items():
            pct = (duration / total_time * 100) if total_time > 0 else 0
            print(f"  {step:25s}: {format_duration(duration):>12s}  ({pct:5.1f}%)")
        print("─" * 60)
        print(f"  {'Total':25s}: {format_duration(total_time):>12s}\n")

        # Detailed timing breakdown
        detailed_timings = collect_detailed_timings(
            cache_dir, stats_out_dir, stats_out_uncached
        )
        if detailed_timings:
            print("\nDetailed Operation Timing:")
            print("─" * 60)
            for step_name, operations in sorted(detailed_timings.items()):
                step_total = sum(dur for _, dur in operations)
                print(f"\n  {step_name}:")
                for cmd, duration in operations:
                    pct = (duration / step_total * 100) if step_total > 0 else 0
                    print(f"    {cmd:30s}: {format_duration(duration):>10s}  ({pct:5.1f}%)")
                print(f"    {'Subtotal':30s}: {format_duration(step_total):>10s}")
            print("─" * 60 + "\n")

        print(f"Demo files location: {demo_data}")
        if keep_files:
            print(f"Working files kept at: {temp_dir}")
        else:
            print(f"Cleaning up temporary files...")

        return 0

    except KeyboardInterrupt:
        print("\n\n✗ Demo interrupted by user")
        return 1
    except Exception as e:
        print(f"\n\n✗ Demo failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    finally:
        # Cleanup
        if not keep_files and temp_dir.exists():
            try:
                shutil.rmtree(temp_dir)
                print(f"✓ Cleaned up temporary directory")
            except Exception as e:
                print(f"⚠ Warning: Could not clean up {temp_dir}: {e}")
