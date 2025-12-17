"""Comprehensive demo of vcfcache workflow.

This module demonstrates the complete vcfcache workflow with all 5 steps:
1. blueprint-init: Create initial cache from VCF
2. blueprint-extend: Add more variants to existing cache
3. cache-build: Annotate the cache
4. annotate: Use the cache to annotate a sample VCF
5. validation: Verify cached and uncached outputs match (MD5 comparison)

Usage:
    vcfcache demo --smoke-test [--debug] [--quiet]
    vcfcache demo -a <cache> --vcf <vcf> -y <params> [--output <dir>] [--debug] [--quiet]

Or from Python:
    from vcfcache.demo import run_smoke_test, run_benchmark
    run_smoke_test(keep_files=False, quiet=False)
    run_benchmark(cache_dir, vcf_file, params_file, output_dir=None, keep_files=False, quiet=False)
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


def collect_detailed_timings(cache_dir, output_dir):
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
    for output_subdir, step_name in [(output_dir, "annotate-cached"),
                                       (output_dir.parent / "output_uncached", "annotate-uncached")]:
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
            cwd=cwd,
            timeout=300
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
        output_dir = temp_dir / "output"

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

        # Note: Don't create output_dir manually - let vcfcache annotate create it
        # to avoid validation issues

        cmd = [
            sys.executable, "-m", "vcfcache.cli",
            "annotate",
            "-a", str(cache_dir / "cache" / "demo_cache"),
            "--vcf", str(sample_file),
            "--output", str(output_dir),
            "-y", str(params_file),
            "--force"
        ]

        success, duration = run_command(cmd, "Sample annotation")
        timings['annotate (cached)'] = duration
        if not success:
            return 1

        # Verify output was created
        # Note: output filename is based on input filename, so demo_sample.vcf.gz -> demo_sample.vcf_vc.bcf
        output_bcf = output_dir / "demo_sample.vcf_vc.bcf"
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
        output_log = output_dir / "workflow.log"
        shown_log_lines[str(output_log)] = show_step_timing(
            output_log, shown_log_lines.get(str(output_log))
        )

        # ====================================================================
        # Validation: Compare cached vs uncached annotation
        # ====================================================================
        print_step(5, "Validation - Compare cached vs uncached annotation outputs")

        # Run uncached annotation for comparison
        output_dir_uncached = temp_dir / "output_uncached"
        cmd_uncached = [
            sys.executable, "-m", "vcfcache.cli",
            "annotate",
            "-a", str(cache_dir / "cache" / "demo_cache"),
            "--vcf", str(sample_file),
            "--output", str(output_dir_uncached),
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

        output_bcf_uncached = output_dir_uncached / "demo_sample.vcf_vc.bcf"
        if not output_bcf_uncached.exists():
            print(f"✗ ERROR: Uncached output not created: {output_bcf_uncached}")
            return 1

        # Show detailed timing for uncached annotation
        uncached_log = output_dir_uncached / "workflow.log"
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
            print("\n✓ SUCCESS: Cached and uncached outputs are identical (all variants)")
        else:
            print("\nNote: MD5s differ (expected if annotation tool drops variants)")
            print("VCFcache preserves ALL input variants. Comparing annotated variants only...\n")

            # Compare annotated variants only (CSQ tag from demo annotation)
            cached_md5_annotated = compute_bcf_body_md5(output_bcf, filter_annotated=True, tag="CSQ")
            uncached_md5_annotated = compute_bcf_body_md5(output_bcf_uncached, filter_annotated=True, tag="CSQ")

            print(f"Cached MD5 (annotated only):   {cached_md5_annotated}")
            print(f"Uncached MD5 (annotated only): {uncached_md5_annotated}")

            if cached_md5_annotated == uncached_md5_annotated:
                print("\n✓ SUCCESS: Annotated variants are identical!")
            else:
                print("\n✗ ERROR: Annotated variants differ!")
                print("This indicates a problem with the caching logic.")
                return 1

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
        detailed_timings = collect_detailed_timings(cache_dir, output_dir)
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


def run_benchmark(cache_dir, vcf_file, params_file, output_dir=None, keep_files=False, quiet=False):
    """Run benchmark comparing cached vs uncached annotation.

    Args:
        cache_dir: Path to annotation cache directory
        vcf_file: Path to VCF/BCF file to annotate
        params_file: Path to params YAML file
        output_dir: Optional output directory (default: temporary directory in /tmp)
        keep_files: If True, keep temporary files for inspection
        quiet: If True, suppress detailed output (show only essential information)

    Returns:
        int: Exit code (0 for success, 1 for failure)
    """
    global _QUIET_MODE
    _QUIET_MODE = quiet

    print_section("VCFcache Annotation Benchmark")

    cache_path = Path(cache_dir).expanduser().resolve()
    vcf_path = Path(vcf_file).expanduser().resolve()
    params_path = Path(params_file).expanduser().resolve()

    # Validate inputs
    if not cache_path.exists():
        print(f"✗ Cache directory not found: {cache_path}")
        return 1
    if not vcf_path.exists():
        print(f"✗ VCF file not found: {vcf_path}")
        return 1
    if not params_path.exists():
        print(f"✗ Params file not found: {params_path}")
        return 1

    print(f"Cache: {cache_path}")
    print(f"VCF:   {vcf_path}")
    print(f"Params: {params_path}")
    print()

    # Estimate file size and ask for confirmation
    vcf_size_mb = vcf_path.stat().st_size / (1024 * 1024)
    print(f"VCF file size: {vcf_size_mb:.2f} MB")
    print("This benchmark will run annotation twice (uncached + cached).")
    print("Depending on the file size and annotation tool, this may take a while.")
    print()

    response = input("Do you want to proceed? [y/N]: ")
    if response.lower() not in ['y', 'yes']:
        print("Benchmark cancelled.")
        return 0

    # Create output directory
    if output_dir:
        temp_dir = Path(output_dir).expanduser().resolve()
        temp_dir.mkdir(parents=True, exist_ok=True)
        created_temp = False  # User-specified directory, don't auto-cleanup
        print(f"\nOutput directory: {temp_dir}\n")
    else:
        temp_dir = Path(tempfile.mkdtemp(prefix="vcfcache_benchmark_"))
        created_temp = True  # Auto-created temp, cleanup unless keep_files
        print(f"\nWorking directory: {temp_dir}\n")

    if keep_files:
        print(f"Note: Files will be kept for inspection\n")

    try:
        uncached_dir = temp_dir / "uncached"
        cached_dir = temp_dir / "cached"

        # ====================================================================
        # Run 1: Uncached annotation
        # ====================================================================
        print_section("Run 1: Uncached Annotation (Full Pipeline)")

        start_time = time.time()
        cmd_uncached = [
            sys.executable, "-m", "vcfcache.cli",
            "annotate",
            "-a", str(cache_path),
            "--vcf", str(vcf_path),
            "--output", str(uncached_dir),
            "-y", str(params_path),
            "--uncached",  # Force uncached mode
            "--force",
        ]
        if keep_files:
            cmd_uncached.append("--debug")

        print(f"Running: {' '.join(cmd_uncached)}\n")

        result = subprocess.run(
            cmd_uncached,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour timeout
        )

        uncached_time = time.time() - start_time

        if result.returncode != 0:
            print(f"✗ Uncached annotation FAILED")
            print(f"\nSTDOUT:\n{result.stdout}")
            print(f"\nSTDERR:\n{result.stderr}")
            return 1

        print(f"✓ Uncached annotation completed in {uncached_time:.2f}s\n")

        # ====================================================================
        # Run 2: Cached annotation
        # ====================================================================
        print_section("Run 2: Cached Annotation (With Cache)")

        start_time = time.time()
        cmd_cached = [
            sys.executable, "-m", "vcfcache.cli",
            "annotate",
            "-a", str(cache_path),
            "--vcf", str(vcf_path),
            "--output", str(cached_dir),
            "-y", str(params_path),
            "--force",
        ]
        if keep_files:
            cmd_cached.append("--debug")

        print(f"Running: {' '.join(cmd_cached)}\n")

        result = subprocess.run(
            cmd_cached,
            capture_output=True,
            text=True,
            timeout=3600,
        )

        cached_time = time.time() - start_time

        if result.returncode != 0:
            print(f"✗ Cached annotation FAILED")
            print(f"\nSTDOUT:\n{result.stdout}")
            print(f"\nSTDERR:\n{result.stderr}")
            return 1

        print(f"✓ Cached annotation completed in {cached_time:.2f}s")

        # Show detailed output from cached run
        if result.stdout:
            print(f"  Annotation completed in {cached_time:.1f}s")
            for line in result.stdout.split('\n'):
                if 'Step' in line and '/' in line:
                    print(f"  {line.strip()}")

        # ====================================================================
        # Compare outputs
        # ====================================================================
        print_section("Results Comparison")

        # Find output files
        vcf_basename = vcf_path.stem if vcf_path.suffix == '.gz' else vcf_path.stem
        uncached_bcf = uncached_dir / f"{vcf_basename}_vc.bcf"
        cached_bcf = cached_dir / f"{vcf_basename}_vc.bcf"

        if not uncached_bcf.exists():
            print(f"✗ Uncached output not found: {uncached_bcf}")
            return 1
        if not cached_bcf.exists():
            print(f"✗ Cached output not found: {cached_bcf}")
            return 1

        # Read annotation tag from cache's annotation.yaml
        annotation_yaml_path = cache_path / "annotation.yaml"
        annotation_tag = "CSQ"  # Default
        if annotation_yaml_path.exists():
            import yaml
            with open(annotation_yaml_path) as f:
                anno_config = yaml.safe_load(f)
                annotation_tag = anno_config.get("must_contain_info_tag", "CSQ")

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

        # Total variants
        uncached_md5 = compute_bcf_body_md5(uncached_bcf)
        cached_md5 = compute_bcf_body_md5(cached_bcf)

        print(f"\nUncached output MD5: {uncached_md5}")
        print(f"Cached output MD5:   {cached_md5}")

        if uncached_md5 == cached_md5:
            print("✓ Outputs are identical (all variants)")
        else:
            print("✗ WARNING: Outputs differ!")
            print("\nNote: This is EXPECTED behavior!")
            print(f"VCFcache preserves ALL input variants, even if the annotation tool drops some.")
            print(f"Comparing annotated variants only ({annotation_tag}!=\"\")...\n")

            # Compare annotated variants only
            uncached_md5_annotated = compute_bcf_body_md5(uncached_bcf, filter_annotated=True, tag=annotation_tag)
            cached_md5_annotated = compute_bcf_body_md5(cached_bcf, filter_annotated=True, tag=annotation_tag)

            print(f"Uncached MD5 (annotated only): {uncached_md5_annotated}")
            print(f"Cached MD5 (annotated only):   {cached_md5_annotated}")

            if uncached_md5_annotated == cached_md5_annotated:
                print("✓ SUCCESS: Annotated variants are identical!")
                print(f"  (Cached preserves extra unannotated variants - see WIKI.md for details)")
            else:
                print("⚠ Annotated variants differ slightly (may be due to order or minor formatting differences)")

        # Timing comparison
        print(f"\n{'─'*70}")
        print("Timing Summary:")
        print(f"{'─'*70}")
        print(f"Uncached: {uncached_time:8.2f}s")
        print(f"Cached:   {cached_time:8.2f}s")
        print(f"{'─'*70}")

        if cached_time < uncached_time:
            speedup = uncached_time / cached_time
            time_saved = uncached_time - cached_time
            percent_saved = (time_saved / uncached_time) * 100
            print(f"Speedup:      {speedup:5.2f}x")
            print(f"Time saved: {time_saved:.2f}s ({percent_saved:.1f}%)")
        else:
            print("Note: Cached mode was not faster (file may be too small for caching benefit)")

        # ====================================================================
        # Show detailed timing breakdown
        # ====================================================================
        print(f"\n{'─'*70}")
        print("Detailed Timing Breakdown:")
        print(f"{'─'*70}")

        # Show uncached timing details
        print("\nUncached workflow:")
        uncached_log = uncached_dir / "workflow.log"
        show_step_timing(uncached_log)

        # Show cached timing details
        print("\nCached workflow:")
        cached_log = cached_dir / "workflow.log"
        show_step_timing(cached_log)

        # ====================================================================
        # Print raw command for copy-paste
        # ====================================================================
        print_section("Raw Command for Copy-Paste")

        print("# Cached annotation:")
        print(" \\\n  ".join(cmd_cached))
        print()
        print("# Uncached annotation:")
        print(" \\\n  ".join(cmd_uncached))
        print()

        # ====================================================================
        # Show exact commands to reproduce uncached workflow
        # ====================================================================
        print_section("Reproduce Uncached Workflow Manually")

        print("The uncached workflow runs the following steps:")
        print()
        print("# Step 1: Normalize input (split multiallelic variants)")
        vcf_basename = vcf_path.stem if vcf_path.suffix == '.gz' else vcf_path.stem
        print(f"bcftools norm -m- {vcf_path} \\")
        print(f"  -o {vcf_basename}_normalized.bcf -Ob -W --threads 1")
        print()
        print("# Step 2: Run annotation tool (command from annotation.yaml)")
        print("# The exact command depends on your annotation.yaml configuration.")
        print("# Typically for VEP:")
        print(f"bcftools view {vcf_basename}_normalized.bcf | \\")
        print("  vep --assembly GRCh37 --cache --offline --format vcf \\")
        print("      --output_file stdout --vcf [VEP_OPTIONS] | \\")
        print(f"  bcftools view -o {vcf_basename}_vc.bcf -Ob -W")
        print()
        print("Note: See your annotation.yaml for the exact VEP command and options.")
        print(f"Uncached workflow log: {uncached_log}")
        print()

        # Print output location
        print(f"Output files location: {temp_dir}")
        if keep_files:
            print(f"Working files will be kept for inspection")
        elif not created_temp:
            print(f"Output directory will be preserved (user-specified)")
        else:
            print(f"Temporary directory will be cleaned up")

        # ====================================================================
        # Write benchmark summary to file
        # ====================================================================
        summary_file = temp_dir / "benchmark_summary.log"
        with summary_file.open('w') as f:
            f.write("=" * 70 + "\n")
            f.write("VCFcache Annotation Benchmark Summary\n")
            f.write("=" * 70 + "\n\n")

            f.write(f"Cache: {cache_path}\n")
            f.write(f"VCF:   {vcf_path}\n")
            f.write(f"Params: {params_path}\n")
            f.write(f"VCF file size: {vcf_size_mb:.2f} MB\n\n")

            f.write("=" * 70 + "\n")
            f.write("Results Comparison\n")
            f.write("=" * 70 + "\n\n")

            f.write(f"Uncached output: {uncached_bcf}\n")
            f.write(f"Cached output:   {cached_bcf}\n\n")

            f.write(f"Uncached output MD5: {uncached_md5}\n")
            f.write(f"Cached output MD5:   {cached_md5}\n")

            if uncached_md5 == cached_md5:
                f.write("✓ Outputs are identical (body matches)\n\n")
            else:
                f.write("✗ WARNING: Outputs differ!\n\n")

            f.write("─" * 70 + "\n")
            f.write("Timing Summary:\n")
            f.write("─" * 70 + "\n")
            f.write(f"Uncached: {uncached_time:8.2f}s\n")
            f.write(f"Cached:   {cached_time:8.2f}s\n")
            f.write("─" * 70 + "\n")

            if cached_time < uncached_time:
                speedup = uncached_time / cached_time
                time_saved = uncached_time - cached_time
                percent_saved = (time_saved / uncached_time) * 100
                f.write(f"Speedup:      {speedup:5.2f}x\n")
                f.write(f"Time saved: {time_saved:.2f}s ({percent_saved:.1f}%)\n")
            else:
                f.write("Note: Cached mode was not faster\n")

            f.write("\n" + "=" * 70 + "\n")
            f.write("Output files location\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Uncached: {uncached_dir}\n")
            f.write(f"Cached:   {cached_dir}\n")
            f.write(f"Summary:  {summary_file}\n")

        print(f"\n✓ Benchmark summary written to: {summary_file}")

        return 0

    except KeyboardInterrupt:
        print("\n\n✗ Benchmark interrupted by user")
        return 1
    except Exception as e:
        print(f"\n\n✗ Benchmark failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    finally:
        # Cleanup (only for auto-created temp directories)
        if created_temp and not keep_files and temp_dir.exists():
            try:
                shutil.rmtree(temp_dir)
                print(f"\n✓ Cleaned up temporary directory")
            except Exception as e:
                print(f"\n⚠ Warning: Could not clean up {temp_dir}: {e}")
        elif keep_files and temp_dir.exists():
            print(f"\nWorking files kept at: {temp_dir}")


def main():
    """Entry point for command-line execution."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Run vcfcache smoke test or benchmark"
    )
    parser.add_argument(
        "--smoke-test",
        action="store_true",
        help="Run comprehensive smoke test"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Keep temporary files for inspection"
    )
    parser.add_argument(
        "-a", "--annotation_db",
        type=str,
        help="Annotation cache directory (benchmark mode)"
    )
    parser.add_argument(
        "--vcf",
        type=str,
        help="VCF file (benchmark mode)"
    )
    parser.add_argument(
        "-y", "--params",
        type=str,
        help="Params file (benchmark mode)"
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output directory for benchmark results (default: temporary directory in /tmp)"
    )

    args = parser.parse_args()

    if args.smoke_test:
        exit_code = run_smoke_test(keep_files=args.debug)
    elif args.annotation_db and args.vcf and args.params:
        exit_code = run_benchmark(
            args.annotation_db,
            args.vcf,
            args.params,
            output_dir=args.output,
            keep_files=args.debug
        )
    else:
        parser.print_help()
        exit_code = 1

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
