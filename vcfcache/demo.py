"""Comprehensive demo of vcfcache workflow.

This module demonstrates the complete vcfcache workflow with all 4 commands:
1. blueprint-init: Create initial cache from VCF
2. blueprint-extend: Add more variants to existing cache
3. cache-build: Annotate the cache
4. annotate: Use the cache to annotate a sample VCF

Usage:
    vcfcache demo [--keep-files]

Or from Python:
    from vcfcache.demo import run_demo
    run_demo()
"""

import sys
import subprocess
import tempfile
import shutil
from pathlib import Path


def print_section(title):
    """Print a section header."""
    print("\n" + "="*70)
    print(f"  {title}")
    print("="*70 + "\n")


def print_step(step_num, description):
    """Print a step header."""
    print(f"\n{'─'*70}")
    print(f"Step {step_num}: {description}")
    print(f"{'─'*70}\n")


def run_command(cmd, description, cwd=None):
    """Run a command and check for success."""
    print(f"Running: {' '.join(cmd)}")
    print()

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=cwd,
            timeout=300
        )

        if result.returncode != 0:
            print(f"✗ {description} FAILED")
            print(f"\nSTDOUT:\n{result.stdout}")
            print(f"\nSTDERR:\n{result.stderr}")
            return False

        print(f"✓ {description} succeeded")

        # Show abbreviated output
        if result.stdout:
            lines = result.stdout.strip().split('\n')
            if len(lines) > 10:
                print(f"\n[Output truncated, showing last 10 lines]")
                print('\n'.join(lines[-10:]))
            else:
                print(f"\n{result.stdout}")

        return True

    except subprocess.TimeoutExpired:
        print(f"✗ {description} TIMED OUT")
        return False
    except Exception as e:
        print(f"✗ {description} FAILED: {e}")
        return False


def get_demo_data_dir():
    """Get the demo_data directory path."""
    import vcfcache
    package_dir = Path(vcfcache.__file__).parent
    return package_dir / "demo_data"


def run_demo(keep_files=False):
    """Run the complete vcfcache demo workflow.

    Args:
        keep_files: If True, keep temporary files for inspection

    Returns:
        int: Exit code (0 for success, 1 for failure)
    """
    print_section("VCFcache Complete Workflow Demo")

    # Get demo data directory
    demo_data = get_demo_data_dir()

    if not demo_data.exists():
        print(f"✗ Demo data directory not found: {demo_data}")
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
        print(f"✗ Missing demo files: {', '.join(missing_files)}")
        return 1

    print(f"✓ Demo data directory: {demo_data}")
    print(f"✓ All required files present\n")

    # Create temporary directory
    temp_dir = Path(tempfile.mkdtemp(prefix="vcfcache_demo_"))
    print(f"Working directory: {temp_dir}")

    if keep_files:
        print(f"Note: Files will be kept for inspection")

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
            "-y", str(params_file),
            "--force"
        ]

        if not run_command(cmd, "Blueprint initialization"):
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

        if not run_command(cmd, "Blueprint extension"):
            return 1

        print(f"\n✓ Blueprint extended with additional variants")

        # Show updated stats
        result = subprocess.run(stats_cmd, capture_output=True, text=True)
        if result.returncode == 0:
            for line in result.stdout.split('\n'):
                if line.startswith('SN') and 'number of records' in line:
                    print(f"  {line.split(':')[1].strip()}")
                    break

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

        if not run_command(cmd, "Cache build"):
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
        if result.returncode == 0 and "DEMO_TAG" in result.stdout:
            print(f"✓ Annotation tag DEMO_TAG present in cache")
        else:
            print(f"⚠ Warning: DEMO_TAG not found in cache header")

        # ====================================================================
        # Step 4: annotate
        # ====================================================================
        print_step(4, "annotate - Use cache to annotate a sample VCF")

        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            sys.executable, "-m", "vcfcache.cli",
            "annotate",
            "-a", str(cache_dir / "cache" / "demo_cache"),
            "--vcf", str(sample_file),
            "--output", str(output_dir),
            "-y", str(params_file),
            "--force"
        ]

        if not run_command(cmd, "Sample annotation"):
            return 1

        # Verify output was created
        output_bcf = output_dir / "demo_sample_vst.bcf"
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
        if result.returncode == 0 and "DEMO_TAG" in result.stdout:
            print(f"✓ Annotation tag DEMO_TAG present in output")

        # ====================================================================
        # Summary
        # ====================================================================
        print_section("Demo Complete!")

        print("✓ All 4 commands executed successfully:\n")
        print("  1. blueprint-init  - Created initial cache")
        print("  2. blueprint-extend - Extended cache with more variants")
        print("  3. cache-build     - Annotated the blueprint")
        print("  4. annotate        - Used cache to annotate sample\n")

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


def main():
    """Entry point for command-line execution."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Run a comprehensive demo of the vcfcache workflow"
    )
    parser.add_argument(
        "--keep-files",
        action="store_true",
        help="Keep temporary files for inspection"
    )

    args = parser.parse_args()

    exit_code = run_demo(keep_files=args.keep_files)
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
