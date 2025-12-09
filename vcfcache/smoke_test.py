"""Smoke test for vcfcache package.

This module provides a minimal smoke test that can be run after installation
to verify that vcfcache is working correctly.

Usage:
    python -m vcfcache.smoke_test

Or from Python:
    from vcfcache.smoke_test import run_smoke_test
    run_smoke_test()
"""

import sys
import subprocess
from pathlib import Path


def test_imports():
    """Test that all core modules can be imported."""
    print("Testing imports...")
    try:
        import vcfcache
        from vcfcache import cli
        from vcfcache.database import base, annotator, initializer, updater
        from vcfcache.utils import validation, paths, logging as vcf_logging
        print("  ✓ All core modules imported successfully")
        return True
    except ImportError as e:
        print(f"  ✗ Import failed: {e}")
        return False


def test_version():
    """Test that version is accessible."""
    print("Testing version...")
    try:
        import vcfcache
        version = vcfcache.__version__
        print(f"  ✓ vcfcache version: {version}")
        return True
    except Exception as e:
        print(f"  ✗ Version check failed: {e}")
        return False


def test_bcftools():
    """Test that bcftools can be detected."""
    print("Testing bcftools detection...")
    try:
        from vcfcache.utils.validation import check_bcftools_installed, MIN_BCFTOOLS_VERSION
        bcftools_path = check_bcftools_installed()
        print(f"  ✓ Found bcftools at: {bcftools_path}")
        print(f"    Minimum required version: {MIN_BCFTOOLS_VERSION}")

        # Get actual version
        result = subprocess.run(
            [bcftools_path, "--version-only"],
            capture_output=True,
            text=True,
            check=True,
            timeout=5
        )
        actual_version = result.stdout.strip()
        print(f"    Installed version: {actual_version}")
        return True
    except FileNotFoundError as e:
        print(f"  ✗ bcftools not found: {e}")
        print("    Install bcftools >= 1.20:")
        print("      Ubuntu/Debian: sudo apt-get install bcftools")
        print("      macOS: brew install bcftools")
        print("      Conda: conda install -c bioconda bcftools")
        return False
    except Exception as e:
        print(f"  ✗ bcftools check failed: {e}")
        return False


def test_cli_help():
    """Test that CLI help command works."""
    print("Testing CLI help command...")

    # Try direct command first
    try:
        result = subprocess.run(
            ["vcfcache", "--help"],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0 and ("vcfcache" in result.stdout.lower() or "annotation" in result.stdout.lower()):
            print("  ✓ CLI help command works")
            return True
        else:
            print(f"  ✗ CLI help failed with return code {result.returncode}")
            return False
    except FileNotFoundError:
        # Fall back to module invocation
        try:
            result = subprocess.run(
                [sys.executable, "-m", "vcfcache.cli", "--help"],
                capture_output=True,
                text=True,
                timeout=10
            )
            # When run as module, it says "cli.py" not "vcfcache"
            if result.returncode == 0 and "annotation" in result.stdout.lower():
                print("  ✓ CLI help command works (via module)")
                return True
            else:
                print(f"  ✗ CLI help failed with return code {result.returncode}")
                return False
        except Exception as e:
            print(f"  ✗ CLI help test failed: {e}")
            return False
    except Exception as e:
        print(f"  ✗ CLI help test failed: {e}")
        return False


def test_recipes():
    """Test that recipe files are accessible."""
    print("Testing recipe files...")
    try:
        # Recipes are inside the vcfcache package
        import vcfcache
        package_dir = Path(vcfcache.__file__).parent
        recipes_dir = package_dir / "recipes"

        if recipes_dir.exists() and recipes_dir.is_dir():
            recipe_count = len(list(recipes_dir.rglob("*.yaml")))
            print(f"  ✓ Found recipes directory with {recipe_count} YAML files")
            print(f"    Location: {recipes_dir}")
            return True
        else:
            print(f"  ✗ Recipes directory not found at: {recipes_dir}")
            return False
    except Exception as e:
        print(f"  ✗ Recipe check failed: {e}")
        return False


def run_smoke_test():
    """Run all smoke tests and return success status."""
    print("\n" + "="*70)
    print("VCFcache Smoke Test")
    print("="*70 + "\n")

    tests = [
        ("Imports", test_imports),
        ("Version", test_version),
        ("bcftools", test_bcftools),
        ("CLI", test_cli_help),
        ("Recipes", test_recipes),
    ]

    results = {}
    for test_name, test_func in tests:
        results[test_name] = test_func()
        print()

    # Summary
    print("="*70)
    print("Summary:")
    print("-"*70)

    passed = sum(results.values())
    total = len(results)

    for test_name, result in results.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"  {status}: {test_name}")

    print("-"*70)
    print(f"Results: {passed}/{total} tests passed")
    print("="*70 + "\n")

    if passed == total:
        print("🎉 All smoke tests passed! vcfcache is ready to use.")
        return 0
    else:
        print("⚠️  Some smoke tests failed. Please check the errors above.")
        return 1


def main():
    """Entry point for command-line execution."""
    exit_code = run_smoke_test()
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
