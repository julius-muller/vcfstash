# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""VCF Annotation Cache.

This script manages a database of genetic variants in BCF/VCF format,
providing functionality to initialize, add to, and annotate variant databases.

Key features:
- Supports BCF/VCF format (uncompressed/compressed)
- Requires pre-indexed input files (CSI/TBI index)
  -> This ensures input is properly sorted and valid
- Maintains database integrity through MD5 checksums
- Provides versioned annotation workflow support
- Includes detailed logging of all operations

Author: Julius Müller, PhD
Organization: GHGA - German Human Genome-Phenome Archive
Date: 16-03-2025
"""

import argparse
import os
import sys
import subprocess
import re
from importlib.metadata import version as pkg_version
from pathlib import Path

import requests
import yaml

from vcfcache.integrations.zenodo import (
    download_doi,
    search_zenodo_records,
    ZenodoError,
)
from vcfcache.database.annotator import DatabaseAnnotator, VCFAnnotator
from vcfcache.database.initializer import DatabaseInitializer
from vcfcache.database.updater import DatabaseUpdater
from vcfcache.utils.logging import log_command, setup_logging
from vcfcache.utils.archive import extract_cache, tar_cache_subset
from vcfcache.utils.paths import get_project_root
from vcfcache.utils.validation import check_bcftools_installed

# Ensure VCFCACHE_ROOT is set (used by packaged resources/recipes)
os.environ.setdefault("VCFCACHE_ROOT", str(get_project_root()))


def _load_dotenv() -> None:
    """Load environment variables from .env file if present.

    Checks for .env in:
    1. User's home directory (~/.env)
    2. Current working directory (./.env) - takes precedence

    Only sets variables that aren't already in os.environ.
    """
    env_files = [
        Path.home() / ".env",
        Path.cwd() / ".env",
    ]

    for env_path in env_files:
        if not env_path.exists():
            continue

        for line in env_path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, val = line.split("=", 1)
            key = key.strip()
            val = val.strip().strip("'").strip('"')
            if key and key not in os.environ:
                os.environ[key] = val


def _show_detailed_timings(workflow_log: Path) -> None:
    """Display detailed timing information from workflow log.

    Args:
        workflow_log: Path to workflow.log file
    """
    import re

    if not workflow_log.exists():
        return

    timing_pattern = re.compile(r'Command completed in ([\d.]+)s: (.+)')
    operations = []

    with workflow_log.open() as f:
        for line in f:
            match = timing_pattern.search(line)
            if match:
                duration = float(match.group(1))
                cmd = match.group(2).strip()
                operations.append((cmd, duration))

    if operations:
        print("\n  Detailed timing:")
        total = sum(dur for _, dur in operations)
        for cmd, duration in operations:
            pct = (duration / total * 100) if total > 0 else 0
            minutes = int(duration // 60)
            secs = duration % 60
            if minutes > 0:
                time_str = f"{minutes}m {secs:.3f}s"
            else:
                time_str = f"{secs:.3f}s"
            print(f"    • {cmd:30s}: {time_str:>10s}  ({pct:5.1f}%)")
        print(f"    {'─' * 55}")
        # Format total time
        total_minutes = int(total // 60)
        total_secs = total % 60
        if total_minutes > 0:
            total_str = f"{total_minutes}m {total_secs:.3f}s"
        else:
            total_str = f"{total_secs:.3f}s"
        print(f"    {'Total':30s}: {total_str:>10s}")


def _print_annotation_command(path_hint: Path, params_override: Path | None = None) -> None:
    """Print cache requirements and the stored annotation command.

    Args:
        path_hint: Path to cache root, cache directory, or specific annotation directory.
    """
    # Try to find the annotation.yaml file
    # First check if path_hint itself has annotation.yaml (specific cache directory)
    params_file = path_hint / "annotation.yaml"
    cache_dir = path_hint

    if not params_file.exists():
        # Try to find cache directory and list available caches
        try:
            cache_dir = _find_cache_dir(path_hint)
            caches = [c for c in cache_dir.iterdir() if c.is_dir()]

            if not caches:
                raise FileNotFoundError(f"No annotation caches found under {cache_dir}")

            if len(caches) == 1:
                # Only one cache, use it
                params_file = caches[0] / "annotation.yaml"
                cache_dir = caches[0]
                if not params_file.exists():
                    raise FileNotFoundError(
                        f"Annotation config not found: {params_file} (cache may be incomplete)"
                    )
            else:
                # Multiple caches, ask user to specify
                print(f"Multiple caches found. Please specify which one:")
                for cache in sorted(caches):
                    status = "" if (cache / "vcfcache_annotated.bcf").exists() else " (incomplete)"
                    print(f"  vcfcache annotate --requirements -a {cache}{status}")
                return
        except Exception as e:
            raise FileNotFoundError(
                f"Could not find annotation.yaml in {path_hint}. "
                f"Please provide path to a specific cache directory. Error: {e}"
            )

    anno_text = params_file.read_text()
    anno = yaml.safe_load(anno_text) or {}
    params = {}
    params_snapshot = cache_dir / "params.snapshot.yaml"
    params_path = params_override if params_override else params_snapshot
    if params_path and params_path.exists():
        params = yaml.safe_load(params_path.read_text()) or {}

    # Try new format (annotation_cmd) first, then fall back to old format (annotation_tool_cmd)
    command = anno.get("annotation_cmd") or anno.get("annotation_tool_cmd")
    if not command:
        raise ValueError(
            "annotation_cmd or annotation_tool_cmd not found in annotation.yaml; cache may be incomplete"
        )

    use_color = sys.stdout.isatty() and os.environ.get("NO_COLOR") is None
    C_RESET = "\033[0m" if use_color else ""
    C_BOLD = "\033[1m" if use_color else ""
    C_CYAN = "\033[36m" if use_color else ""
    C_GREEN = "\033[32m" if use_color else ""
    C_RED = "\033[31m" if use_color else ""

    def _hdr(text: str) -> str:
        return f"{C_BOLD}{C_CYAN}{text}{C_RESET}"

    def _ok(text: str) -> str:
        return f"{C_GREEN}{text}{C_RESET}"

    def _bad(text: str) -> str:
        return f"{C_RED}{text}{C_RESET}"

    # Extract required ${params.*} keys from annotation_cmd
    required_keys = sorted(set(re.findall(r"\$\{params\.([A-Za-z0-9_\\.]+)\}", command)))

    print(_hdr("Cache annotation recipe"))
    print(f"  annotation.yaml: {params_file}")
    print("\n" + anno_text.strip() + "\n")

    print(_hdr("Cache requirements (from annotation.yaml)"))
    if anno.get("genome_build"):
        print(f"  genome_build: {anno.get('genome_build')}")
    if anno.get("must_contain_info_tag"):
        print(f"  must_contain_info_tag: {anno.get('must_contain_info_tag')}")
    if anno.get("required_tool_version"):
        print(f"  required_tool_version: {anno.get('required_tool_version')}")

    print("\n" + _hdr("Required params (extracted from annotation_cmd)"))
    if required_keys:
        for key in required_keys:
            print(f"  {key}")
    else:
        print("  (none)")

    print("\n" + _hdr("Params evaluation"))
    if params_path and params_path.exists():
        if params_override is None:
            print(f"  params.yaml: {params_path} (from cache snapshot; no -y provided)")
        else:
            print(f"  params.yaml: {params_path}")
    else:
        print("  params.yaml: (missing)")

    if required_keys:
        for key in required_keys:
            parts = key.split(".")
            cur = params
            for part in parts:
                if not isinstance(cur, dict) or part not in cur:
                    cur = None
                    break
                cur = cur[part]
            if cur is None:
                print(f"  {key}: {_bad('<missing>')}")
            else:
                print(f"  {key}: {cur}")
    else:
        print("  (no params required)")

    print("\n" + _hdr("Tool checks"))
    bcftools_cmd = params.get("bcftools_cmd") or "bcftools"
    try:
        res = subprocess.run([str(bcftools_cmd), "--version-only"], capture_output=True, text=True)
        if res.returncode == 0:
            print(f"  bcftools: {res.stdout.strip()}  {_ok('✓')}")
        else:
            print(f"  bcftools: ERROR ({res.stderr.strip() or 'failed'})  {_bad('✗')}")
    except Exception as exc:
        print(f"  bcftools: ERROR ({exc})  {_bad('✗')}")

    tool_version_cmd = params.get("tool_version_command")
    if tool_version_cmd:
        try:
            res = subprocess.run(tool_version_cmd, shell=True, capture_output=True, text=True)
            if res.returncode == 0:
                print(f"  annotation tool: {res.stdout.strip()}  {_ok('✓')}")
            else:
                print(f"  annotation tool: ERROR ({res.stderr.strip() or 'failed'})  {_bad('✗')}")
        except Exception as exc:
            print(f"  annotation tool: ERROR ({exc})  {_bad('✗')}")
    else:
        print(f"  annotation tool: {_bad('no tool_version_command provided')}  {_bad('✗')}")

    def _substitute_params(text: str, params_map: dict) -> str:
        def _replace(match: re.Match) -> str:
            key = match.group(1)
            parts = key.split(".")
            cur = params_map
            for part in parts:
                if not isinstance(cur, dict) or part not in cur:
                    return match.group(0)
                cur = cur[part]
            return str(cur)
        return re.sub(r"\$\{params\.([A-Za-z0-9_\\.]+)\}", _replace, text)

    substituted = _substitute_params(command, params)
    print("\n" + _hdr("Annotation command (with params substituted)"))
    print(substituted)


def _find_cache_dir(path_hint: Path) -> Path:
    """Resolve various user inputs to the cache directory.

    Accepts either the cache root, the cache directory itself, or a specific
    annotation directory (e.g., /cache/db/cache/vep_gnomad). Returns the path to
    the cache directory that contains annotation subfolders.
    """

    if (path_hint / "cache").exists():
        return path_hint / "cache"

    if path_hint.name == "cache" and path_hint.exists():
        return path_hint

    annotation_dir = path_hint
    if (annotation_dir / "vcfcache_annotated.bcf").exists():
        return annotation_dir.parent

    raise FileNotFoundError(
        "Could not locate a cache directory. Provide -a pointing to a cache root, "
        "cache directory, or an annotation directory containing vcfcache_annotated.bcf."
    )


def _list_annotation_caches(path_hint: Path) -> list[str]:
    """Return sorted annotation cache names under the given path hint.

    Marks incomplete caches (still building) with ' (incomplete)' suffix.
    """
    cache_dir = _find_cache_dir(path_hint)
    names = []
    for child in cache_dir.iterdir():
        if not child.is_dir():
            continue
        # Check if cache is complete (has annotated BCF file index)
        is_complete = (child / "vcfcache_annotated.bcf.csi").exists()
        cache_name = child.name
        if not is_complete:
            cache_name += " (incomplete)"
        names.append(cache_name)
    return sorted(names)


def main() -> None:
    """Main entry point for the vcfcache command-line interface.

    Parses command-line arguments and executes the appropriate command.
    """
    # Load .env file if present (for VCFCACHE_DIR, ZENODO_TOKEN, etc.)
    _load_dotenv()

    parser = argparse.ArgumentParser(
        description="Speed up VCF annotation by using pre-cached common variants.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # Get version, fallback to __init__.py if package not installed
    try:
        version_str = pkg_version("vcfcache")
    except Exception:
        from vcfcache import __version__
        version_str = __version__

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=version_str,
        help="Show version and exit",
    )

    # Create parent parser for shared arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=1,
        help="Increase verbosity: default is INFO, -v for DEBUG",
    )
    parent_parser.add_argument(
        "-q",
        "--quiet",
        action="store_const",
        const=0,
        dest="verbose",
        help="Quiet mode: only show warnings and errors",
    )
    parent_parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Debug mode, keeping intermediate files such as the work directory",
    )
    # Define params in parent parser but don't set required
    parent_parser.add_argument(
        "-y",
        "--yaml",
        dest="params",
        required=False,
        help="(optional) Path to a params YAML containing environment variables related to paths and resources. Defaults to cache's params.snapshot.yaml if not provided.",
    )

    subparsers = parser.add_subparsers(
        dest="command", required=True, title="Available commands", metavar="command"
    )

    # Minimal parent parser for blueprint-init (no config/yaml/manifest)
    init_parent_parser = argparse.ArgumentParser(add_help=False)
    init_parent_parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="count",
        default=1,
        help="(optional) Increase verbosity: default is INFO, -v for DEBUG",
    )
    init_parent_parser.add_argument(
        "-q",
        "--quiet",
        action="store_const",
        const=0,
        dest="verbose",
        help="(optional) Quiet mode: only show warnings and errors",
    )
    init_parent_parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="(optional) Keep intermediate work directory for debugging. "
             "Also uses Zenodo sandbox instead of production for list/download operations.",
    )

    # init command
    init_parser = subparsers.add_parser(
        "blueprint-init",
        help="Initialize blueprint from VCF or Zenodo",
        parents=[init_parent_parser],
        description=(
            "Initialize a blueprint from either a local VCF/BCF file or by downloading from Zenodo. "
            "When creating from VCF: removes genotypes and INFO fields, splits multiallelic sites. "
            "When downloading from Zenodo: extracts blueprint to specified directory."
        )
    )

    # Create mutually exclusive group for source
    source_group = init_parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "-i",
        "--vcf",
        dest="i",
        metavar="VCF",
        help="Input VCF/BCF file to create blueprint from (must be indexed with .csi)"
    )
    source_group.add_argument(
        "--doi",
        dest="doi",
        metavar="DOI",
        help="Zenodo DOI to download blueprint from (e.g., 10.5281/zenodo.XXXXX)"
    )

    init_parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="./cache",
        metavar="DIR",
        help="(optional) Output directory (default: ./cache)"
    )
    init_parser.add_argument(
        "-y",
        "--yaml",
        dest="params",
        required=False,
        metavar="YAML",
        help="(optional) Params YAML used for local blueprint operations",
    )
    init_parser.add_argument(
        "-n",
        "--normalize",
        dest="normalize",
        action="store_true",
        default=False,
        help="(optional) Split multiallelic variants during blueprint creation",
    )
    init_parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="(optional) Force overwrite if output directory exists"
    )

    # blueprint-extend command
    extend_parser = subparsers.add_parser(
        "blueprint-extend",
        help="Add variants to existing blueprint",
        parents=[init_parent_parser],
        description="Extend an existing blueprint by adding variants from a new VCF/BCF file."
    )
    extend_parser.add_argument(
        "-d",
        "--db",
        dest="db",
        required=True,
        metavar="DIR",
        help="Path to existing blueprint directory"
    )
    extend_parser.add_argument(
        "-i",
        "--vcf",
        dest="i",
        required=True,
        metavar="VCF",
        help="Input VCF/BCF file to add (must be indexed with .csi)"
    )
    extend_parser.add_argument(
        "-n",
        "--normalize",
        dest="normalize",
        action="store_true",
        default=False,
        help="(optional) Split multiallelic variants when extending blueprint",
    )
    extend_parser.add_argument(
        "-y",
        "--yaml",
        dest="params",
        required=False,
        metavar="YAML",
        help="(optional) Params YAML used for local blueprint operations",
    )

    # cache-build command
    cache_build_parser = subparsers.add_parser(
        "cache-build",
        help="Build or download annotated cache",
        parents=[init_parent_parser],
        description=(
            "Build an annotated cache from a blueprint, or download a pre-built cache from Zenodo. "
            "\n\n"
            "Two modes:\n"
            "1. Build from blueprint (local or Zenodo): Requires -a/--anno-config to define annotation workflow.\n"
            "2. Download pre-built cache from Zenodo: DOI points to cache, -a forbidden, -n optional."
        )
    )
    # Source: local blueprint directory OR Zenodo DOI (blueprint or cache)
    source = cache_build_parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "-d",
        "--db",
        dest="db",
        metavar="DIR",
        help="Path to existing blueprint directory (requires -a)"
    )
    source.add_argument(
        "--doi",
        dest="doi",
        metavar="DOI",
        help="Zenodo DOI (blueprint or cache). If blueprint: requires -a. If cache: forbids -a."
    )
    cache_build_parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=False,
        metavar="DIR",
        help=(
            "(optional) Base directory for downloaded caches/blueprints when using --doi "
            "(overrides VCFCACHE_DIR for this command)."
        ),
    )
    cache_build_parser.add_argument(
        "-n",
        "--name",
        dest="name",
        required=False,
        metavar="NAME",
        help="(Optional) Name for the cache. Required when building from blueprint. Ignored when downloading pre-built cache."
    )
    cache_build_parser.add_argument(
        "-a",
        "--anno-config",
        dest="anno_config",
        required=False,
        metavar="YAML",
        help="(Optional) Annotation config YAML. Required when building cache from blueprint."
    )
    cache_build_parser.add_argument(
        "-y",
        "--params",
        dest="params",
        required=False,
        metavar="YAML",
        help="(optional) Params YAML file with tool paths and resources. Auto-generated if not provided."
    )
    cache_build_parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="(optional) Force overwrite if cache already exists"
    )
    # Main functionality, apply to user vcf
    vcf_parser = subparsers.add_parser(
        "annotate",
        help="Annotate VCF using pre-built cache",
        parents=[parent_parser],
        description=(
            "Annotate a sample VCF file using a pre-built annotation cache. "
            "The cache enables rapid annotation by reusing annotations from common variants "
            "and only annotating novel variants not found in the cache."
        ),
    )
    vcf_parser.add_argument(
        "-a",
        "--annotation_db",
        dest="a",
        required=True,
        metavar="DIR",
        help=(
            "Path to specific annotation cache directory (e.g., cache_root/cache/vep_gnomad). "
            "Use --list to see available caches."
        ),
    )
    vcf_parser.add_argument(
        "-i",
        "--vcf",
        dest="i",
        required=False,
        metavar="VCF",
        help="Input VCF/BCF file to annotate (required unless using --list or --requirements)",
    )
    vcf_parser.add_argument(
        "-o",
        "--output",
        "--output-file",
        dest="output_file",
        required=False,
        metavar="FILE",
        help=(
            "Output BCF file (required unless using --list or --requirements). "
            "Use '-' or 'stdout' to stream to stdout. "
            "If no extension is provided, '.bcf' is appended."
        ),
    )
    stats_group = vcf_parser.add_mutually_exclusive_group()
    stats_group.add_argument(
        "--stats-dir",
        dest="stats_dir",
        required=False,
        metavar="DIR",
        help=(
            "(optional) Directory to store annotation logs, workflow files, and auxiliary outputs. "
            "If provided, files are written under <stats_dir>/<input_basename>_vcstats. "
            "If omitted, stats are written to <cwd>/<input_basename>_vcstats."
        ),
    )
    stats_group.add_argument(
        "--no-stats",
        action="store_true",
        default=False,
        help="(optional) Disable stats/logs output (disables vcfcache compare).",
    )
    vcf_parser.add_argument(
        "--md5-all",
        action="store_true",
        default=False,
        help=(
            "(optional) Compute MD5 of all variants (no header) and store in stats. "
            "WARNING: can be slow for large files and may differ between runs due to upstream tool quirks."
        ),
    )
    vcf_parser.add_argument(
        "--uncached",
        action="store_true",
        default=False,
        help="(optional) Skip cache, annotate all variants directly. For benchmarking only (default: False)",
    )
    vcf_parser.add_argument(
        "--preserve-unannotated",
        action="store_true",
        default=False,
        help="(optional) Preserve variants without annotation in output. By default, vcfcache mirrors annotation tool behavior (default: False)",
    )
    vcf_parser.add_argument(
        "--skip-split-multiallelic",
        action="store_true",
        default=False,
        help="(optional) Skip splitting multiallelic variants. Use ONLY if certain input has no multiallelic variants (ALT field contains no commas). Provides small speedup (~6%% of runtime) but may cause format inconsistencies if multiallelic variants are present (default: False)",
    )
    vcf_parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="(optional) Force overwrite if output file or stats directory exists (default: False)",
    )
    vcf_parser.add_argument(
        "-p",
        "--parquet",
        dest="parquet",
        action="store_true",
        default=False,
        help="(optional) Also convert output to Parquet format for DuckDB access (default: False)",
    )
    vcf_parser.add_argument(
        "-r",
        "--requirements",
        action="store_true",
        default=False,
        help=(
            "(optional) Show cache requirements (tool, versions, genome build) "
            "and the stored annotation command, then exit."
        ),
    )
    vcf_parser.add_argument(
        "--list",
        action="store_true",
        default=False,
        help=(
            "(optional) List available annotation caches in the specified directory. "
            "Incomplete caches are marked."
        ),
    )

    list_parser = subparsers.add_parser(
        "list",
        help="List or inspect blueprints and caches",
        parents=[init_parent_parser],
        description=(
            "List available vcfcache blueprints and caches from Zenodo or local storage. "
            "Can also inspect individual items for detailed information."
        ),
    )
    list_parser.add_argument(
        "selector",
        choices=["blueprints", "caches"],
        help="What to list: 'blueprints' or 'caches'",
    )
    list_parser.add_argument(
        "--genome",
        metavar="GENOME",
        help="(optional) Filter by genome build (e.g., GRCh38, GRCh37)",
    )
    list_parser.add_argument(
        "--source",
        metavar="SOURCE",
        help="(optional) Filter by data source (e.g., gnomad)",
    )
    list_parser.add_argument(
        "--local",
        nargs="?",
        const="",
        default=None,
        metavar="PATH",
        help=(
            "(optional) List locally available items instead of querying Zenodo. "
            "Optionally specify a custom path (default: VCFCACHE_DIR or ~/.cache/vcfcache)."
        ),
    )
    list_parser.add_argument(
        "--inspect",
        metavar="PATH",
        help=(
            "(optional) Inspect a specific local cache/blueprint. "
            "Accepts a filesystem path or a name under VCFCACHE_DIR/{caches,blueprints}."
        ),
    )

    # push command
    push_parser = subparsers.add_parser(
        "push",
        help="Upload cache to remote storage",
        parents=[init_parent_parser],
        description=(
            "Upload a cache directory to remote storage as a versioned, citable dataset. "
            "Auto-detects blueprint vs cache and generates appropriate naming: "
            "bp_{name}.tar.gz for blueprints, cache_{name}.tar.gz for caches. "
            "Requires ZENODO_TOKEN environment variable (or ZENODO_SANDBOX_TOKEN for --test mode)."
        )
    )
    push_parser.add_argument(
        "--cache-dir",
        required=True,
        metavar="DIR",
        help=(
            "Path to a vcfcache base directory (blueprint upload) or a specific cache "
            "directory under <base>/cache/<cache_name> (cache upload)."
        ),
    )
    push_parser.add_argument(
        "--dest",
        choices=["zenodo"],
        default="zenodo",
        metavar="DEST",
        help="(optional) Upload destination: zenodo (default: zenodo)"
    )
    push_parser.add_argument(
        "--test",
        action="store_true",
        help=(
            "(optional) Upload to test/sandbox environment instead of production. "
            "Uses ZENODO_SANDBOX_TOKEN instead of ZENODO_TOKEN. "
            "Test uploads do not affect production and can be safely deleted."
        )
    )
    push_parser.add_argument(
        "--metadata",
        required=False,
        metavar="FILE",
        help=(
            "(optional) Path to YAML/JSON file with Zenodo metadata for the upload. "
            "If not provided, minimal metadata is auto-generated when using --publish."
        ),
    )
    push_parser.add_argument(
        "--yes",
        action="store_true",
        help="(optional) Skip confirmation prompt and proceed with upload.",
    )
    push_parser.add_argument(
        "--publish",
        action="store_true",
        help=(
            "(optional) Publish the dataset immediately after upload. "
            "If not set, upload will remain as a draft for manual review. "
            "WARNING: Published datasets cannot be deleted, only versioned."
        )
    )

    # demo command
    demo_parser = subparsers.add_parser(
        "demo",
        help="Run comprehensive smoke test of all vcfcache commands",
        parents=[parent_parser],
        description=(
            "Run vcfcache smoke test to verify installation.\n\n"
            "This command tests all 4 main commands (blueprint-init, blueprint-extend,\n"
            "cache-build, annotate) using bundled demo data and verifies that cached\n"
            "and uncached annotation produce identical results.\n\n"
            "For comparing existing annotation runs, use 'vcfcache compare' instead.\n\n"
            "Examples:\n"
            "  vcfcache demo              # Run smoke test\n"
            "  vcfcache demo --debug      # Keep temporary files for inspection\n"
            "  vcfcache demo -q           # Quiet mode (minimal output)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # Note: --debug and -q inherited from parent_parser

    # ===== compare command =====
    compare_parser = subparsers.add_parser(
        "compare",
        help="Compare two vcfcache annotate runs (e.g., cached vs uncached)",
        description=(
            "Compare two successful vcfcache annotate runs and display performance metrics.\n\n"
            "This command requires completion flags (.vcfcache_complete) in both stats directories,\n"
            "which are created by vcfcache annotate.\n\n"
            "Examples:\n"
            "  vcfcache compare run1/uncached_stats run1/cached_stats\n"
            "  vcfcache compare /path/to/stats1 /path/to/stats2"
        ),
        parents=[init_parent_parser],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    compare_parser.add_argument(
        "dir1",
        type=str,
        help="First annotate stats directory",
    )
    compare_parser.add_argument(
        "dir2",
        type=str,
        help="Second annotate stats directory",
    )

    args = parser.parse_args(args=None if sys.argv[1:] else ["--help"])

    show_command_only = args.command == "annotate" and getattr(
        args, "requirements", False
    )
    list_only = args.command == "annotate" and getattr(args, "list", False)

    if show_command_only and list_only:
        parser.error("--requirements and --list cannot be used together")

    if args.command == "annotate" and not (show_command_only or list_only):
        if not args.i or not args.output_file:
            parser.error(
                "annotate command requires -i/--vcf and -o/--output unless --requirements is used"
            )

    # Setup logging with verbosity
    logger = setup_logging(args.verbose)
    log_command(logger)

    # Check bcftools once early (skip for pure manifest ops)
    bcftools_path = None
    if not (show_command_only or list_only or args.command in ["list", "push"]):
        from vcfcache.utils.validation import MIN_BCFTOOLS_VERSION
        logger.debug(f"Minimum required bcftools version: {MIN_BCFTOOLS_VERSION}")
        bcftools_path = check_bcftools_installed()

    try:
        if args.command == "blueprint-init":
            if args.doi:
                # Download blueprint from Zenodo
                zenodo_env = "sandbox" if args.debug else "production"
                logger.info(f"Downloading blueprint from Zenodo ({zenodo_env}) DOI: {args.doi}")
                if args.output == "./cache":
                    base = Path(os.environ.get("VCFCACHE_DIR", "~/.cache/vcfcache")).expanduser()
                    output_dir = (base / "blueprints").resolve()
                else:
                    output_dir = Path(args.output).expanduser().resolve()
                output_dir.mkdir(parents=True, exist_ok=True)

                # Download to temporary tarball
                import tempfile
                with tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False) as tmp:
                    tar_path = Path(tmp.name)

                # Use sandbox if --debug is provided
                download_doi(args.doi, tar_path, sandbox=args.debug)
                logger.info(f"Downloaded to: {tar_path}")

                # Extract
                extracted = extract_cache(tar_path, output_dir)
                logger.info(f"Blueprint extracted to: {extracted}")

                # Clean up tarball
                tar_path.unlink()

            else:
                # Create blueprint from VCF
                logger.debug(f"Creating blueprint from VCF: {Path(args.output)}")

                initializer = DatabaseInitializer(
                    input_file=Path(args.i),
                    params_file=Path(args.params) if getattr(args, "params", None) else None,
                    output_dir=Path(args.output),
                    verbosity=args.verbose,
                    force=args.force,
                    debug=args.debug,
                    bcftools_path=bcftools_path,
                    normalize=args.normalize,
                )
                initializer.initialize()

                # Show detailed timing if --debug is enabled
                if args.debug:
                    workflow_log = Path(args.output) / "blueprint" / "workflow.log"
                    _show_detailed_timings(workflow_log)

        elif args.command == "blueprint-extend":
            logger.debug(f"Adding to blueprint: {args.db}")
            updater = DatabaseUpdater(
                db_path=args.db,
                input_file=args.i,
                params_file=Path(args.params) if getattr(args, "params", None) else None,
                verbosity=args.verbose,
                debug=args.debug,
                bcftools_path=bcftools_path,
                normalize=args.normalize,
            )
            updater.add()

            # Show detailed timing if --debug is enabled
            if args.debug:
                workflow_log = Path(args.db) / "blueprint" / "workflow.log"
                _show_detailed_timings(workflow_log)

        elif args.command == "cache-build":
            # Helper to detect if directory is blueprint or cache
            def is_blueprint(directory: Path) -> bool:
                """Return True if directory is a blueprint, False if cache.

                A cache has both blueprint/ and cache/ with annotation subdirectories.
                A blueprint has only blueprint/ (cache/ is empty or absent).
                """
                blueprint_marker = directory / "blueprint" / "vcfcache.bcf"
                cache_dir = directory / "cache"

                has_blueprint = blueprint_marker.exists()
                has_cache = cache_dir.exists() and cache_dir.is_dir() and any(cache_dir.iterdir())

                # If has both blueprint and cache → it's an annotated cache
                if has_blueprint and has_cache:
                    return False

                # If has only blueprint → it's a blueprint
                if has_blueprint:
                    return True

                # If has only cache → it's a cache
                if has_cache:
                    return False

                # If neither, assume blueprint (for error messaging)
                return True

            # Handle source: local directory or Zenodo DOI
            skip_annotation = False  # Initialize - may be set to True for existing caches

            if args.doi:
                zenodo_env = "sandbox" if args.debug else "production"
                logger.info(f"Downloading from Zenodo ({zenodo_env}) DOI: {args.doi}")

                # Download to appropriate cache directory
                import tempfile
                with tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False) as tmp:
                    tar_path = Path(tmp.name)

                download_doi(args.doi, tar_path, sandbox=args.debug)

                # Determine cache root (CLI flag > env var > default)
                if getattr(args, "output", None):
                    cache_base = Path(args.output).expanduser()
                    logger.info(f"Using download output directory: {cache_base}")
                else:
                    vcfcache_root = os.environ.get("VCFCACHE_DIR")
                    if vcfcache_root:
                        cache_base = Path(vcfcache_root)
                        logger.info(f"Using VCFCACHE_DIR: {cache_base}")
                    else:
                        cache_base = Path.home() / ".cache/vcfcache"

                # Extract to temporary location to detect type
                temp_extract = cache_base / "temp"
                temp_extract.mkdir(parents=True, exist_ok=True)
                extracted_dir = extract_cache(tar_path, temp_extract)
                tar_path.unlink()

                # Detect if blueprint or cache
                if is_blueprint(extracted_dir):
                    logger.info(f"Downloaded blueprint: {extracted_dir}")

                    # Blueprint: requires -a and -n
                    if not args.anno_config:
                        raise ValueError(
                            "DOI points to a blueprint. -a/--anno-config is required to build cache."
                        )
                    if not args.name:
                        raise ValueError(
                            "DOI points to a blueprint. -n/--name is required to name the cache."
                        )

                    # Move to blueprints cache
                    blueprint_store = cache_base / "blueprints"
                    blueprint_store.mkdir(parents=True, exist_ok=True)
                    final_dir = blueprint_store / extracted_dir.name
                    if final_dir.exists():
                        import shutil
                        shutil.rmtree(final_dir)
                    extracted_dir.rename(final_dir)
                    temp_extract.rmdir()

                    db_path = final_dir
                    is_prebuilt_cache = False
                else:
                    logger.info(f"Downloaded pre-built cache: {extracted_dir}")

                    # Pre-built cache: forbid -a
                    if args.anno_config:
                        raise ValueError(
                            "DOI points to a pre-built cache. -a/--anno-config must not be provided. "
                            "The cache is already annotated and ready to use."
                        )

                    # Move to caches directory
                    cache_store = cache_base / "caches"
                    cache_store.mkdir(parents=True, exist_ok=True)
                    final_dir = cache_store / extracted_dir.name
                    if final_dir.exists():
                        import shutil
                        shutil.rmtree(final_dir)
                    extracted_dir.rename(final_dir)
                    temp_extract.rmdir()

                    logger.info(f"Pre-built cache ready at: {final_dir}")
                    logger.info(
                        f"Use with: vcfcache annotate -a {final_dir} -i sample.vcf -o sample_vc.bcf [--stats-dir output/]"
                    )
                    is_prebuilt_cache = True
            else:
                # Local blueprint/cache directory
                if not args.name:
                    raise ValueError(
                        "Local directory requires -n/--name to specify cache name."
                    )

                db_path = Path(args.db)
                if not args.anno_config:
                    raise ValueError(
                        "Local directory requires -a/--anno-config to build cache."
                    )
                logger.debug(f"Using local blueprint: {db_path}")
                is_prebuilt_cache = False

            # If pre-built cache from Zenodo, we're done
            if is_prebuilt_cache:
                return

            # Build cache from blueprint (unless skipping)
            if not is_prebuilt_cache:
                logger.debug(f"Running annotation workflow on blueprint: {db_path}")

                annotator = DatabaseAnnotator(
                    annotation_name=args.name,
                    db_path=db_path,
                    anno_config_file=Path(args.anno_config),
                    params_file=Path(args.params) if args.params else None,
                    verbosity=args.verbose,
                    force=args.force,
                    debug=args.debug,
                    bcftools_path=bcftools_path,
                )
                annotator.annotate()

            # Show detailed timing if --debug is enabled
            if args.debug:
                workflow_log = Path(db_path) / "cache" / args.name / "workflow.log"
                _show_detailed_timings(workflow_log)

        elif args.command == "annotate":
            if args.requirements:
                params_override = Path(args.params).expanduser().resolve() if getattr(args, "params", None) else None
                _print_annotation_command(Path(args.a), params_override=params_override)
                return

            if args.list:
                names = _list_annotation_caches(Path(args.a) if args.a else Path.cwd())
                if not names:
                    print("No cached annotations found.")
                else:
                    print("Available cached annotations:")
                    for name in names:
                        print(f"- {name}")
                return

            # Always show what we're doing (even in default mode)
            input_name = Path(args.i).name
            mode = "uncached" if args.uncached else "cached"
            print(f"Annotating {input_name} ({mode} mode)...")

            vcf_annotator = VCFAnnotator(
                annotation_db=args.a,
                input_vcf=args.i,
                params_file=Path(args.params) if args.params else None,
                output_file=args.output_file,
                stats_dir=args.stats_dir,
                no_stats=args.no_stats,
                verbosity=args.verbose,
                force=args.force,
                debug=args.debug,
                bcftools_path=bcftools_path,
            )

            preserve_unannotated = getattr(args, 'preserve_unannotated', False)
            skip_split_multiallelic = getattr(args, 'skip_split_multiallelic', False)
            vcf_annotator.annotate(
                uncached=args.uncached,
                convert_parquet=args.parquet,
                preserve_unannotated=preserve_unannotated,
                skip_split_multiallelic=skip_split_multiallelic,
                md5_all=args.md5_all,
            )

            # Show detailed timing if --debug is enabled
            if args.debug and args.stats_dir:
                output_name = Path(args.output_file).name if args.output_file not in ("-", "stdout") else "stdout"
                workflow_log = Path(args.stats_dir) / f"{output_name}_vcstats" / "workflow.log"
                _show_detailed_timings(workflow_log)

        elif args.command == "list":
            # Get selector (required positional argument)
            item_type = args.selector

            def _inspect_local(path_or_alias: str) -> None:
                import re

                def resolve_path(s: str) -> Path:
                    p = Path(s).expanduser()
                    if p.exists():
                        return p.resolve()

                    cache_base = Path(os.environ.get("VCFCACHE_DIR", "~/.cache/vcfcache")).expanduser()
                    for sub in ("caches", "blueprints"):
                        cand = (cache_base / sub / s).expanduser()
                        if cand.exists():
                            return cand.resolve()
                    raise FileNotFoundError(
                        f"Could not find '{s}'. Provide an existing path, or a name under "
                        f"{cache_base}/caches or {cache_base}/blueprints."
                    )

                def is_cache_root(root: Path) -> bool:
                    return (root / "blueprint" / "vcfcache.bcf").exists() and (root / "cache").is_dir()

                def is_annotation_dir(p: Path) -> bool:
                    return (p / "annotation.yaml").exists() and (p / "vcfcache_annotated.bcf").exists()

                def parse_required_params(annotation_text: str) -> list[str]:
                    keys = set(re.findall(r"\$\{params\.([A-Za-z0-9_]+)\}", annotation_text))
                    return sorted(keys)

                target = resolve_path(path_or_alias)

                # Normalize target into either a cache root or a specific annotation dir
                cache_root: Path | None = None
                annotation_dirs: list[Path] = []

                if is_annotation_dir(target):
                    cache_root = target.parent.parent
                    annotation_dirs = [target]
                elif is_cache_root(target):
                    cache_root = target
                    cache_dir = cache_root / "cache"
                    annotation_dirs = [p for p in cache_dir.iterdir() if p.is_dir()]
                elif (target / "blueprint" / "vcfcache.bcf").exists():
                    # Blueprint root with no cache/
                    cache_root = target
                    annotation_dirs = []
                else:
                    raise ValueError(
                        f"Path does not look like a vcfcache cache/blueprint: {target}"
                    )

                print("\nVCFcache inspect")
                print("=" * 80)
                print(f"Path: {cache_root}")

                blueprint_bcf = cache_root / "blueprint" / "vcfcache.bcf"
                if blueprint_bcf.exists():
                    print(f"Blueprint: {blueprint_bcf}")

                if not annotation_dirs:
                    print("Cache: (none found)")
                    print("This looks like a blueprint-only directory.")
                    print("Build a cache with: vcfcache cache-build --db <dir> -a <annotation.yaml> -n <name>")
                    print("=" * 80)
                    return

                if len(annotation_dirs) > 1 and not is_annotation_dir(target):
                    print("Cache annotations:")
                    for p in sorted(annotation_dirs):
                        print(f"- {p.name} ({p})")
                    print("\nInspect a specific annotation dir with:")
                    print(f"  vcfcache list --inspect {cache_root}/cache/<name>")
                    print("=" * 80)
                    return

                anno_dir = annotation_dirs[0]
                anno_yaml = anno_dir / "annotation.yaml"
                params_snapshot = anno_dir / "params.snapshot.yaml"

                print(f"Annotation dir: {anno_dir}")
                if params_snapshot.exists():
                    print(f"Params snapshot: {params_snapshot}")
                else:
                    print("Params snapshot: (missing) params.snapshot.yaml")

                anno_text = anno_yaml.read_text(encoding="utf-8")
                anno_cfg = yaml.safe_load(anno_text) or {}

                required_keys = parse_required_params(anno_text)
                params = {}
                if params_snapshot.exists():
                    params = yaml.safe_load(params_snapshot.read_text(encoding="utf-8")) or {}

                def get_nested(d: dict, dotted: str):
                    cur = d
                    for part in dotted.split("."):
                        if not isinstance(cur, dict) or part not in cur:
                            return None
                        cur = cur[part]
                    return cur

                print("\nAnnotation requirements")
                required_tool_version = anno_cfg.get("required_tool_version")
                must_tag = anno_cfg.get("must_contain_info_tag")
                if required_tool_version:
                    print(f"- required_tool_version: {required_tool_version}")
                if must_tag:
                    print(f"- must_contain_info_tag: {must_tag}")

                print("\nRequired params.yaml keys referenced by annotation.yaml")
                if not required_keys:
                    print("(none)")
                else:
                    for k in required_keys:
                        val = get_nested(params, k)
                        if val is None:
                            print(f"- {k}: MISSING")
                        else:
                            print(f"- {k}: {val}")

                print("\nMinimal params.yaml template")
                if required_keys:
                    print("params:")
                    for k in required_keys:
                        print(f"  {k}: {get_nested(params, k) if get_nested(params, k) is not None else '<fill-me>'}")
                print("=" * 80)

            def _dir_size_mb(root: Path) -> float:
                total = 0
                for p in root.rglob("*"):
                    try:
                        if p.is_file():
                            total += p.stat().st_size
                    except FileNotFoundError:
                        continue
                return total / (1024 * 1024)

            def _is_valid_blueprint_root(root: Path) -> bool:
                """Check if directory is a valid blueprint using .vcfcache_complete."""
                complete_file = root / ".vcfcache_complete"
                if not complete_file.exists():
                    return False
                try:
                    complete_data = yaml.safe_load(complete_file.read_text())
                    return (
                        complete_data.get("completed") is True
                        and complete_data.get("mode") == "blueprint-init"
                    )
                except Exception:
                    return False

            def _is_valid_cache_root(root: Path) -> bool:
                """Check if directory contains valid caches using .vcfcache_complete."""
                cache_dir = root / "cache"
                if not cache_dir.is_dir():
                    return False
                for p in cache_dir.iterdir():
                    if not p.is_dir():
                        continue
                    complete_file = p / ".vcfcache_complete"
                    if not complete_file.exists():
                        continue
                    try:
                        complete_data = yaml.safe_load(complete_file.read_text())
                        if (
                            complete_data.get("completed") is True
                            and complete_data.get("mode") == "cache-build"
                        ):
                            return True
                    except Exception:
                        continue
                return False

            def _get_genome_build(root: Path, item_type: str) -> str | None:
                """Extract genome_build from appropriate YAML file."""
                try:
                    if item_type == "blueprints":
                        # Read from workflow/init.yaml
                        yaml_file = root / "workflow" / "init.yaml"
                    else:
                        # Read from cache/<cache_name>/params.snapshot.yaml
                        cache_dir = root / "cache"
                        if not cache_dir.is_dir():
                            return None
                        # Find first valid cache subdirectory
                        for p in cache_dir.iterdir():
                            if p.is_dir():
                                yaml_file = p / "params.snapshot.yaml"
                                if yaml_file.exists():
                                    break
                        else:
                            return None

                    if yaml_file.exists():
                        data = yaml.safe_load(yaml_file.read_text())
                        return data.get("genome_build")
                except Exception:
                    pass
                return None

            def _matches_filters(root: Path, item_type: str, genome_filter: str | None, source_filter: str | None) -> bool:
                """Check if item matches genome and source filters."""
                from vcfcache.utils.naming import CacheName

                # Apply genome filter
                if genome_filter:
                    genome_build = _get_genome_build(root, item_type)
                    if not genome_build or genome_build.lower() != genome_filter.lower():
                        return False

                # Apply source filter
                if source_filter:
                    try:
                        parsed = CacheName.parse(root.name)
                        if parsed.source.lower() != source_filter.lower():
                            return False
                    except Exception:
                        return False

                return True

            def _cache_relevant_size_mb(root: Path) -> float:
                # Cheap “semantic size” to avoid crawling arbitrary large folders.
                files: list[Path] = []
                files.extend(
                    [
                        root / "blueprint" / "vcfcache.bcf",
                        root / "blueprint" / "vcfcache.bcf.csi",
                        root / "blueprint" / "sources.info",
                    ]
                )
                cache_dir = root / "cache"
                if cache_dir.is_dir():
                    for p in cache_dir.iterdir():
                        if not p.is_dir():
                            continue
                        files.extend(
                            [
                                p / "vcfcache_annotated.bcf",
                                p / "vcfcache_annotated.bcf.csi",
                                p / "annotation.yaml",
                                p / "params.snapshot.yaml",
                            ]
                        )
                total = 0
                for f in files:
                    try:
                        if f.exists() and f.is_file():
                            total += f.stat().st_size
                    except FileNotFoundError:
                        continue
                return total / (1024 * 1024)

            def _blueprint_variant_count(root: Path) -> int | None:
                from vcfcache.utils.validation import find_bcftools
                bcftools = find_bcftools()
                if not bcftools:
                    return None
                bcf = root / "blueprint" / "vcfcache.bcf"
                if not bcf.exists():
                    return None
                try:
                    index_path = bcf.with_suffix(bcf.suffix + ".csi")
                    target = index_path if index_path.exists() else bcf
                    result = subprocess.run(
                        [bcftools, "index", "-n", str(target)],
                        capture_output=True,
                        text=True,
                    )
                    if result.returncode != 0:
                        return None
                    return int(result.stdout.strip())
                except Exception:
                    return None

            def _cache_dir_size_mb(root: Path) -> float:
                cache_dir = root / "cache"
                if cache_dir.is_dir():
                    return _dir_size_mb(cache_dir)
                return 0.0

            def _blueprint_dir_size_mb(root: Path) -> float:
                bp_dir = root / "blueprint"
                if bp_dir.is_dir():
                    return _dir_size_mb(bp_dir)
                return 0.0

            def _cache_variant_count(root: Path) -> int | None:
                from vcfcache.utils.validation import find_bcftools
                bcftools = find_bcftools()
                if not bcftools:
                    return None
                cache_dir = root / "cache"
                if not cache_dir.is_dir():
                    return None
                total = 0
                any_found = False
                for p in cache_dir.iterdir():
                    if not p.is_dir():
                        continue
                    bcf = p / "vcfcache_annotated.bcf"
                    if not bcf.exists():
                        continue
                    try:
                        index_path = bcf.with_suffix(bcf.suffix + ".csi")
                        target = index_path if index_path.exists() else bcf
                        result = subprocess.run(
                            [bcftools, "index", "-n", str(target)],
                            capture_output=True,
                            text=True,
                        )
                        if result.returncode != 0:
                            continue
                        total += int(result.stdout.strip())
                        any_found = True
                    except Exception:
                        continue
                return total if any_found else None
                if not bcf.exists():
                    return None
                try:
                    result = subprocess.run(
                        [bcftools, "index", "-n", str(bcf)],
                        capture_output=True,
                        text=True,
                        check=True,
                    )
                    return int(result.stdout.strip())
                except Exception:
                    return None

            def _list_local(item_type: str, base_dir: str | None, genome_filter: str | None = None, source_filter: str | None = None) -> None:
                default_base = Path(os.environ.get("VCFCACHE_DIR", "~/.cache/vcfcache")).expanduser()
                base = Path(base_dir).expanduser() if base_dir else default_base

                # If the user provided a path, treat it as:
                # - a base dir (contains caches/ and/or blueprints/) -> use <base>/<item_type>
                # - an item dir (already points at caches/ or blueprints/, or directly contains items) -> use it as-is
                search_dir = base

                if base_dir:
                    if base.name in ("caches", "blueprints"):
                        search_dir = base
                    elif (base / "caches").is_dir() or (base / "blueprints").is_dir():
                        search_dir = base / item_type
                    else:
                        # User likely passed the directory that directly contains cache/blueprint roots.
                        search_dir = base
                else:
                    search_dir = default_base / item_type

                if not search_dir.exists():
                    print(f"No local {item_type} found at: {search_dir}")
                    return

                candidates = [p for p in sorted(search_dir.iterdir()) if p.is_dir()]
                entries: list[Path] = []
                for root in candidates:
                    if item_type == "blueprints":
                        if not _is_valid_blueprint_root(root):
                            continue
                    else:
                        if not _is_valid_cache_root(root):
                            continue

                    # Apply genome and source filters
                    if not _matches_filters(root, item_type, genome_filter, source_filter):
                        continue

                    if item_type == "blueprints":
                        size_mb = _blueprint_dir_size_mb(root)
                    else:
                        size_mb = _cache_dir_size_mb(root)
                    if size_mb < 1.0:
                        continue
                    entries.append(root)
                if not entries:
                    filter_msg = ""
                    if genome_filter or source_filter:
                        filters = []
                        if genome_filter:
                            filters.append(f"genome={genome_filter}")
                        if source_filter:
                            filters.append(f"source={source_filter}")
                        filter_msg = f" (with filters: {', '.join(filters)})"
                    print(f"No local {item_type} found at: {search_dir}{filter_msg}")
                    return

                print(f"\nLocal vcfcache {item_type}:")
                print("=" * 80)
                for root in entries:
                    if item_type == "blueprints":
                        size_mb = _blueprint_dir_size_mb(root)
                    else:
                        size_mb = _cache_dir_size_mb(root)
                    mtime = ""
                    try:
                        mtime = __import__("datetime").datetime.fromtimestamp(root.stat().st_mtime).date().isoformat()
                    except Exception:
                        mtime = "Unknown"

                    title = root.name
                    # If the directory name is an alias, format it like Zenodo titles.
                    fake_record = {"title": "", "keywords": [root.name]}
                    title = _display_title(fake_record, "caches" if item_type == "caches" else "blueprints")

                    print(f"\n{title}")
                    if item_type == "blueprints":
                        variants = _blueprint_variant_count(root)
                        variants_str = f"{variants:,}" if variants is not None else "N/A"
                        print(f"  Path: {root} | Updated: {mtime} | Size: {size_mb:.1f} MB | Variants: {variants_str}")
                    else:
                        variants = _cache_variant_count(root)
                        variants_str = f"{variants:,}" if variants is not None else "N/A"
                        print(f"  Path: {root} | Updated: {mtime} | Size: {size_mb:.1f} MB | Variants: {variants_str}")

                print(f"\n{'=' * 80}")
                print(f"Total: {len(entries)} {item_type} found")
                if item_type == "caches":
                    print("Inspect: vcfcache list --inspect <path-to-cache-root-or-annotation-dir>")

            def _pretty_source(s: str) -> str:
                if s.lower() == "gnomad":
                    return "gnomAD"
                return s

            def _pretty_tool(t: str) -> str:
                if t.lower() == "vep":
                    return "VEP"
                return t

            def _format_release(release: str) -> str:
                import re
                m = re.match(r"^(\d+(?:\.\d+)*)(.*)$", release)
                if not m:
                    return release
                version, suffix = m.group(1), m.group(2)
                if suffix:
                    suffix = suffix.lstrip("-")
                    return f"v{version} {suffix}"
                return f"v{version}"

            def _format_af_filter(filt: str) -> str:
                # Convention: AF#### where ####/1000 is allele frequency threshold.
                # Example: AF0100 -> 0.100 (10%)
                import re
                m = re.match(r"^AF(\d+)$", filt)
                if not m:
                    return filt
                af = int(m.group(1)) / 1000.0
                af_str = f"{af:.2f}" if af >= 0.01 else f"{af:.3f}"
                return f"AF \u2265 {af_str}"

            def _display_title(record: dict, item_type: str) -> str:
                from vcfcache.utils.naming import CacheName

                title = record.get("title")
                if isinstance(title, str) and title.strip() and title.strip() != "Unknown":
                    return title.strip()

                keywords = record.get("keywords") or []
                alias = next((k for k in keywords if isinstance(k, str) and (k.startswith("cache-") or k.startswith("bp-"))), None)
                if not alias:
                    return record.get("title", "Unknown")
                try:
                    parsed = CacheName.parse(alias)
                except Exception:
                    return record.get("title", "Unknown")

                source = _pretty_source(parsed.source)
                genome = parsed.genome
                release = _format_release(parsed.release)
                filt = _format_af_filter(parsed.filt)

                if item_type == "caches" and parsed.is_cache:
                    tool = _pretty_tool(parsed.tool or "")
                    tool_version = parsed.tool_version or ""
                    preset = (parsed.preset or "").strip()
                    preset_part = f" {preset}" if preset else ""
                    return f"{source} {release} {genome} {filt} annotated with {tool} {tool_version}{preset_part} annotations".strip()

                # Blueprint (or non-standard record): keep it compact but descriptive.
                return f"{source} {release} {genome} {filt} blueprint".strip()

            # Handle --inspect mode
            if args.inspect:
                _inspect_local(args.inspect)
                return

            # Handle --local mode
            if args.local is not None:
                # --local was provided, with optional path argument
                local_path = args.local if args.local else None
                genome_filter = args.genome if hasattr(args, "genome") else None
                source_filter = args.source if hasattr(args, "source") else None
                _list_local(item_type, local_path, genome_filter, source_filter)
                return

            zenodo_env = "sandbox" if args.debug else "production"
            logger.info(f"Searching Zenodo ({zenodo_env}) for vcfcache {item_type}...")

            # Search Zenodo (use sandbox if --debug is provided)
            records = search_zenodo_records(
                item_type=item_type,
                genome=args.genome if hasattr(args, "genome") else None,
                source=args.source if hasattr(args, "source") else None,
                sandbox=args.debug,
                min_size_mb=1.0,
            )

            if not records:
                zenodo_msg = "Zenodo Sandbox" if args.debug else "Zenodo"
                print(f"No {item_type} found on {zenodo_msg}.")
                return

            # Display results
            print(f"\nAvailable vcfcache {item_type} on Zenodo:")
            print("=" * 80)

            for record in records:
                title = _display_title(record, item_type)
                doi = record.get("doi", "Unknown")
                created = record.get("created", "Unknown")
                size_mb = record.get("size_mb", 0)

                print(f"\n{title}")
                print(f"  DOI: {doi} | Created: {created} | Size: {size_mb:.1f} MB")

            print(f"\n{'=' * 80}")
            print(f"Total: {len(records)} {item_type} found")

            # Show appropriate download instructions based on type
            cache_location = os.environ.get("VCFCACHE_DIR", "~/.cache/vcfcache")
            debug_flag = " --debug" if args.debug else ""
            if item_type == "blueprints":
                print(f"Download: vcfcache blueprint-init --doi <DOI> -o <output_dir>{debug_flag}")
                print(
                    f"Or build cache: vcfcache cache-build --doi <DOI> -a <annotation.yaml> -n <name>{debug_flag}"
                )
                print(f"  (downloads to {cache_location}/blueprints/)\n")
            else:  # caches
                print(f"Download: vcfcache cache-build --doi <DOI>{debug_flag}")
                print(f"  (downloads to {cache_location}/caches/)")
                print(
                    f"Then use: vcfcache annotate -a {cache_location}/caches/<cache_name> -i sample.vcf -o sample_vc.bcf [--stats-dir output/]"
                )
                print(f"\nTip: Set VCFCACHE_DIR=/path/to/large/disk to change download location\n")

        elif args.command == "push":
            from vcfcache.integrations import zenodo
            from vcfcache.utils.archive import file_md5
            import json

            # Use --test flag to determine sandbox mode
            sandbox = args.test
            token = (
                os.environ.get("ZENODO_SANDBOX_TOKEN")
                if sandbox
                else os.environ.get("ZENODO_TOKEN")
            )
            if not token:
                raise RuntimeError(
                    "ZENODO_SANDBOX_TOKEN environment variable required for --test mode"
                    if sandbox
                    else "ZENODO_TOKEN environment variable required for push"
                )

            cache_dir = Path(args.cache_dir).expanduser().resolve()
            if not cache_dir.exists():
                raise FileNotFoundError(f"Cache directory not found: {cache_dir}")

            selected_cache_name = None
            base_dir = cache_dir

            if cache_dir.name == "cache" and cache_dir.is_dir():
                is_cache_root = (cache_dir.parent / "blueprint").is_dir() and not (cache_dir / "blueprint").is_dir()
                if is_cache_root:
                    raise ValueError(
                        f"Cache directory {cache_dir} is the cache root. "
                        "Provide a specific cache at <base>/cache/<cache_name> or the base directory."
                    )

            if cache_dir.parent.name == "cache":
                selected_cache_name = cache_dir.name
                base_dir = cache_dir.parent.parent

            has_blueprint_file = (base_dir / "blueprint" / "vcfcache.bcf").exists()
            has_cache_dir = (base_dir / "cache").is_dir()
            has_cache_content = has_cache_dir and any((base_dir / "cache").iterdir())

            if not has_blueprint_file:
                raise ValueError(
                    f"Base directory {base_dir} is missing blueprint/vcfcache.bcf."
                )

            is_cache = selected_cache_name is not None
            is_blueprint = not is_cache

            def _assert_complete(path: Path, expected_mode: str) -> None:
                complete_file = path / ".vcfcache_complete"
                if not complete_file.exists():
                    raise ValueError(f"Missing .vcfcache_complete in {path}")
                try:
                    complete_data = yaml.safe_load(complete_file.read_text()) or {}
                except Exception as exc:
                    raise ValueError(f"Invalid .vcfcache_complete in {path}: {exc}") from exc
                if complete_data.get("completed") is not True:
                    raise ValueError(f"Incomplete run in {path} (.vcfcache_complete completed != true)")
                if complete_data.get("mode") != expected_mode:
                    raise ValueError(
                        f"Unexpected .vcfcache_complete mode in {path}: "
                        f"{complete_data.get('mode')} (expected {expected_mode})"
                    )

            _assert_complete(base_dir, "blueprint-init")
            if is_cache:
                cache_path = base_dir / "cache" / selected_cache_name
                if not cache_path.is_dir():
                    raise ValueError(f"Cache directory not found: {cache_path}")
                _assert_complete(cache_path, "cache-build")

            dir_name = base_dir.name
            prefix = "bp" if is_blueprint else "cache"
            tar_name = f"{prefix}_{dir_name}.tar.gz"
            tar_path = base_dir.parent / tar_name

            if is_blueprint and has_cache_content:
                logger.info(
                    "Blueprint upload selected: existing cache contents will be excluded."
                )

            upload_kind = "blueprint-only" if is_blueprint else f"cache + blueprint ({selected_cache_name})"
            print("Upload plan:")
            print(f"  Base directory: {base_dir}")
            print(f"  Upload type: {upload_kind}")
            if is_cache:
                print(f"  Cache included: cache/{selected_cache_name}")
            else:
                print("  Cache included: (none)")
            print(f"  Tarball: {tar_path}")
            if args.metadata:
                print(f"  Metadata file: {Path(args.metadata).expanduser().resolve()}")
            else:
                print("  Metadata file: (none)")

            if not args.yes:
                resp = input("Proceed with upload? [y/N]: ").strip().lower()
                if resp not in ("y", "yes"):
                    print("Upload cancelled.")
                    return

            logger.info(f"Detected {'blueprint' if is_blueprint else 'cache'}: {dir_name}")
            logger.info(f"Creating archive: {tar_name}")

            tar_cache_subset(
                base_dir,
                tar_path,
                include_blueprint=True,
                include_cache_name=selected_cache_name,
                include_empty_cache_dir=True,
            )
            md5 = file_md5(tar_path)

            logger.info(f"Archive MD5: {md5}")

            dep = zenodo.create_deposit(token, sandbox=sandbox)

            metadata = {}
            if args.metadata:
                mpath = Path(args.metadata).expanduser().resolve()
                text = mpath.read_text()
                metadata = (
                    json.loads(text)
                    if text.strip().startswith("{")
                    else yaml.safe_load(text)
                )

            # Always ensure our deposits are discoverable by API search.
            keywords = ["vcfcache", "blueprint" if is_blueprint else "cache", dir_name]
            try:
                from vcfcache.utils.naming import CacheName

                parsed = CacheName.parse(dir_name)
                keywords.extend([parsed.genome, parsed.source, parsed.release, parsed.filt])
                if parsed.tool:
                    keywords.append(parsed.tool)
                if parsed.tool_version:
                    keywords.append(parsed.tool_version)
                if parsed.preset:
                    keywords.append(parsed.preset)
            except Exception:
                pass
            keywords = sorted({k for k in keywords if k})

            if metadata:
                existing = metadata.get("keywords")
                if isinstance(existing, list):
                    merged = sorted({*existing, *keywords})
                    metadata["keywords"] = merged
                elif existing is None:
                    metadata["keywords"] = keywords

            if args.publish and not metadata:
                # Zenodo requires minimal metadata before publishing.
                item_type = "blueprint" if is_blueprint else "annotated cache"
                metadata = {
                    "title": f"VCFcache {item_type}: {dir_name}",
                    "upload_type": "dataset",
                    "description": (
                        f"VCFcache {item_type} uploaded as {tar_name}. "
                        f"{'This is a test/sandbox record.' if sandbox else ''}"
                    ),
                    "creators": [{"name": "vcfcache"}],
                    "keywords": keywords,
                }

            if metadata:
                zenodo_url = (
                    f"{zenodo._api_base(sandbox)}/deposit/depositions/{dep['id']}"
                )
                resp = requests.put(
                    zenodo_url,
                    params={"access_token": token},
                    json={"metadata": metadata},
                    timeout=30,
                )
                if not resp.ok:
                    error_msg = f"Failed to update metadata: {resp.status_code} {resp.reason}"
                    try:
                        error_detail = resp.json()
                        error_msg += f"\nZenodo error: {error_detail}"
                    except Exception:
                        error_msg += f"\nResponse: {resp.text[:500]}"
                    raise RuntimeError(error_msg)
                resp.raise_for_status()

            zenodo.upload_file(dep, tar_path, token, sandbox=sandbox)
            if args.publish:
                dep = zenodo.publish_deposit(dep, token, sandbox=sandbox)
            print(
                f"Upload complete. Deposition ID: {dep.get('id', 'unknown')} "
                f"DOI: {dep.get('doi', 'draft')} MD5: {md5}"
            )

        elif args.command == "demo":
            from vcfcache.demo import run_smoke_test

            # Run smoke test (only mode now)
            # Derive quiet mode from verbosity: quiet when verbosity is 0 (via -q flag)
            quiet = (args.verbose == 0)
            exit_code = run_smoke_test(keep_files=args.debug, quiet=quiet)
            sys.exit(exit_code)

        elif args.command == "compare":
            from vcfcache.compare import compare_runs

            try:
                compare_runs(Path(args.dir1), Path(args.dir2))
            except (FileNotFoundError, ValueError) as e:
                print(f"Error: {e}")
                sys.exit(1)

    except ZenodoError as e:
        logger.error(str(e))
        raise SystemExit(2)
    except KeyboardInterrupt:
        logger.error("Interrupted.")
        raise SystemExit(130)
    except Exception as e:
        # Keep traceback for unexpected failures (useful for bug reports).
        logger.error(f"Error during execution: {e}")
        raise


if __name__ == "__main__":
    main()
