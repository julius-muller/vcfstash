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
from importlib.metadata import version as pkg_version
from pathlib import Path

import requests
import yaml

from vcfcache import EXPECTED_BCFTOOLS_VERSION
from vcfcache.integrations.zenodo import (
    download_doi,
    resolve_zenodo_alias,
    search_zenodo_records,
)
from vcfcache.database.annotator import DatabaseAnnotator, VCFAnnotator
from vcfcache.database.initializer import DatabaseInitializer
from vcfcache.database.updater import DatabaseUpdater
from vcfcache.utils.logging import log_command, setup_logging
from vcfcache.utils.archive import extract_cache, tar_cache
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


def _print_annotation_command(path_hint: Path) -> None:
    """Print the stored annotation_tool_cmd from an annotation cache.

    Args:
        path_hint: Path to cache root, cache directory, or specific annotation directory.
    """
    # Try to find the annotation.yaml file
    # First check if path_hint itself has annotation.yaml (specific cache directory)
    params_file = path_hint / "annotation.yaml"

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
                if not params_file.exists():
                    raise FileNotFoundError(
                        f"Annotation config not found: {params_file} (cache may be incomplete)"
                    )
            else:
                # Multiple caches, ask user to specify
                print(f"Multiple caches found. Please specify which one:")
                for cache in sorted(caches):
                    status = "" if (cache / "vcfcache_annotated.bcf").exists() else " (incomplete)"
                    print(f"  vcfcache annotate --show-command -a {cache}{status}")
                return
        except Exception as e:
            raise FileNotFoundError(
                f"Could not find annotation.yaml in {path_hint}. "
                f"Please provide path to a specific cache directory. Error: {e}"
            )

    params = yaml.safe_load(params_file.read_text()) or {}

    # Try new format (annotation_cmd) first, then fall back to old format (annotation_tool_cmd)
    command = params.get("annotation_cmd") or params.get("annotation_tool_cmd")

    if not command:
        raise ValueError(
            "annotation_cmd or annotation_tool_cmd not found in annotation.yaml; cache may be incomplete"
        )

    print("Annotation command recorded in cache:")
    print(command)


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
        default=0,
        help="Increase verbosity (can be used multiple times, e.g. -vv)",
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
        help="Path to a params YAML containing environment variables related to paths and resources",
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
        default=0,
        help="(optional) Increase verbosity: -v for INFO, -vv for DEBUG",
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
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=1,
        metavar="N",
        help="(optional) Number of threads for bcftools when creating from VCF (default: 1)"
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
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=1,
        metavar="N",
        help="(optional) Number of threads for bcftools (default: 1)"
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
        help="(Optional) Annotation config YAML. Required when source is blueprint. Forbidden when source is pre-built cache."
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
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=1,
        metavar="N",
        help="(optional) Number of threads for bcftools (default: 1)"
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
            "Path to annotation cache directory or cache root. "
            "Use --list to see available caches."
        ),
    )
    vcf_parser.add_argument(
        "-i",
        "--vcf",
        dest="i",
        required=False,
        metavar="VCF",
        help="Input VCF/BCF file to annotate (required unless using --list or --show-command)",
    )
    vcf_parser.add_argument(
        "-o",
        "--output",
        required=False,
        metavar="DIR",
        help="Output directory for annotated VCF (required unless using --list or --show-command)",
    )
    vcf_parser.add_argument(
        "--uncached",
        action="store_true",
        default=False,
        help="(optional) Skip cache, annotate all variants directly. For benchmarking only (default: False)",
    )
    vcf_parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="(optional) Force overwrite if output directory exists (default: False)",
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
        "--show-command",
        action="store_true",
        default=False,
        help=(
            "(optional) Display the annotation command stored in the cache and exit. "
            "Does not run annotation."
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
        help="List available blueprints and caches from Zenodo",
        parents=[init_parent_parser],
        description="Query Zenodo to discover available vcfcache blueprints and caches.",
    )
    list_parser.add_argument(
        "item_type",
        nargs="?",
        choices=["blueprints", "caches"],
        default="blueprints",
        help="Type of items to list: blueprints or caches (default: blueprints)",
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
        "--inspect",
        metavar="PATH_OR_ALIAS",
        help=(
            "(optional) Inspect a local cache/blueprint on disk (no Zenodo query). "
            "Accepts a path, or a directory name under VCFCACHE_DIR/{caches,blueprints}."
        ),
    )
    list_parser.add_argument(
        "--local",
        action="store_true",
        default=False,
        help=(
            "(optional) List locally available caches/blueprints instead of querying Zenodo. "
            "Uses VCFCACHE_DIR by default."
        ),
    )
    list_parser.add_argument(
        "--path",
        metavar="DIR",
        help=(
            "(optional) Base directory to search when using --local. "
            "Defaults to VCFCACHE_DIR or ~/.cache/vcfcache."
        ),
    )

    # Convenience aliases (so `vcfcache caches --local` works)
    caches_parser = subparsers.add_parser(
        "caches",
        help="List caches (Zenodo or local)",
        parents=[init_parent_parser],
        description="Shortcut for: vcfcache list caches",
    )
    caches_parser.add_argument("--genome", metavar="GENOME", help="(optional) Filter by genome build")
    caches_parser.add_argument("--source", metavar="SOURCE", help="(optional) Filter by data source")
    caches_parser.add_argument("--local", action="store_true", default=False, help="(optional) List local caches")
    caches_parser.add_argument("--path", metavar="DIR", help="(optional) Base directory for --local")
    caches_parser.add_argument(
        "--inspect",
        metavar="PATH_OR_ALIAS",
        help="(optional) Inspect a local cache/blueprint on disk (see: vcfcache list --inspect).",
    )

    blueprints_parser = subparsers.add_parser(
        "blueprints",
        help="List blueprints (Zenodo or local)",
        parents=[init_parent_parser],
        description="Shortcut for: vcfcache list blueprints",
    )
    blueprints_parser.add_argument("--genome", metavar="GENOME", help="(optional) Filter by genome build")
    blueprints_parser.add_argument("--source", metavar="SOURCE", help="(optional) Filter by data source")
    blueprints_parser.add_argument("--local", action="store_true", default=False, help="(optional) List local blueprints")
    blueprints_parser.add_argument("--path", metavar="DIR", help="(optional) Base directory for --local")
    blueprints_parser.add_argument(
        "--inspect",
        metavar="PATH_OR_ALIAS",
        help="(optional) Inspect a local cache/blueprint on disk (see: vcfcache list --inspect).",
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
        help="Cache directory to upload (blueprint or annotated cache)"
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
            "(optional) Path to YAML/JSON file with Zenodo metadata. "
            "Should contain: title, description, creators (name, affiliation, orcid), "
            "keywords, upload_type. If not provided, minimal metadata will be auto-generated."
        )
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
        help="Run demo workflow or benchmark cached vs uncached annotation",
        parents=[parent_parser],
    )
    demo_parser.add_argument(
        "--smoke-test",
        action="store_true",
        help="Run comprehensive smoke test of all 4 commands (blueprint-init, blueprint-extend, cache-build, annotate)",
    )
    demo_parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress detailed output (show only essential information)",
    )
    demo_parser.add_argument(
        "-a",
        "--annotation_db",
        type=str,
        help="Path to annotation cache directory (for benchmark mode)",
    )
    demo_parser.add_argument(
        "--vcf",
        type=str,
        help="Path to VCF/BCF file to annotate (for benchmark mode)",
    )
    demo_parser.add_argument(
        "--output",
        type=str,
        help="Output directory for benchmark results (default: temporary directory in /tmp)",
    )
    # Note: -y/--params and --debug inherited from parent_parser

    args = parser.parse_args(args=None if sys.argv[1:] else ["--help"])

    show_command_only = args.command == "annotate" and getattr(
        args, "show_command", False
    )
    list_only = args.command == "annotate" and getattr(args, "list", False)

    if show_command_only and list_only:
        parser.error("--show-command and --list cannot be used together")

    if args.command == "annotate" and not (show_command_only or list_only):
        if not args.i or not args.output:
            parser.error(
                "annotate command requires -i/--vcf and -o/--output unless --show-command is used"
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
                    threads=args.threads,
                    normalize=args.normalize,
                )
                initializer.initialize()

        elif args.command == "blueprint-extend":
            logger.debug(f"Adding to blueprint: {args.db}")
            updater = DatabaseUpdater(
                db_path=args.db,
                input_file=args.i,
                params_file=Path(args.params) if getattr(args, "params", None) else None,
                verbosity=args.verbose,
                debug=args.debug,
                bcftools_path=bcftools_path,
                threads=args.threads,
                normalize=args.normalize,
            )
            updater.add()

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
                    logger.info(f"Use with: vcfcache annotate -a {final_dir} -i sample.vcf -o output/")
                    is_prebuilt_cache = True
            else:
                # Local blueprint directory: always requires -a and -n
                if not args.anno_config:
                    raise ValueError(
                        "Local blueprint directory requires -a/--anno-config to build cache."
                    )
                if not args.name:
                    raise ValueError(
                        "Local blueprint directory requires -n/--name to name the cache."
                    )

                db_path = args.db
                logger.debug(f"Using local blueprint: {db_path}")
                is_prebuilt_cache = False

            # If pre-built cache, we're done
            if is_prebuilt_cache:
                return

            # Build cache from blueprint
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

        elif args.command == "annotate":
            if args.show_command:
                _print_annotation_command(Path(args.a))
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
                output_dir=args.output,
                verbosity=args.verbose,
                force=args.force,
                debug=args.debug,
                bcftools_path=bcftools_path,
            )

            vcf_annotator.annotate(uncached=args.uncached, convert_parquet=args.parquet)

        elif args.command in ("list", "caches", "blueprints"):
            item_type = getattr(args, "item_type", None) or args.command
            if item_type not in ("blueprints", "caches"):
                item_type = "blueprints"

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

            def _list_local(item_type: str, base_dir: str | None) -> None:
                base = Path(base_dir or os.environ.get("VCFCACHE_DIR", "~/.cache/vcfcache")).expanduser()
                search_dir = base
                if base.name != item_type:
                    search_dir = base / item_type

                if not search_dir.exists():
                    print(f"No local {item_type} found at: {search_dir}")
                    return

                entries = [p for p in sorted(search_dir.iterdir()) if p.is_dir()]
                if not entries:
                    print(f"No local {item_type} found at: {search_dir}")
                    return

                print(f"\nLocal vcfcache {item_type}:")
                print("=" * 80)
                for root in entries:
                    size_mb = _dir_size_mb(root)
                    mtime = ""
                    try:
                        mtime = __import__("datetime").datetime.fromtimestamp(root.stat().st_mtime).date().isoformat()
                    except Exception:
                        mtime = "Unknown"

                    title = root.name
                    # If the directory name is an alias, format it like Zenodo titles.
                    fake_record = {"title": root.name, "keywords": [root.name]}
                    title = _display_title(fake_record, "caches" if item_type == "caches" else "blueprints")

                    print(f"\n{title}")
                    print(f"  Path: {root} | Updated: {mtime} | Size: {size_mb:.1f} MB")

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
                from vcfcache.naming import CacheName

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

            if args.inspect:
                _inspect_local(args.inspect)
                return

            if getattr(args, "local", False):
                _list_local(item_type, getattr(args, "path", None))
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
            if item_type == "blueprints":
                print(f"Download: vcfcache blueprint-init --doi <DOI> -o <output_dir>")
                print(f"Or build cache: vcfcache cache-build --doi <DOI> -a <annotation.yaml> -n <name>")
                print(f"  (downloads to {cache_location}/blueprints/)\n")
            else:  # caches
                print(f"Download: vcfcache cache-build --doi <DOI>")
                print(f"  (downloads to {cache_location}/caches/)")
                print(f"Then use: vcfcache annotate -a {cache_location}/caches/<cache_name> -i sample.vcf -o output/")
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

            # Auto-detect blueprint vs cache and generate appropriate name
            # Blueprint: has blueprint/ dir with vcfcache.bcf, cache/ dir is empty or absent
            # Cache: has cache/ dir with subdirectories (named annotation caches)
            has_blueprint_dir = (cache_dir / "blueprint").is_dir()
            has_blueprint_file = (cache_dir / "blueprint" / "vcfcache.bcf").exists()
            cache_subdir = cache_dir / "cache"
            has_cache_content = cache_subdir.is_dir() and any(cache_subdir.iterdir())

            is_blueprint = has_blueprint_dir and has_blueprint_file and not has_cache_content
            is_cache = has_cache_content

            if not is_blueprint and not is_cache:
                raise ValueError(
                    f"Directory {cache_dir} does not appear to be a valid blueprint or cache. "
                    "Expected 'blueprint/vcfcache.bcf' (blueprint) or non-empty 'cache/' directory (cache)."
                )

            dir_name = cache_dir.name
            prefix = "bp" if is_blueprint else "cache"
            tar_name = f"{prefix}_{dir_name}.tar.gz"
            tar_path = cache_dir.parent / tar_name

            logger.info(f"Detected {'blueprint' if is_blueprint else 'cache'}: {dir_name}")
            logger.info(f"Creating archive: {tar_name}")

            tar_cache(cache_dir, tar_path)
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
                from vcfcache.naming import CacheName

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
            from vcfcache.demo import run_smoke_test, run_benchmark

            # Validate mode selection
            if args.smoke_test:
                # Smoke test mode
                exit_code = run_smoke_test(keep_files=args.debug, quiet=args.quiet)
                sys.exit(exit_code)

            elif args.annotation_db or args.vcf:
                # Benchmark mode - validate required arguments
                if not args.annotation_db:
                    print("Error: -a/--annotation_db is required when using --vcf")
                    print("Usage: vcfcache demo -a <cache> --vcf <file> -y <params> [--output <dir>] [--debug]")
                    sys.exit(1)
                if not args.vcf:
                    print("Error: --vcf is required when using -a/--annotation_db")
                    print("Usage: vcfcache demo -a <cache> --vcf <file> -y <params> [--output <dir>] [--debug]")
                    sys.exit(1)
                if not args.params:
                    print("Error: --params (-y) is required for benchmark mode")
                    print("Usage: vcfcache demo -a <cache> --vcf <file> -y <params> [--output <dir>] [--debug]")
                    sys.exit(1)

                # All required arguments provided, run benchmark
                exit_code = run_benchmark(
                    cache_dir=args.annotation_db,
                    vcf_file=args.vcf,
                    params_file=args.params,
                    output_dir=args.output,
                    keep_files=args.debug,
                    quiet=args.quiet,
                )
                sys.exit(exit_code)

            else:
                # No mode selected, show help
                demo_parser.print_help()
                sys.exit(0)

    except Exception as e:
        # Only log the top-level error without traceback - it will be shown by the raise
        logger.error(f"Error during execution: {e}")
        raise  # This will show the full traceback


if __name__ == "__main__":
    main()
