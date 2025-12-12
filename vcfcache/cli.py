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


def _print_annotation_command(annotation_dir: Path) -> None:
    """Print the stored annotation_tool_cmd from an annotation cache.

    Args:
        annotation_dir: Path to the cache/<annotation_name> directory.
    """

    params_file = annotation_dir / "annotation.yaml"
    if not params_file.exists():
        raise FileNotFoundError(
            f"Params file not found in annotation cache: {params_file}"
        )

    params = yaml.safe_load(params_file.read_text()) or {}
    command = params.get("annotation_tool_cmd")

    if not command:
        raise ValueError(
            "annotation_tool_cmd not found in annotation.yaml; cache may be incomplete"
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
    """Return sorted annotation cache names under the given path hint."""

    cache_dir = _find_cache_dir(path_hint)
    names = []
    for child in cache_dir.iterdir():
        if not child.is_dir():
            continue
        if (child / "vcfcache_annotated.bcf").exists():
            names.append(child.name)
    return sorted(names)


def main() -> None:
    """Main entry point for the vcfcache command-line interface.

    Parses command-line arguments and executes the appropriate command.
    """
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
        help="Annotate blueprint to create cache",
        parents=[init_parent_parser],
        description=(
            "Run annotation workflow on a blueprint to create an annotated cache. "
            "Requires an annotation configuration YAML file that defines the annotation "
            "command (e.g., VEP, SnpEff, bcftools +split-vep) to apply to the blueprint. "
            "The resulting cache can be used for rapid sample annotation."
        )
    )
    cache_build_parser.add_argument(
        "-d",
        "--db",
        dest="db",
        required=True,
        metavar="DIR",
        help="Path to existing blueprint directory"
    )
    cache_build_parser.add_argument(
        "-n",
        "--name",
        dest="name",
        required=True,
        metavar="NAME",
        help="Name for the cache (will create cache/{name}/ subdirectory)"
    )
    cache_build_parser.add_argument(
        "-a",
        "--anno-config",
        dest="anno_config",
        required=True,
        metavar="YAML",
        help="Annotation config YAML file (defines annotation_cmd and must_contain_info_tag)"
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
        "annotate", help="Run a cached VCF annotation", parents=[parent_parser]
    )
    vcf_parser.add_argument(
        "-a",
        "--annotation_db",
        dest="a",
        required=True,
        help="Path to the annotation database directory",
    )
    vcf_parser.add_argument(
        "-i", "--vcf", dest="i", required=False, help="Input VCF file to annotate"
    )
    vcf_parser.add_argument("-o", "--output", required=False, help="Output directory")
    vcf_parser.add_argument(
        "--uncached",
        action="store_true",
        help="Do not use the database (makes only sense for benchmarking)",
    )
    vcf_parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        help="Force overwrite of existing annotation directory",
        default=False,
    )
    vcf_parser.add_argument(
        "-p",
        "--parquet",
        dest="parquet",
        action="store_true",
        help="Convert the final bcf file to parquet format optimized for duck.db access",
        default=False,
    )
    vcf_parser.add_argument(
        "--show-command",
        action="store_true",
        help=(
            "Show the annotation_tool_cmd recorded in the annotation cache and exit. "
            "Skips running any annotation jobs."
        ),
    )
    vcf_parser.add_argument(
        "--list",
        action="store_true",
        help=(
            "List available cached annotation names. Provide -a pointing to a cache root, "
            "a cache directory, or an annotation directory."
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
            logger.debug(f"Running annotation workflow on blueprint: {args.db}")

            annotator = DatabaseAnnotator(
                annotation_name=args.name,
                db_path=args.db,
                anno_config_file=Path(args.anno_config),
                params_file=Path(args.params) if args.params else None,
                verbosity=args.verbose,
                force=args.force,
                debug=args.debug,
                bcftools_path=bcftools_path,
            )
            annotator.annotate()

        elif args.command == "annotate":
            # If annotation_db isn't a path, treat it as a Zenodo alias or DOI.
            alias_or_path = Path(args.a)
            if not alias_or_path.exists():
                sandbox = os.environ.get("ZENODO_SANDBOX", "0") == "1"
                doi, resolved_alias = resolve_zenodo_alias(
                    args.a, item_type="caches", sandbox=sandbox
                )

                cache_store = Path.home() / ".cache/vcfcache/caches"
                cache_store.mkdir(parents=True, exist_ok=True)
                tar_dest = cache_store / f"{resolved_alias}.tar.gz"
                print(
                    f"Downloading cache for alias '{resolved_alias}' from DOI {doi}..."
                )
                download_doi(doi, tar_dest, sandbox=sandbox)
                cache_dir = extract_cache(tar_dest, cache_store)

                # Canonical layout: <root>/cache/<alias>
                candidate = cache_dir / "cache" / resolved_alias
                if candidate.exists():
                    alias_or_path = candidate
                elif (cache_dir / "annotation.yaml").exists():
                    alias_or_path = cache_dir
                elif (cache_dir / "cache").exists():
                    subdirs = [p for p in (cache_dir / "cache").iterdir() if p.is_dir()]
                    if len(subdirs) == 1:
                        alias_or_path = subdirs[0]
                    else:
                        raise FileNotFoundError(
                            f"Could not locate extracted annotation cache for alias '{resolved_alias}' "
                            f"under {cache_dir}"
                        )
                else:
                    raise FileNotFoundError(
                        f"Could not locate extracted annotation cache for alias '{resolved_alias}' "
                        f"under {cache_dir}"
                    )

                args.a = str(alias_or_path)

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

        elif args.command == "list":
            item_type = args.item_type
            zenodo_env = "sandbox" if args.debug else "production"
            logger.info(f"Searching Zenodo ({zenodo_env}) for vcfcache {item_type}...")

            # Search Zenodo (use sandbox if --debug is provided)
            records = search_zenodo_records(
                item_type=item_type,
                genome=args.genome if hasattr(args, "genome") else None,
                source=args.source if hasattr(args, "source") else None,
                sandbox=args.debug,
            )

            if not records:
                zenodo_msg = "Zenodo Sandbox" if args.debug else "Zenodo"
                print(f"No {item_type} found on {zenodo_msg}.")
                return

            # Display results
            print(f"\nAvailable vcfcache {item_type} on Zenodo:")
            print("=" * 80)

            for record in records:
                title = record.get("title", "Unknown")
                doi = record.get("doi", "Unknown")
                created = record.get("created", "Unknown")
                size_mb = record.get("size_mb", 0)

                print(f"\n{title}")
                print(f"  DOI: {doi} | Created: {created} | Size: {size_mb:.1f} MB")

            print(f"\n{'=' * 80}")
            print(f"Total: {len(records)} {item_type} found")
            print(f"Download: vcfcache blueprint-init --doi <DOI> -o ./blueprints\n")

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
                requests.put(
                    zenodo_url,
                    params={"access_token": token},
                    json={"metadata": metadata},
                    timeout=30,
                ).raise_for_status()

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
