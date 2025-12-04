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
import sys
from importlib.metadata import version as pkg_version
from pathlib import Path

import yaml

from vcfstash import EXPECTED_BCFTOOLS_VERSION
from vcfstash.database.annotator import DatabaseAnnotator, VCFAnnotator
from vcfstash.database.initializer import DatabaseInitializer
from vcfstash.database.updater import DatabaseUpdater
from vcfstash.utils.logging import log_command, setup_logging
from vcfstash.utils.validation import check_bcftools_installed


def _print_annotation_command(annotation_dir: Path) -> None:
    """Print the stored annotation_tool_cmd from an annotation stash.

    Args:
        annotation_dir: Path to the stash/<annotation_name> directory.
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


def _find_stash_dir(path_hint: Path) -> Path:
    """Resolve various user inputs to the stash directory.

    Accepts either the stash root, the stash directory itself, or a specific
    annotation directory (e.g., /cache/db/stash/vep_gnomad). Returns the path to
    the stash directory that contains annotation subfolders.
    """

    if (path_hint / "stash").exists():
        return path_hint / "stash"

    if path_hint.name == "stash" and path_hint.exists():
        return path_hint

    annotation_dir = path_hint
    if (annotation_dir / "vcfstash_annotated.bcf").exists():
        return annotation_dir.parent

    raise FileNotFoundError(
        "Could not locate a stash directory. Provide -a pointing to a stash root, "
        "stash directory, or an annotation directory containing vcfstash_annotated.bcf."
    )


def _list_annotation_caches(path_hint: Path) -> list[str]:
    """Return sorted annotation cache names under the given path hint."""

    stash_dir = _find_stash_dir(path_hint)
    names = []
    for child in stash_dir.iterdir():
        if not child.is_dir():
            continue
        if (child / "vcfstash_annotated.bcf").exists():
            names.append(child.name)
    return sorted(names)


def main() -> None:
    """Main entry point for the vcfstash command-line interface.

    Parses command-line arguments and executes the appropriate command.
    """
    parser = argparse.ArgumentParser(
        description="Speed up VCF annotation by using pre-cached common variants.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=pkg_version("vcfstash"),
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
        help="Debug mode, keeping intermediate files such as the nextflow work directory",
    )
    parent_parser.add_argument(
        "-c",
        "--config",
        dest="config",
        required=False,
        help="Path to an optional nextflow config containing only a process definition",
    )
    # Define params in parent parser but don't set required
    parent_parser.add_argument(
        "-y",
        "--yaml",
        dest="params",
        required=False,
        help="Path to a nextflow yaml containing environment variables related to paths and resources",
    )
    parent_parser.add_argument(
        "--nf",
        "--use-nextflow",
        dest="use_nextflow",
        action="store_true",
        default=False,
        help="Use legacy Nextflow backend (requires Java). Default is pure Python.",
    )

    subparsers = parser.add_subparsers(
        dest="command", required=True, title="Available commands", metavar="command"
    )

    # init command
    init_parser = subparsers.add_parser(
        "stash-init", help="Initialize VCF stash blueprint", parents=[parent_parser]
    )
    init_parser.add_argument(
        "-i", "--vcf", dest="i", required=True, help="CSI-indexed BCF file"
    )
    init_parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default=".",
        help="Output directory, the name of the database will be the top level directory",
    )
    init_parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        help="Force overwrite of existing database directory",
        default=False,
    )
    init_parser.add_argument(
        "-n",
        "--normalize",
        dest="normalize",
        action="store_true",
        help="Apply normalization steps (add chr prefix, filter chromosomes, split multiallelic sites)",
        default=False,
    )

    # add command
    add_parser = subparsers.add_parser(
        "stash-add", help="Add new VCF to the blueprint", parents=[parent_parser]
    )
    add_parser.add_argument(
        "-d", "--db", required=True, help="Path to the existing database directory"
    )
    add_parser.add_argument(
        "-i", "--vcf", dest="i", help="Path to the VCF file to be added", required=True
    )
    add_parser.add_argument(
        "-n",
        "--normalize",
        dest="normalize",
        action="store_true",
        help="Apply normalization steps (add chr prefix, filter chromosomes, split multiallelic sites)",
        default=False,
    )

    # annotate command
    annotate_parser = subparsers.add_parser(
        "stash-annotate",
        help="Run annotation workflow on blueprint and instantiate a stash",
        parents=[parent_parser],
    )
    annotate_parser.add_argument(
        "-n",
        "--name",
        dest="name",
        required=True,
        help="Name of the stash instance, will be the directory name",
    )
    annotate_parser.add_argument(
        "-d", "--db", required=True, help="Path to the existing database directory"
    )
    annotate_parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        help="Force overwrite of existing stash directory",
        default=False,
    )
    annotate_parser.add_argument(
        "-a",
        "--anno_config",
        dest="anno_config",
        required=True,
        help="Path to an nextflow config file containing the annotation commands to run",
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
            "List available cached annotation names. Provide -a pointing to a stash root, "
            "a stash directory, or an annotation directory."
        ),
    )

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

    # Check if required args exists based on command
    if args.command == "stash-init" and not args.params:
        parser.error("stash-init command requires -y/--yaml parameter")

    # Setup logging with verbosity
    logger = setup_logging(args.verbose)
    log_command(logger)

    # Check bcftools if params file is provided (required for stash-init)
    # For other commands, we'll use params from the database or fall back to init.yaml
    bcftools_path = None
    if not (show_command_only or list_only):
        logger.debug(f"Expected bcftools version: {EXPECTED_BCFTOOLS_VERSION}")
        if args.params:
            logger.debug(
                f"Checking bcftools installation using params file: {args.params}"
            )
            bcftools_path = check_bcftools_installed(Path(args.params))
        elif args.command in ["stash-add", "stash-annotate", "annotate"]:
            # For these commands, try to get bcftools path from the workflow directory
            workflow_dir = None
            if args.command == "stash-add" or args.command == "stash-annotate":
                workflow_dir = Path(args.db) / "workflow"
            elif args.command == "annotate":
                workflow_dir = Path(args.a).parent.parent / "workflow"

            if workflow_dir and workflow_dir.exists():
                logger.debug(
                    f"Checking bcftools installation using init.yaml from: {workflow_dir}"
                )
                bcftools_path = check_bcftools_installed(workflow_dir=workflow_dir)
            else:
                logger.warning(f"Workflow directory not found: {workflow_dir}")
                bcftools_path = check_bcftools_installed()

    try:
        if args.command == "stash-init":
            logger.debug(f"Initializing database: {Path(args.output).parent}")

            initializer = DatabaseInitializer(
                input_file=Path(args.i),
                config_file=Path(args.config) if args.config else None,
                params_file=Path(args.params),
                output_dir=Path(args.output),
                verbosity=args.verbose,
                force=args.force,
                debug=args.debug,
                bcftools_path=bcftools_path,
                normalize=args.normalize,
                use_nextflow=args.use_nextflow,
            )
            initializer.initialize()

        elif args.command == "stash-add":
            logger.debug(f"Adding to database: {args.db}")
            updater = DatabaseUpdater(
                db_path=args.db,
                input_file=args.i,
                config_file=Path(args.config) if args.config else None,
                params_file=Path(args.params) if args.params else None,
                verbosity=args.verbose,
                debug=args.debug,
                bcftools_path=bcftools_path,
                normalize=args.normalize,
                use_nextflow=args.use_nextflow,
            )
            updater.add()

        elif args.command == "stash-annotate":
            logger.debug(f"Running annotation workflow on database: {args.db}")

            annotator = DatabaseAnnotator(
                annotation_name=args.name,
                db_path=args.db,
                config_file=Path(args.config) if args.config else None,
                anno_config_file=Path(args.anno_config),
                params_file=Path(args.params) if args.params else None,
                verbosity=args.verbose,
                force=args.force,
                debug=args.debug,
                bcftools_path=bcftools_path,
                use_nextflow=args.use_nextflow,
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
                config_file=Path(args.config) if args.config else None,
                params_file=Path(args.params) if args.params else None,
                output_dir=args.output,
                verbosity=args.verbose,
                force=args.force,
                debug=args.debug,
                bcftools_path=bcftools_path,
                use_nextflow=args.use_nextflow,
            )

            vcf_annotator.annotate(uncached=args.uncached, convert_parquet=args.parquet)

    except Exception as e:
        # Only log the top-level error without traceback - it will be shown by the raise
        logger.error(f"Error during execution: {e}")
        raise  # This will show the full traceback


if __name__ == "__main__":
    main()
