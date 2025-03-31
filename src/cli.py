
"""
VEP Annotation Cache

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
from pathlib import Path
from src.utils.logging import setup_logging, log_command
from src.database.initializer import DatabaseInitializer
from src.database.updater import DatabaseUpdater
from src.database.annotator import DatabaseAnnotator, VCFAnnotator
from src.utils.validation import check_bcftools_installed


def main() -> None:

    parser = argparse.ArgumentParser(
        description="Speed up VEP annotation by using pre-cached common variants", add_help=False
    )

    # Create parent parser for shared arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("-v", "--verbose", action="count", default=0,
                              help="Increase verbosity (can be used multiple times, e.g. -vv)")
    parent_parser.add_argument("-c", "--config", dest="config",
                              required=False, help="Path to a nextflow.config")
    parent_parser.add_argument("--test", action="store_true",
                                 help="Use self contained test cases with reduced resource requirements")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # init command
    init_parser = subparsers.add_parser("stash-init", help="Initialize VEP database", parents=[parent_parser])
    init_parser.add_argument("-n", "--name", dest="name", required=True, help="Name of the new database directory")
    init_parser.add_argument("-i", "--vcf", dest="i", required=True, help="CSI-indexed BCF file")
    init_parser.add_argument("-o", "--output", dest="output", default=".", help="Output directory")
    init_parser.add_argument("-f", "--force", dest="force", action="store_true",
                              help="Force overwrite of existing database directory", default=False)

    # add command
    add_parser = subparsers.add_parser("stash-add", help="Add new VCF to the database", parents=[parent_parser])
    add_parser.add_argument("-d", "--db", required=True, help="Path to the existing database directory")
    add_parser.add_argument("-i", "--vcf", dest="i", help="Path to the VCF file to be added")

    # annotate command
    annotate_parser = subparsers.add_parser("stash-annotate", help="Run annotation workflow on database", parents=[parent_parser])
    annotate_parser.add_argument("-n", "--name", dest="name", required=True,
                             help="Name of the annotation run, will be the directory name")
    annotate_parser.add_argument("-d", "--db", required=True, help="Path to the existing database directory")
    annotate_parser.add_argument("-f", "--force", dest="force", action="store_true",
                             help="Force overwrite of existing annotation directory", default=False)
    annotate_parser.add_argument("-a", "--anno_config", dest="anno_config",
                              required=False, help="Path to a .config file with annotation relevant settings")

    # Main functionality, apply to user vcf
    vcf_parser = subparsers.add_parser("annotate", help="Annotate a VCF file using the database", parents=[parent_parser])
    vcf_parser.add_argument("-a", "--annotation_db", dest='a', required=True, help="Path to the annotation database directory")
    vcf_parser.add_argument("-i", "--vcf", dest="i", required=True, help="Input VCF file to annotate")
    vcf_parser.add_argument("-o", "--output", required=True, help="Output directory")
    vcf_parser.add_argument("--uncached", action="store_true", help="Do not use the database (makes only sense for benchmarking)")
    vcf_parser.add_argument("-f", "--force", dest="force", action="store_true",
                             help="Force overwrite of existing annotation directory", default=False)

    args = parser.parse_args()
    # Setup logging with verbosity
    logger = setup_logging(args.verbose)
    log_command(logger)

    check_bcftools_installed()

    try:
        if args.command == "stash-init":
            logger.info(f"Initializing database: {args.name}")

            initializer = DatabaseInitializer(
                name=args.name,
                input_file=Path(args.i) if args.i else None,
                config_file=Path(args.config) if args.config else None,
                output_dir=Path(args.output),
                verbosity=args.verbose,
                force=args.force,
                test_mode=args.test
            )
            initializer.initialize()

        elif args.command == "stash-add":
            logger.info(f"Adding to database: {args.db}")
            updater = DatabaseUpdater(
                db_path=args.db,
                input_file=args.i,
                config_file=Path(args.config) if args.config else None,
                verbosity=args.verbose,
                test_mode=args.test
            )
            updater.add()

        elif args.command == "stash-annotate":
            logger.info(f"Running annotation workflow on database: {args.db}")

            annotator = DatabaseAnnotator(
                annotation_name = args.name,
                db_path=args.db,
                config_file=Path(args.config) if args.config else None,
                anno_config_file=Path(args.anno_config) if args.anno_config else None,
                verbosity=args.verbose,
                force = args.force,
                test_mode=args.test
            )
            annotator.annotate()

        elif args.command == "annotate":
            logger.info(f"Annotating VCF: {args.i}")

            annotator = VCFAnnotator(
                annotation_db=args.a,
                input_vcf=args.i,
                config_file=Path(args.config) if args.config else None,
                output_dir=args.output,
                verbosity=args.verbose,
                force=args.force,
                uncached=args.uncached
            )

            annotator.annotate()

    except Exception as e:
        logger.error(f"Error during execution: {e}", exc_info=True)
        raise  # This will show the full traceback


def generate_test_command(
        vepstash_path="~/projects/vepstash/vepstash.py",
        stash_name="nftest",
        vcf_path="~/projects/vepstash/tests/data/nodata/crayz_db.bcf",
        output_dir="/home/j380r/tmp/test/test_out",
        config_path="/home/j380r/projects/vepstash/tests/config/nextflow_test.config",
        add_vcf_path="~/projects/vepstash/tests/data/nodata/crayz_db2.bcf",
        annotate_name="testor",
        test_mode=True,
        force=True
):
    """
    Generate a nicely formatted test command string for vepstash operations.

    Returns:
        str: A copy-pastable command string with proper formatting
    """
    # Build each part of the command
    cmd_init = (
        f"{vepstash_path} stash-init "
        f"--name {stash_name} "
        f"--vcf {vcf_path} "
        f"--output {output_dir} "
        f"-c {config_path} "
        f"{'-f' if force else ''} "
        f"{'--test' if test_mode else ''}"
    ).strip()

    cmd_add = (
        f"{vepstash_path} stash-add "
        f"--db {output_dir}/{stash_name} "
        f"-i {add_vcf_path} "
        f"{'--test' if test_mode else ''}"
    ).strip()

    cmd_annotate = (
        f"{vepstash_path} stash-annotate "
        f"--name {annotate_name} "
        f"--db {output_dir}/{stash_name} "
        f"{'--test' if test_mode else ''}"
    ).strip()

    # Combine commands
    full_cmd = f"{cmd_init} ; {cmd_add} ; {cmd_annotate}"

    # Also create a nicely formatted display version for easier reading
    formatted_cmds = f"""
# INITIALIZE
{cmd_init}

# ADD
{cmd_add}

# ANNOTATE
{cmd_annotate}

# COMBINED COMMAND (for copy-paste)
alias stx="{full_cmd}"
"""

    print(formatted_cmds)
    return full_cmd


if __name__ == "__main__":
    main()

    # generate_test_command()


# Example usage into test dir in repo:
# ~/projects/vepstash/vepstash.py stash-init --name nftest --vcf ~/projects/vepstash/tests/data/nodata/crayz_db.bcf --output /home/j380r/tmp/test/test_out -f --test -vv
# ~/projects/vepstash/vepstash.py stash-add --db ~/projects/vepstash/tests/data/test_out/nftest -i ~/projects/vepstash/tests/data/nodata/crayz_db2.bcf --test -vv
# ~/projects/vepstash/vepstash.py stash-annotate --name testor --db ~/projects/vepstash/tests/data/test_out/nftest --test -vv -f
# ... or locally
# ~/projects/vepstash/vepstash.py stash-init --name nftest --vcf ~/projects/vepstash/tests/data/nodata/crayz_db.bcf --output test_out --test -f -vv
# ~/projects/vepstash/vepstash.py stash-add --db test_out/nftest -i ~/projects/vepstash/tests/data/nodata/crayz_db2.bcf --test -vv
# ~/projects/vepstash/vepstash.py stash-annotate --name testor --db test_out/nftest --test -vv -f
# ~/projects/vepstash/vepstash.py annotate --a ~/tmp/test/test_out/nftest/annotations/testor --vcf ~/projects/vepstash/tests/data/nodata/sample4.bcf --output ~/tmp/test/aout --test -f

# as one:
# alias stx="~/projects/vepstash/vepstash.py stash-init --name nftest --vcf ~/projects/vepstash/tests/data/nodata/crayz_db.bcf --output test_out --test -f ; ~/projects/vepstash/vepstash.py stash-add --db test_out/nftest -i ~/projects/vepstash/tests/data/nodata/crayz_db2.bcf --test ; ~/projects/vepstash/vepstash.py stash-annotate --name testor --db test_out/nftest --test "

# ~/projects/vepstash/vepstash.py annotate --a ~/tmp/test/test_out/nftest/annotations/testor --vcf ~/projects/vepstash/tests/data/nodata/sample1.vcf --output ~/tmp/test/aout --test -f