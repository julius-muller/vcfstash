
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
import sys

from src.database.initializer import DatabaseInitializer
from src.database.updater import DatabaseUpdater
from src.database.annotator import DatabaseAnnotator, VCFAnnotator
from src.utils.validation import check_bcftools_installed

def main() -> None:
    check_bcftools_installed()
    parser = argparse.ArgumentParser(
        description="Manage VEP database with BCFtools."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # init command
    init_parser = subparsers.add_parser("stash-init", help="Initialize VEP database")
    init_parser.add_argument("-n", "--name", dest="name", required=True, help="Name of the new database directory")
    init_parser.add_argument("-i", "--vcf", dest="i", help="CSI-indexed BCF file")
    init_parser.add_argument("-t", "--threads", dest="t", help="Use multithreading with <int> worker threads [4]", default=4)
    init_parser.add_argument("-f", "--fasta", dest="fasta", required=True, help="Reference FASTA file for variant normalization")
    init_parser.add_argument("-o", "--output", dest="output", default=".", help="Output directory for the database (default: current directory)")

    # add command
    add_parser = subparsers.add_parser("stash-add", help="Add new VCF to the database")
    add_parser.add_argument("-d", "--db", required=True, help="Path to the existing database directory")
    add_parser.add_argument("-i", "--vcf", dest="i", help="Path to the VCF file to be added")
    add_parser.add_argument("-t", "--threads", dest="t", help="Use multithreading with <int> worker threads [4]", default=4)
    add_parser.add_argument("-f", "--fasta", dest="fasta", required=True, help="Reference FASTA file for variant normalization")

    # annotate command
    annotate_parser = subparsers.add_parser("stash-annotate", help="Run annotation workflow on database")
    annotate_parser.add_argument("-d", "--db", required=True, help="Path to the existing database directory")
    annotate_parser.add_argument("-w", "--workflow", required=True, help="Directory containing workflow files (main.nf, nextflow.config)")
    annotate_parser.add_argument("-p", "--params", help="Optional parameters file for nextflow workflow")
    # Add support for passing through Nextflow args
    annotate_parser.add_argument('nextflow_args', nargs=argparse.REMAINDER, help='Additional arguments for Nextflow')

    # Main functionality, apply to user vcf
    vcf_parser = subparsers.add_parser("annotate", help="Annotate a VCF file using the database")
    vcf_parser.add_argument("-d", "--db", required=True, help="Path to the database directory")
    vcf_parser.add_argument("-i", "--vcf", dest="i", required=True, help="Input VCF file to annotate")
    vcf_parser.add_argument("-o", "--output", required=True, help="Output directory")
    vcf_parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")
    # vcf_parser.add_argument("--no-parquet", action="store_true", help="Skip Parquet conversion")

    args = parser.parse_args()

    try:
        if args.command == "stash-init":
            initializer = DatabaseInitializer(
                name=args.name,
                input_file=args.i,
                fasta_ref=args.fasta,
                output_dir=args.output,
                threads=args.t
            )
            initializer.initialize()

        elif args.command == "stash-add":
                    updater = DatabaseUpdater(
                        db_path=args.db,
                        input_file=args.i,
                        fasta_ref=args.fasta,
                        threads=args.t
                    )
                    updater.add()

        elif args.command == "stash-annotate":
            annotator = DatabaseAnnotator(
                db_path=args.db,
                workflow_dir=args.workflow,
                params_file=args.params,
                nextflow_args=args.nextflow_args
            )
            annotator.annotate()

        elif args.command == "annotate":
            annotator = VCFAnnotator(
                db_path=args.db,
                input_vcf=args.i,
                output_dir=args.output,
                threads=args.threads
            )
            annotator.annotate()

    except Exception as e:
        raise  # This will show the full traceback

if __name__ == "__main__":
    main()


# ./vepstash.py stash-init --name test_db --vcf tests/data/nodata/gnomad_test.bcf --fasta tests/data/nodata/reference.fasta --output tests/data/test_out/gtest --threads 2
# ./vepstash.py stash-add --db tests/data/test_out/gtest/test_db/ --vcf tests/data/nodata/dbsnp_test.bcf --fasta tests/data/nodata/reference.fasta --threads 2
# ./vepstash.py stash-annotate --db tests/data/test_out/gtest/test_db --workflow workflow --params config/nextflow_test.yml
# ./vepstash.py annotate --db tests/data/test_out/gtest --vcf tests/data/nodata/sample1.vcf --output tests/data/test_out/gtest/annotated --threads 2