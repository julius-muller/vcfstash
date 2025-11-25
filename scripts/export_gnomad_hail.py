#!/usr/bin/env python3
"""Export gnomAD variants from Hail table to BCF format.

This script queries the gnomAD browser Hail table directly from GCS,
filters by allele frequency, and exports to BCF format. This is much
more efficient than downloading the full VCF files.

Requires:
- hail
- Google Cloud credentials with requester-pays enabled
"""

import argparse
import os
import sys
import urllib.request
from pathlib import Path
import hail as hl


def download_gcs_connector() -> str:
    """Download GCS connector JAR if not already cached.

    Returns:
        Path to the GCS connector JAR file.
    """
    jar_dir = Path.home() / ".vcfstash" / "jars"
    jar_dir.mkdir(parents=True, exist_ok=True)

    # Use version 2.2.30 (newer and more stable)
    jar_name = "gcs-connector-hadoop3-2.2.30-shaded.jar"
    jar_path = jar_dir / jar_name

    if jar_path.exists():
        print(f"Using cached GCS connector: {jar_path}")
        return str(jar_path)

    # Download from Maven Central
    url = "https://repo1.maven.org/maven2/com/google/cloud/bigdataoss/gcs-connector/hadoop3-2.2.30/gcs-connector-hadoop3-2.2.30-shaded.jar"

    print(f"Downloading GCS connector v2.2.30 (~35MB)...")
    print(f"URL: {url}")

    try:
        urllib.request.urlretrieve(url, jar_path)
        print(f"✓ Downloaded to: {jar_path}")
        return str(jar_path)
    except Exception as e:
        print(f"Failed to download shaded connector: {e}")
        print("Trying unshaded version as fallback...")

        # Try unshaded version as fallback
        jar_name = "gcs-connector-hadoop3-2.2.30.jar"
        jar_path = jar_dir / jar_name
        url = "https://repo1.maven.org/maven2/com/google/cloud/bigdataoss/gcs-connector/hadoop3-2.2.30/gcs-connector-hadoop3-2.2.30.jar"

        urllib.request.urlretrieve(url, jar_path)
        print(f"✓ Downloaded unshaded connector to: {jar_path}")
        return str(jar_path)


def export_gnomad_bcf(
    output_path: str,
    af_threshold: float = 0.10,
    genome_build: str = "GRCh38",
    chromosomes: list[str] | None = None,
    gnomad_version: str = "4.1",
    use_public_bucket: bool = True,
    gcs_project: str | None = None,
) -> None:
    """Export gnomAD variants filtered by AF to BCF format.

    Args:
        output_path: Path where the BCF file will be written
        af_threshold: Minimum allele frequency threshold (e.g., 0.10 for 10%)
        genome_build: Genome build (GRCh38 or GRCh37)
        chromosomes: List of chromosomes to include (e.g., ['chr1', 'chr2']).
                     If None, includes all chromosomes.
        gnomad_version: gnomAD version (e.g., "4.1")
        use_public_bucket: Use public bucket (no auth required) vs requester-pays
        gcs_project: GCS project for requester-pays access (only if use_public_bucket=False)
    """
    print(f"Initializing Hail with genome build: {genome_build}")

    # For public buckets, prevent metadata server lookups by setting env var early
    if use_public_bucket:
        os.environ['GCE_METADATA_HOST'] = 'metadata.google.internal.invalid'
        print("Disabled GCE metadata server for anonymous access")

    # Download GCS connector JAR if needed
    gcs_jar = download_gcs_connector()

    # Configure Spark to use GCS connector
    spark_conf = {
        'spark.hadoop.fs.gs.impl': 'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem',
        'spark.hadoop.fs.AbstractFileSystem.gs.impl': 'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS',
        'spark.jars': gcs_jar,  # Use local JAR instead of downloading via Maven
    }

    # Configure authentication based on bucket type
    if use_public_bucket:
        # For public buckets, disable all authentication
        print("Configuring Spark for anonymous GCS access...")
        spark_conf.update({
            'spark.hadoop.google.cloud.auth.service.account.enable': 'false',
            'spark.hadoop.google.cloud.auth.null.enable': 'true',
            'spark.hadoop.fs.gs.auth.service.account.enable': 'false',
            'spark.hadoop.fs.gs.impl': 'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem',
            'spark.hadoop.fs.AbstractFileSystem.gs.impl': 'com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS',
        })
    else:
        # For requester-pays buckets, use application default credentials
        print("Configuring authenticated access for requester-pays bucket...")
        spark_conf['spark.hadoop.google.cloud.auth.service.account.enable'] = 'true'
        if gcs_project:
            spark_conf['spark.hadoop.fs.gs.project.id'] = gcs_project
            # Configure requester pays
            if gcs_project:
                hl.utils.hadoop_copy._configure_requester_pays_gcs(gcs_project)

    # Initialize Hail with appropriate reference genome and GCS support
    hl.init(
        default_reference=genome_build,
        log="/tmp/hail.log",
        spark_conf=spark_conf,
        quiet=False,
    )

    # Choose bucket based on use_public_bucket flag
    if use_public_bucket:
        bucket = "gcp-public-data--gnomad"
        print(f"Using public gnomAD bucket")
    else:
        bucket = "gnomad-public-requester-pays"
        print(f"Using requester-pays bucket")
        if not gcs_project:
            print("WARNING: No GCS project provided for requester-pays bucket")

    # Construct path to gnomAD Hail table (joint WES+WGS dataset)
    ht_path = f"gs://{bucket}/release/{gnomad_version}/ht/joint/gnomad.joint.v{gnomad_version}.sites.ht"

    print(f"Reading gnomAD v{gnomad_version} Hail table from: {ht_path}")
    ht = hl.read_table(ht_path)

    print(f"Filtering variants with AF >= {af_threshold}")
    # Filter for variants where ANY population has AF >= threshold
    # Joint dataset stores frequencies in 'joint.freq' array
    ht_filtered = ht.filter(
        hl.any(lambda freq: freq.AF >= af_threshold, ht.joint.freq)
    )

    # Filter by chromosomes if specified
    if chromosomes:
        print(f"Filtering to chromosomes: {', '.join(chromosomes)}")
        ht_filtered = ht_filtered.filter(
            hl.literal(set(chromosomes)).contains(ht_filtered.locus.contig)
        )

    # Count variants
    n_variants = ht_filtered.count()
    print(f"Found {n_variants:,} variants matching criteria")

    if n_variants == 0:
        print("WARNING: No variants found matching the criteria!")
        sys.exit(1)

    # Select only essential fields for VCF export
    # VCF can't handle complex nested structures, so we keep it simple
    print("Selecting essential fields for VCF export...")

    # Get weighted average AF across all samples from the joint dataset
    # joint.freq[0] contains the combined frequency across all samples (WES+WGS)
    af = hl.or_else(
        ht_filtered.joint.freq[0].AF,
        0.0
    )

    # Get allele count from the first freq element (combined WES+WGS)
    # AC = allele count (number of times this allele was observed)
    allele_count = hl.or_else(
        ht_filtered.joint.freq[0].AC,
        0
    )

    # Keep only locus, alleles, AF, and AC in INFO
    # For VCF export, INFO fields must be in an 'info' struct
    ht_simple = ht_filtered.select(
        info=hl.struct(
            AF=af,
            AC=allele_count
        ),
        filters=hl.empty_set(hl.tstr)  # PASS for all variants
    )

    # Export directly to VCF from Table (sites-only VCF)
    vcf_path = output_path.replace(".bcf", ".vcf.bgz")
    print(f"Exporting to VCF: {vcf_path}")

    # Export the simplified table to VCF format
    hl.export_vcf(ht_simple, vcf_path, tabix=True)

    print(f"✓ Successfully exported {n_variants:,} variants to {vcf_path}")
    print(f"✓ Index file created: {vcf_path}.tbi")

    # Convert to BCF if bcftools is available
    import subprocess
    import shutil

    if shutil.which("bcftools"):
        bcf_path = output_path
        print(f"Converting to BCF: {bcf_path}")

        subprocess.run(
            ["bcftools", "view", "-Ob", "-o", bcf_path, vcf_path],
            check=True
        )
        subprocess.run(
            ["bcftools", "index", bcf_path],
            check=True
        )

        print(f"✓ Successfully converted to BCF: {bcf_path}")
        print(f"✓ BCF index created: {bcf_path}.csi")
    else:
        print("WARNING: bcftools not found, skipping BCF conversion")
        print(f"VCF file available at: {vcf_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Export gnomAD variants from Hail table to BCF format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Export common variants (AF >= 10%) to BCF (no auth required!)
  python export_gnomad_hail.py --output gnomad_af0.10.bcf --af 0.10

  # Export only chromosome 1 with AF >= 5%
  python export_gnomad_hail.py --output gnomad_chr1_af0.05.bcf --af 0.05 --chr chr1

  # Export multiple chromosomes
  python export_gnomad_hail.py --output gnomad_chr1-2.bcf --af 0.10 --chr chr1 chr2

  # Advanced: Use requester-pays bucket (requires auth)
  python export_gnomad_hail.py --output gnomad.bcf --af 0.10 --use-requester-pays --gcs-project my-project
        """
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output BCF file path"
    )

    parser.add_argument(
        "--af",
        type=float,
        default=0.10,
        help="Minimum allele frequency threshold (default: 0.10)"
    )

    parser.add_argument(
        "--genome",
        choices=["GRCh38", "GRCh37"],
        default="GRCh38",
        help="Genome build (default: GRCh38)"
    )

    parser.add_argument(
        "--chr", "--chromosomes",
        nargs="+",
        help="Chromosomes to include (e.g., chr1 chr2). If not specified, includes all."
    )

    parser.add_argument(
        "--gnomad-version",
        default="4.1",
        help="gnomAD version (default: 4.1)"
    )

    parser.add_argument(
        "--use-requester-pays",
        action="store_true",
        help="Use requester-pays bucket instead of public bucket (requires authentication)"
    )

    parser.add_argument(
        "--gcs-project",
        help="GCS project for requester-pays access (only needed with --use-requester-pays)"
    )

    args = parser.parse_args()

    try:
        export_gnomad_bcf(
            output_path=args.output,
            af_threshold=args.af,
            genome_build=args.genome,
            chromosomes=args.chr,
            gnomad_version=args.gnomad_version,
            use_public_bucket=not args.use_requester_pays,
            gcs_project=args.gcs_project,
        )
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
