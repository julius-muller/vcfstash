#!/usr/bin/env python3
import argparse
import time
from datetime import datetime
from pathlib import Path
import subprocess
import warnings
import sys

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pysam
from multiprocessing import Pool

# Global variables for ease of modification
INFO_FIELDS = ['GT', 'DP', 'AF', "gnomadg_af", "gnomade_af", "gnomadg_ac", "gnomade_ac",'clinvar_clnsig', "deeprvat_score"]
TRANSCRIPT_KEYS = [
    'SYMBOL', 'Feature', 'Consequence', 'HGVS_OFFSET', 'HGVSc', 'HGVSp','IMPACT','DISTANCE','PICK','VARIANT_CLASS'
]
BASES = {"A","C","G","T"}
#Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|
# Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|
# HGNC_ID|CANONICAL|MANE|MANE_SELECT|MANE_PLUS_CLINICAL|APPRIS|ENSP|SWISSPROT|TREMBL|UNIPARC|
# UNIPROT_ISOFORM|DOMAINS|HGVS_OFFSET
# "DP", "AF", "VKX", "gnomadg_af", "gnomade_af", "clinvar_clnsig", "gnomadg_ac", "gnomade_ac"

# Weighted average function
def wavg(f1: float | None, f2: float | None, n1: int, n2: int) -> float | None:
    """Weighted average for Allele Frequencies."""
    total_weight = n1 + n2
    if total_weight == 0:
        return None  # No valid weighted average if total weight is zero
    if f1 is not None and f2 is not None:
        return (f1 * n1 + f2 * n2) / total_weight
    elif f1 is None and f2 is None:
        return None
    elif f1 is None:
        return f2
    elif f2 is None:
        return f1


def parse_vep_info(vep_data:list):
    """Parses VEP INFO field and expands transcript consequences.

    vep_data=record.info["CSQ"]
    transcript_keys=transcript_keys

    """

    def convert_vepstr(value):
        if value is None or value == "":
            return None
        try:
            if isinstance(value, str) and ("." in value or "e" in value or "E" in value):
                return float(value)
            return int(value)
        except ValueError:
            return value

    expanded_data = []
    for tn, transcript in enumerate(vep_data):
        expanded_data.append({})
        for entry in TRANSCRIPT_KEYS:
            if entry not in vep_data[tn]:
                raise ValueError(f"Did not find key={entry} in CSQ INFO tag")
            if entry == "PICK":
                expanded_data[tn][entry] = vep_data[tn][entry] == "1"
            else:
                expanded_data[tn][entry] = convert_vepstr(vep_data[tn][entry])

    return expanded_data


def process_region(args, canchr:bool = True):
    """Processes a genomic region (chromosome or chunk) and returns a DataFrame."""
    bcf_path, region = args
    # bcf_path, region = ['/home/j380r/Downloads/MG1311834374_pop13_sub1_norm_final.bcf', 'chr22']

    vcf = pysam.VariantFile(bcf_path)
    records = []
    variant_count = 0
    excluded_count = 0
    vep_csq_fields = vcf.header.info['CSQ'].description.split(' ')[-1].split('|')
    for record in vcf.fetch(region=region):
        variant_count += 1
        try:
            # Extract basic variant fields
            chrom = record.chrom
            if canchr and chrom[:3] != 'chr':
                continue
            pos = record.pos
            ref = record.ref
            alt = record.alts[0]  # Assuming single ALT for simplicity
            if not all ([x in BASES for x in alt]):
                continue
            # variant_id = record.info.get("VKX", "Unknown")

            # Extract AD and DP from FORMAT field
            if len(record.samples) == 0:
                print(f"Error: No samples found in record at {record.chrom}:{pos}")
                continue

            sample = record.samples[0]  # Get first sample
            ad = sample.get('AD', None)
            dp = sample.get('DP', None)

            # Calculate AF from AD
            if ad and len(ad) >= 2:
                ref_depth = ad[0]
                alt_depth = ad[1]
                af = alt_depth / (ref_depth + alt_depth) if (ref_depth + alt_depth) > 0 else None
            else:
                af = None

            info = {key: None for key in INFO_FIELDS}
            info |= {key: record.info.get(key, None) for key in INFO_FIELDS if key in record.info}
            try:
                info |= {
                'GT': sample.get('GT', None),
                'AD': ad[1] if ad and len(ad) > 1 else None,  # Safely access alternate allele depth
                'DP': dp,
                'AF': af
                }
            except Exception as e:
                print(f"Error processing variant at {record.chrom}:{record.pos}: {e}")
                print(sample)
                raise e

            # Parse clinvar_clnsig as a single comma-separated string
            clinvar_clnsig = info.get("clinvar_clnsig", None)
            if clinvar_clnsig and clinvar_clnsig[0] != 'null':
                info["clinvar_clnsig"] = ", ".join(clinvar_clnsig)
            else:
                info["clinvar_clnsig"] = None

                # Clean gnomadg and gnomade fields
            gnomad_fields = ["gnomadg_ac", "gnomade_ac", "gnomadg_af", "gnomade_af"]
            for field in gnomad_fields:
                value = info.get(field, None)
                if isinstance(value, (int, float)) and value < 0:
                    value = None
                info[field] = value

            # Calculate gnomad_af (weighted average)
            gnomadg_ac = float(info["gnomadg_af"]) if info.get("gnomadg_ac", None) else 0
            gnomade_ac = int(info.get("gnomade_ac", 0)) if info.get("gnomade_ac", None) else 0
            gnomadg_af = float(info["gnomadg_af"]) if info.get("gnomadg_af", None) else None
            gnomade_af = float(info["gnomade_af"]) if info.get("gnomade_af", None) else None

            if gnomadg_af is not None and gnomade_af is None:
                info["gnomad_af"] = gnomadg_af
            elif gnomade_af is not None and gnomadg_af is None:
                info["gnomad_af"] = gnomade_af
            elif gnomadg_af is not None and gnomade_af is not None:
                info["gnomad_af"] = wavg(gnomadg_af, gnomade_af, gnomadg_ac, gnomade_ac)
            else:
                info["gnomad_af"] = None

            info = {k:v for k,v in info.items() if k not in gnomad_fields}

            vep_csqs = [dict(zip(vep_csq_fields,x.split("|"))) for x in record.info["CSQ"]]

            # Expand VEP annotations
            expanded_transcripts = parse_vep_info(vep_csqs)
            # assert not (len(expanded_transcripts) > 1 and  expanded_transcripts[0]['IMPACT'] == "HIGH"), "VVV"

            for transcript in expanded_transcripts:
                row = {
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    # "VKX": variant_id,
                    **info,
                    **transcript
                }
                records.append(row)
        except Exception as e:
            excluded_count += 1
            print(f"Error processing variant at {record.chrom}:{record.pos}: {e}")
            raise

    return pd.DataFrame(records)

def store_workflow_dag(run_dir:Path, cmd:list):
    """Generate and store workflow DAG visualization"""

    try:
        del cmd[cmd.index("-with-trace")]
    except ValueError:
        pass  # Element not in list

    try:
        # Generate DAG in dot format
        dag_cmd = cmd + [ "-preview", "-with-dag", str(run_dir / "flowchart.html") ]

        subprocess.run(dag_cmd, check=True, capture_output=True)


    except subprocess.CalledProcessError as e:
        print(run_dir / "annotation.info", f"Warning: Failed to generate workflow DAG: {e}")


def validate_vcf_format(vcf_path: Path) -> tuple[bool, str]:
    """Validate VCF format fields.

    Args:
        vcf_path: Path to the VCF file

    Returns:
        Tuple of (is_valid, error_message)
    """
    try:
        vcf = pysam.VariantFile(str(vcf_path))
        required_formats = {'AD', 'DP', 'GT'}
        available_formats = set(vcf.header.formats.keys())

        missing_formats = required_formats - available_formats
        if missing_formats:
            return False, f"Missing required FORMAT fields: {', '.join(missing_formats)}"

        return True, None
    except Exception as e:
        return False, f"Error reading VCF file: {e}"

def run_vep_workflow(input_vcf: Path, output_dir: Path, workflow_dir: Path | None, cache_dir: Path | None,
                     threads: int) -> str:
    """Run VEP workflow on input VCF file.

    Args:
        input_vcf: Path to input VCF/BCF file
        output_dir: Output directory for results
        workflow_dir: Optional directory containing main.nf
        cache_dir: Optional directory containing vep_db.bcf and workflow files
        threads: Number of threads to use

    Returns:
        Path to annotated BCF file
    """
    if not input_vcf.exists():
        raise FileNotFoundError(f"Input VCF file not found: {input_vcf}")

    # Validate VCF format fields
    is_valid, error = validate_vcf_format(input_vcf)
    if not is_valid:
        raise ValueError(f"Invalid VCF file: {error}")

    sample_name = input_vcf.stem.split('.')[0]
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine workflow directory
    if workflow_dir:
        nf_path = workflow_dir / "main.nf"
    elif cache_dir:
        nf_path = cache_dir / "main.nf"
    else:
        raise ValueError("Neither workflow_dir nor cache_dir provided")

    if not nf_path.exists():
        raise FileNotFoundError(f"Workflow file not found: {nf_path}")

    # Verify database file if cache provided
    db_bcf = cache_dir / "vep_db.bcf" if cache_dir else None
    if cache_dir and not db_bcf.exists():
        raise FileNotFoundError(f"Database file not found: {db_bcf}")

    cmd = [
        "nextflow", "run", str(nf_path),
        "--input", str(input_vcf),
        "--output", str(output_dir),
        "--db_bcf", str(db_bcf) if cache_dir else "false",
        "--db_mode", "false",
        "--vep_max_forks", str(min(threads, 4)),
        "--vep_max_chr_parallel", str(min(threads, 4)),
        "-with-trace"
    ]

    print(f"[{datetime.now()}] Starting VEP workflow...")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"VEP workflow failed: {e}")

    # Store workflow DAG
    store_workflow_dag(output_dir, cmd)

    final_bcf = output_dir / f"{sample_name}_norm_final.bcf"
    if not final_bcf.exists():
        raise FileNotFoundError(f"VEP workflow output not found: {final_bcf}")

    return str(final_bcf)


def main():
    parser = argparse.ArgumentParser(description="Convert VCF to optimized Parquet format with VEP annotations")
    parser.add_argument("vcf", help="Input VCF file path")
    parser.add_argument("-c", "--cache", type=Path, help="Directory containing VEP cache and database")
    parser.add_argument("-w", "--workflow", type=Path, help="Directory containing main.nf (required if no cache)")
    parser.add_argument("-o", "--output", type=Path, help="Output directory (default: current directory)")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    args = parser.parse_args()

    start_time = time.time()

    # Input validation
    input_vcf = Path(args.vcf)
    if not input_vcf.exists():
        sys.exit(f"Error: Input VCF file not found: {input_vcf}")

    if not args.cache and not args.workflow:
        sys.exit("Error: Either --cache or --workflow must be provided")

    # Set output directory
    output_dir = args.output or Path.cwd()
    workflow_dir = args.workflow

    try:
        annotated_bcf = run_vep_workflow(input_vcf, output_dir, workflow_dir, args.cache, args.threads)
        vep_time = time.time()
        print(f"[{datetime.now()}] VEP workflow completed in {vep_time - start_time:.1f}s")
    except Exception as e:
        sys.exit(f"Error running VEP workflow: {e}")

    print(f"[{datetime.now()}] Processing BCF: {annotated_bcf}")
    vcf = pysam.VariantFile(annotated_bcf)
    regions = list(vcf.header.contigs.keys())
    args_list = [(annotated_bcf, region) for region in regions]

    print(f"[{datetime.now()}] Processing {len(regions)} regions with {args.threads} threads...")
    with Pool(args.threads) as pool:
        dataframes = pool.map(process_region, args_list)

    process_time = time.time()
    print(f"[{datetime.now()}] Region processing completed in {process_time - vep_time:.1f}s")

    # Filter out empty dataframes
    dataframes = [df for df in dataframes if not df.empty and not df.isna().all(axis=None)]
    if not dataframes:
        sys.exit("Error: No valid variants found in input file")

    # Combine and write output
    combined_df = pd.concat(dataframes, ignore_index=True)
    output_file = output_dir / f"{input_vcf.stem}.parquet"

    print(f"[{datetime.now()}] Writing optimized Parquet file: {output_file}")
    table = pa.Table.from_pandas(combined_df)
    pq.write_table(
        table,
        output_file,
        compression="snappy",
        use_dictionary=True,
        row_group_size=100000,
        data_page_size=65536,
        write_statistics=True,
    )

    end_time = time.time()
    print(f"[{datetime.now()}] Parquet conversion completed in {end_time - process_time:.1f}s")
    print(f"[{datetime.now()}] Total execution time: {end_time - start_time:.1f}s")


if __name__ == "__main__":
    main()
