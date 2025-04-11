#!/bin/bash
# custom_annotate.sh
#
# This script is executed inside the Nextflow process "Annotate"
# It receives the following parameters:
#   1. input_bcf        - Path to the input BCF file.
#   2. input_bcf_index  - Path to the BCF index file.
#   3. bcftools_cmd     - Command to run bcftools.
#
# The script should create an annotated BCF file named "annootated.bcf"
#
# Example usage:
#   ./custom_annotate.sh input.bcf input.bcf.csi bcftools

# Exit on error
set -euo pipefail

# Read input parameters.
input_bcf="$1"
input_bcf_index="$2"
bcftools_cmd="$3"
threads="$4"
buffer_size="$5"

echo "Starting annotation with:
  BCF: ${input_bcf},
  BCFTools: ${bcftools_cmd},
  Threads: ${threads},
  Buffer Size: ${buffer_size}"

### CUSTOM ANNOTATION PROCESS --->

$bcftools_cmd view $input_bcf | \
docker run --user "$(id -u)":"$(id -g)" -i \
  -v /home/j380r/references/ensembl-vcf/113/cachedir:/home/j380r/references/ensembl-vcf/113/cachedir \
  -v /mnt/data/resources/reference/ucsc/:/mnt/data/resources/reference/ucsc/ \
  --rm ensemblorg/ensembl-vcf:release_113.0 vep \
  -a GRCh38 \
  --transcript_version \
  --total_length \
  --flag_pick \
  --exclude_predicted \
  --hgvs \
  --hgvsg \
  --spdi \
  --variant_class \
  --uniprot \
  --gene_version \
  --protein \
  --symbol \
  --canonical \
  --appris \
  --mane \
  --biotype \
  --domains \
  --refseq \
  --format vcf \
  --cache \
  --dir_cache /home/j380r/references/ensembl-vcf/113/cachedir \
  --fa /home/j380r/projects/vcfstash/tests/data/references/reference.fasta \
  --offline \
  --buffer_size $buffer_size \
  --fork $threads \
  --vcf \
  --no_stats \
  -i stdin \
  -o stdout | \
$bcftools_cmd view -Ob -o annotated.bcf --write-index

### <--- END CUSTOM ANNOTATION PROCESS

# Check if output files exist
if [[ ! -f "annotated.bcf" ]] || [[ ! -f "annotated.bcf.csi" ]]; then
    echo "Error: Output files not found"
    exit 1
fi

echo "Annotation complete. Output file: annotated.bcf"