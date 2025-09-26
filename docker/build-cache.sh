#!/usr/bin/env bash
set -euo pipefail

# Args (with defaults)
AF="0.10"
CACHE_DIR="/opt/vcfstash-cache"
THREADS="8"
TOOL_VER="115.1"
CNAME="vep_gnomad"
GENOME="GRCh38"
PARAMS_FILE=""
CONFIG_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gnomad-af) AF="$2"; shift 2;;
    --cache-dir) CACHE_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --tool-version) TOOL_VER="$2"; shift 2;;
    --cache-name) CNAME="$2"; shift 2;;
    --genome) GENOME="$2"; shift 2;;
    --params) PARAMS_FILE="$2"; shift 2;;
    --config) CONFIG_FILE="$2"; shift 2;;
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

if [[ -z "${PARAMS_FILE}" || -z "${CONFIG_FILE}" ]]; then
  echo "ERROR: --params and --config are required"
  exit 1
fi

mkdir -p "${CACHE_DIR}"

# --- Download gnomAD genomes for specified build and filter by AF -----------
# Adjust URL if you prefer exomes or a different mirror.
G_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/${GENOME}/gnomad.genomes.vcf.bgz"
G_SRC="/tmp/gnomad.${GENOME}.vcf.bgz"
echo "Downloading gnomAD from: ${G_URL}"
curl -L "${G_URL}" -o "${G_SRC}"
tabix -p vcf "${G_SRC}"

echo "Filtering gnomAD for AF >= ${AF}"
bcftools view -i "AF>=${AF}" "${G_SRC}" -Ob -o /tmp/gnomad_af.bcf --threads "${THREADS}"
bcftools index /tmp/gnomad_af.bcf

# --- Build blueprint and annotate with VEP 115.1 ---------------------------
echo "Running vcfstash stash-init"
vcfstash stash-init \
  --vcf /tmp/gnomad_af.bcf \
  --output "${CACHE_DIR}" \
  -y "${PARAMS_FILE}"

echo "Running vcfstash stash-annotate"
vcfstash stash-annotate \
  --db "${CACHE_DIR}" \
  --name "${CNAME}" \
  -a "${CONFIG_FILE}"

echo "Cache ready at: ${CACHE_DIR}"