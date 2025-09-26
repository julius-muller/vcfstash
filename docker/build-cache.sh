#!/usr/bin/env bash
set -euo pipefail

# --- argument parser -------------------------------------------------------
AF="0.10"
CACHE_DIR="/home/micromamba/cache"
THREADS="8"
TOOL_VER="115.1"
CNAME="vep_gnomad"
GENOME="GRCh38"
GNOMAD_URL="${GNOMAD_URL:-}"   # will come from Docker ARG/ENV, may be overridden here

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gnomad-af)    AF="$2"; shift 2;;
    --cache-dir)    CACHE_DIR="$2"; shift 2;;
    --threads)      THREADS="$2"; shift 2;;
    --tool-version) TOOL_VER="$2"; shift 2;;
    --cache-name)   CNAME="$2"; shift 2;;
    --genome)       GENOME="$2"; shift 2;;
    --params)       PARAMS_FILE="$2"; shift 2;;
    --config)       CONFIG_FILE="$2"; shift 2;;
    --url|-u)       GNOMAD_URL="$2"; shift 2;;
    *) echo "Unknown arg $1"; exit 1;;
  esac
done

# --------------------------------------------------------------------------
# 1. pick a default URL if none was passed (tiny file good for CI)
if [[ -z "${GNOMAD_URL}" ]]; then
  GNOMAD_URL="https://raw.githubusercontent.com/vcftools/vcftools/master/examples/sequence.vcf.gz"
fi

mkdir -p "${CACHE_DIR}"

# --------------------------------------------------------------------------
# 2. Download & index
G_SRC="/tmp/gnomad.${GENOME}.vcf.gz"
echo "Downloading gnomAD slice from: ${GNOMAD_URL}"
curl -L "${GNOMAD_URL}" -o "${G_SRC}"
tabix -p vcf "${G_SRC}"

echo "Filtering for AF >= ${AF}"
bcftools view -i "AF>=${AF}" "${G_SRC}" \
        -Ob -o /tmp/gnomad_af.bcf --threads "${THREADS}"
bcftools index /tmp/gnomad_af.bcf

# --------------------------------------------------------------------------
# 3. Build blueprint & annotate
vcfstash stash-init \
        --vcf /tmp/gnomad_af.bcf \
        --output "${CACHE_DIR}" \
        -y "${PARAMS_FILE}"

vcfstash stash-annotate \
        --db    "${CACHE_DIR}" \
        --name  "${CNAME}" \
        -a      "${CONFIG_FILE}"

echo "Cache created in ${CACHE_DIR}"