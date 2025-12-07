#!/usr/bin/env bash
set -euo pipefail

# --- argument parser -------------------------------------------------------
AF="0.10"
CACHE_DIR="/home/micromamba/cache"
THREADS="8"
TOOL_VER="115.2"
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
# 1. Download and filter real gnomAD data

echo "Setting up real gnomAD data for AF threshold: ${AF}"

# Use proper VEP configuration (not test config)
PARAMS_FILE="${PARAMS_FILE:-/tmp/recipe/params.yaml}"
CONFIG_FILE="${CONFIG_FILE:-/tmp/recipe/annotation.config}"

# Create working directory for gnomAD data
WORK_DIR="/tmp/gnomad_work"
mkdir -p "${WORK_DIR}"

# --------------------------------------------------------------------------  
# 2. Download gnomAD data filtered by allele frequency

if [[ -n "${GNOMAD_URL}" ]]; then
    echo "Using provided gnomAD URL: ${GNOMAD_URL}"
    G_SRC="${WORK_DIR}/gnomad_filtered.vcf.gz"
    
    # Download and filter gnomAD data
    echo "Downloading and filtering gnomAD data (AF >= ${AF})..."
    curl -L "${GNOMAD_URL}" | \
    bcftools view -i "INFO/AF >= ${AF}" -Oz -o "${G_SRC}"
    
    # Index the filtered file
    tabix -p vcf "${G_SRC}"
    
    echo "✓ Downloaded and filtered gnomAD data: ${G_SRC}"
else
    echo "WARNING: No GNOMAD_URL provided, falling back to test data"
    echo "This will result in a cache with minimal coverage!"
    G_SRC="/build/tests/data/nodata/gnomad_test.bcf"
    PARAMS_FILE="/build/tests/config/test_params.yaml"  
    CONFIG_FILE="/build/tests/config/test_annotation.config"
fi

echo "Configuration details:"
echo "  - Params: ${PARAMS_FILE}"
echo "  - Config: ${CONFIG_FILE}"  
echo "  - VCF Source: ${G_SRC}"
echo "  - AF Threshold: ${AF}"

# Verify the source file exists
if [[ ! -f "${G_SRC}" ]]; then
    echo "ERROR: Source file ${G_SRC} not found!"
    exit 1
fi

# --- build blueprint & annotate -------------------------------------------
DB_DIR="${CACHE_DIR}/db"
rm -rf "${DB_DIR}"

vcfcache blueprint-init   --force \
        --vcf "${G_SRC}" \
        --output "${DB_DIR}" \
        -y "${PARAMS_FILE}"

vcfcache cache-build \
        --db   "${DB_DIR}" \
        --name "${CNAME}" \
        -a     "${CONFIG_FILE}"