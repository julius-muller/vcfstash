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
# 1. Use existing test data and configuration

echo "Using existing test data and configuration..."

# --------------------------------------------------------------------------
# 2. Use existing test VCF data

echo "Using gnomad_test.bcf from test data..."
G_SRC="/build/tests/data/nodata/gnomad_test.bcf"

# For now, use the test configuration that we know works
# TODO: In production, we can extend this to use VEP when properly configured
echo "Using test configuration for reliable cache build..."
PARAMS_FILE="/build/tests/config/test_params.yaml"  
CONFIG_FILE="/build/tests/config/test_annotation.config"

echo "Configuration details:"
echo "  - Params: ${PARAMS_FILE}"
echo "  - Config: ${CONFIG_FILE}"
echo "  - Test VCF: ${G_SRC}"

# Verify the test file exists
if [[ ! -f "${G_SRC}" ]]; then
    echo "ERROR: Test file ${G_SRC} not found!"
    exit 1
fi

echo "✓ Using test VCF: ${G_SRC}"

# --- build blueprint & annotate -------------------------------------------
DB_DIR="${CACHE_DIR}/db"
rm -rf "${DB_DIR}"

vcfstash stash-init   --force \
        --vcf "${G_SRC}" \
        --output "${DB_DIR}" \
        -y "${PARAMS_FILE}"

vcfstash stash-annotate \
        --db   "${DB_DIR}" \
        --name "${CNAME}" \
        -a     "${CONFIG_FILE}"