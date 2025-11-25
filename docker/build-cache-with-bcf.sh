#!/usr/bin/env bash
set -euo pipefail

# --- argument parser -------------------------------------------------------
AF="0.10"
CACHE_DIR="/home/micromamba/cache"
THREADS="8"
TOOL_VER="115.2"
CNAME="vep_gnomad"
GENOME="GRCh38"
BCF_FILE=""  # Pre-generated BCF file

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
    --bcf-file)     BCF_FILE="$2"; shift 2;;
    *) echo "Unknown arg $1"; exit 1;;
  esac
done

# --------------------------------------------------------------------------
# Setup configuration files

echo "Setting up cache with pre-generated BCF file"

# Use proper VEP configuration (not test config)
PARAMS_FILE="${PARAMS_FILE:-/tmp/recipe/params.yaml}"
CONFIG_FILE="${CONFIG_FILE:-/tmp/recipe/annotation.config}"

# Create working directory
WORK_DIR="/tmp/gnomad_work"
mkdir -p "${WORK_DIR}"

# --------------------------------------------------------------------------
# Use pre-generated BCF file

if [[ -n "${BCF_FILE}" ]] && [[ -f "${BCF_FILE}" ]]; then
    echo "Using pre-generated BCF file: ${BCF_FILE}"
    G_SRC="${BCF_FILE}"

    # Verify index exists
    if [[ ! -f "${G_SRC}.csi" ]]; then
        echo "Indexing BCF file..."
        bcftools index "${G_SRC}"
    fi

    echo "✓ Using BCF file: ${G_SRC}"
else
    echo "ERROR: No BCF file provided or file not found: ${BCF_FILE}"
    echo "Please provide a BCF file using --bcf-file option"
    exit 1
fi

echo "Configuration details:"
echo "  - Params: ${PARAMS_FILE}"
echo "  - Config: ${CONFIG_FILE}"
echo "  - BCF Source: ${G_SRC}"
echo "  - AF Threshold: ${AF}"

# Verify the source file exists
if [[ ! -f "${G_SRC}" ]]; then
    echo "ERROR: Source file ${G_SRC} not found!"
    exit 1
fi

# --- build blueprint & annotate -------------------------------------------
DB_DIR="${CACHE_DIR}/db"
rm -rf "${DB_DIR}"

echo "Initializing VCFstash cache..."
vcfstash stash-init --force \
        --vcf "${G_SRC}" \
        --output "${DB_DIR}" \
        -y "${PARAMS_FILE}"

echo "Annotating cache..."
vcfstash stash-annotate \
        --db   "${DB_DIR}" \
        --name "${CNAME}" \
        -a     "${CONFIG_FILE}"

echo "✓ Cache build complete!"
echo "  - Cache directory: ${DB_DIR}"
echo "  - Annotation name: ${CNAME}"
