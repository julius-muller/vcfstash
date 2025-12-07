#!/usr/bin/env bash
set -euo pipefail

# --- argument parser -------------------------------------------------------
AF="0.10"
CACHE_DIR="/opt/vcfcache-cache"
THREADS="12"
CNAME="vep_gnomad"
GENOME="GRCh38"
BCF_FILE=""
PARAMS_FILE=""
ANNOTATION_CONFIG=""
ANNOTATION_NAME="vep_gnomad"
VEP_CACHE="/opt/vep/.vep"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gnomad-af)         AF="$2"; shift 2;;
    --cache-dir)         CACHE_DIR="$2"; shift 2;;
    --threads)           THREADS="$2"; shift 2;;
    --cache-name)        CNAME="$2"; shift 2;;
    --genome)            GENOME="$2"; shift 2;;
    --params)            PARAMS_FILE="$2"; shift 2;;
    --annotation-config) ANNOTATION_CONFIG="$2"; shift 2;;
    --annotation-name)   ANNOTATION_NAME="$2"; shift 2;;
    --bcf-file)          BCF_FILE="$2"; shift 2;;
    --vep-cache)         VEP_CACHE="$2"; shift 2;;
    *) echo "Unknown arg $1"; exit 1;;
  esac
done

# --------------------------------------------------------------------------
# Validate inputs

if [[ -z "${BCF_FILE}" ]] || [[ ! -f "${BCF_FILE}" ]]; then
    echo "ERROR: BCF file required and must exist: ${BCF_FILE}"
    exit 1
fi

if [[ ! -f "${BCF_FILE}.csi" ]]; then
    echo "Indexing BCF file..."
    if [ -x "/opt/bcftools/bin/bcftools" ]; then
        /opt/bcftools/bin/bcftools index "${BCF_FILE}"
    else
        bcftools index "${BCF_FILE}"
    fi
fi

if [[ -z "${PARAMS_FILE}" ]] || [[ ! -f "${PARAMS_FILE}" ]]; then
    echo "ERROR: params file required and must exist: ${PARAMS_FILE}"
    exit 1
fi

if [[ -z "${ANNOTATION_CONFIG}" ]] || [[ ! -f "${ANNOTATION_CONFIG}" ]]; then
    echo "ERROR: annotation config required and must exist: ${ANNOTATION_CONFIG}"
    exit 1
fi

# Verify VEP cache exists
if [[ ! -d "${VEP_CACHE}" ]]; then
    echo "ERROR: VEP cache directory not found: ${VEP_CACHE}"
    exit 1
fi

echo "========================================="
echo "Building Annotated VCFcache Cache"
echo "========================================="
echo "Configuration:"
echo "  - BCF Source: ${BCF_FILE}"
echo "  - AF Threshold: ${AF}"
echo "  - Genome: ${GENOME}"
echo "  - Cache Name: ${CNAME}"
echo "  - Annotation Name: ${ANNOTATION_NAME}"
echo "  - Params: ${PARAMS_FILE}"
echo "  - Annotation Config: ${ANNOTATION_CONFIG}"
echo "  - VEP Cache: ${VEP_CACHE}"
echo "  - Threads: ${THREADS}"
echo "========================================="

# --------------------------------------------------------------------------
# Step 1: Build blueprint with blueprint-init

DB_DIR="${CACHE_DIR}/db"
rm -rf "${DB_DIR}"

echo ""
echo "Step 1: Initializing VCFcache blueprint..."
vcfcache blueprint-init --force \
    --vcf "${BCF_FILE}" \
    --output "${DB_DIR}" \
    -y "${PARAMS_FILE}"

echo "✓ Blueprint created"

# --------------------------------------------------------------------------
# Step 2: Annotate blueprint with VEP (cache-build)

echo ""
echo "Step 2: Annotating blueprint with VEP..."
echo "This may take 30-60 minutes for full genome data..."

vcfcache cache-build \
    --db "${DB_DIR}" \
    --name "${ANNOTATION_NAME}" \
    -a "${ANNOTATION_CONFIG}" \
    -y "${PARAMS_FILE}"

echo "✓ Annotation complete"

# --------------------------------------------------------------------------
# Verify the annotated cache

ANNOTATED_BCF="${DB_DIR}/cache/${ANNOTATION_NAME}/vcfcache_annotated.bcf"
if [[ ! -f "${ANNOTATED_BCF}" ]]; then
    echo "ERROR: Annotated BCF not found at ${ANNOTATED_BCF}"
    exit 1
fi

echo ""
echo "========================================="
echo "✓ Annotated cache build complete!"
echo "========================================="
echo "Cache location:"
echo "  - Blueprint: ${DB_DIR}/blueprint/"
echo "  - Annotated: ${DB_DIR}/cache/${ANNOTATION_NAME}/"
echo "========================================="

# Show some stats
echo ""
echo "Cache statistics:"
# Use compiled bcftools 1.22 if available, otherwise use bcftools in PATH
if [ -x "/opt/bcftools/bin/bcftools" ]; then
    /opt/bcftools/bin/bcftools stats "${ANNOTATED_BCF}" | grep "number of records:"
else
    bcftools stats "${ANNOTATED_BCF}" | grep "number of records:"
fi
echo ""
