#!/usr/bin/env bash
set -euo pipefail

#===============================================================================
# Subset Master BCF to Different AF Thresholds
#===============================================================================
# This script takes the master BCF (AF >= 0.001) and creates subsets for
# different allele frequency thresholds. This is much faster than regenerating
# with Hail.
#
# Usage:
#   ./02-subset-af-thresholds.sh MASTER_BCF [OPTIONS]
#
# Arguments:
#   MASTER_BCF              Path to master BCF file
#
# Options:
#   --af-thresholds "X Y Z" Space-separated AF thresholds (default: "0.01 0.05 0.10")
#   --output-dir DIR        Output directory (default: ./data)
#   --threads N             Number of threads (default: 4)
#===============================================================================

# Check arguments
if [ $# -lt 1 ]; then
  grep "^#" "$0" | grep -v "#!/" | sed 's/^# \?//'
  exit 1
fi

MASTER_BCF="$1"
shift

# Default configuration
AF_THRESHOLDS="0.01 0.05 0.10"
OUTPUT_DIR="/mnt/data/vcfstash_data/gnomad"
THREADS=4

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --af-thresholds) AF_THRESHOLDS="$2"; shift 2;;
    --output-dir)    OUTPUT_DIR="$2"; shift 2;;
    --threads)       THREADS="$2"; shift 2;;
    --help)
      grep "^#" "$0" | grep -v "#!/" | sed 's/^# \?//'
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Validate master BCF
if [ ! -f "${MASTER_BCF}" ]; then
  echo "❌ ERROR: Master BCF not found: ${MASTER_BCF}"
  exit 1
fi

if [ ! -f "${MASTER_BCF}.csi" ]; then
  echo "⚠️  Master BCF index not found, creating..."
  bcftools index "${MASTER_BCF}"
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Extract base filename (without path and extension)
BASENAME=$(basename "${MASTER_BCF}" .bcf)
# Remove "_master_af0001" suffix to get base name
BASEPREFIX=$(echo "${BASENAME}" | sed 's/_master_af[0-9]*$//')

echo "==============================================================================="
echo "Subsetting Master BCF to AF Thresholds"
echo "==============================================================================="
echo "Master BCF:      ${MASTER_BCF}"
echo "Output dir:      ${OUTPUT_DIR}"
echo "AF thresholds:   ${AF_THRESHOLDS}"
echo "Threads:         ${THREADS}"
echo "==============================================================================="
echo ""

# Get master BCF stats
echo "Master BCF statistics:"
bcftools stats "${MASTER_BCF}" | grep "number of records:"
echo ""

# Process each AF threshold
for AF in ${AF_THRESHOLDS}; do
  AF_CLEAN=$(echo "${AF}" | tr -d '.')
  OUTPUT_FILE="${OUTPUT_DIR}/${BASEPREFIX}_af${AF_CLEAN}.bcf"

  echo "───────────────────────────────────────────────────────────────────────────"
  echo "Processing AF >= ${AF} ..."
  echo "───────────────────────────────────────────────────────────────────────────"

  # Check if output exists
  if [ -f "${OUTPUT_FILE}" ]; then
    echo "⚠️  Output exists: ${OUTPUT_FILE}"
    read -p "Overwrite? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
      echo "Skipping."
      echo ""
      continue
    fi
    rm -f "${OUTPUT_FILE}" "${OUTPUT_FILE}.csi"
  fi

  # Subset with bcftools
  echo "Filtering variants with AF >= ${AF}..."
  START_TIME=$(date +%s)

  bcftools view \
    --threads "${THREADS}" \
    -i "INFO/AF>=${AF}" \
    -Ob \
    -o "${OUTPUT_FILE}" \
    "${MASTER_BCF}"

  # Index the output
  echo "Indexing..."
  bcftools index "${OUTPUT_FILE}"

  END_TIME=$(date +%s)
  DURATION=$((END_TIME - START_TIME))

  # Show stats
  echo "✅ Created: ${OUTPUT_FILE}"
  bcftools stats "${OUTPUT_FILE}" | grep "number of records:"
  echo "   Duration: ${DURATION} seconds"
  echo ""
done

echo "==============================================================================="
echo "✅ All subsets created successfully!"
echo "==============================================================================="
echo ""
echo "Created files:"
for AF in ${AF_THRESHOLDS}; do
  AF_CLEAN=$(echo "${AF}" | tr -d '.')
  OUTPUT_FILE="${OUTPUT_DIR}/${BASEPREFIX}_af${AF_CLEAN}.bcf"
  if [ -f "${OUTPUT_FILE}" ]; then
    SIZE=$(du -h "${OUTPUT_FILE}" | cut -f1)
    echo "  - ${OUTPUT_FILE} (${SIZE})"
  fi
done

echo ""
echo "Next steps:"
echo "  1. Run 03-build-blueprint.sh to build blueprint Docker images"
echo "  2. Run 04-build-annotated.sh to build annotated Docker images"
echo ""
