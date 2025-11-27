#!/usr/bin/env bash
set -euo pipefail

#===============================================================================
# Generate Master BCF from gnomAD with minimal AF filtering
#===============================================================================
# This script generates a master BCF file with minimal allele frequency
# filtering (AF >= 0.001 = 0.1%). This master file can then be quickly
# subsetted to different AF thresholds without re-running Hail.
#
# Usage:
#   ./01-generate-master-bcf.sh [OPTIONS]
#
# Options:
#   --genome GENOME         Genome build (default: GRCh38)
#   --gnomad-version VER    gnomAD version (default: 4.1)
#   --type TYPE             Dataset type: joint/genomes/exomes (default: joint)
#   --chromosomes CHR       Chromosomes to include: all/chrY/chr22 (default: all)
#   --output-dir DIR        Output directory (default: ./data)
#   --master-af AF          Master AF threshold (default: 0.001)
#===============================================================================

# Default configuration
GENOME="GRCh38"
GNOMAD_VERSION="4.1"
TYPE="joint"
CHROMOSOMES="all"
OUTPUT_DIR="./data"
MASTER_AF="0.001"

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome)         GENOME="$2"; shift 2;;
    --gnomad-version) GNOMAD_VERSION="$2"; shift 2;;
    --type)           TYPE="$2"; shift 2;;
    --chromosomes)    CHROMOSOMES="$2"; shift 2;;
    --output-dir)     OUTPUT_DIR="$2"; shift 2;;
    --master-af)      MASTER_AF="$2"; shift 2;;
    --help)
      grep "^#" "$0" | grep -v "#!/" | sed 's/^# \?//'
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use --help for usage information"
      exit 1
      ;;
  esac
done

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Generate output filename
GENOME_LOWER=$(echo "${GENOME}" | tr '[:upper:]' '[:lower:]')
AF_CLEAN=$(echo "${MASTER_AF}" | tr -d '.')
if [ "${CHROMOSOMES}" = "all" ]; then
  OUTPUT_FILE="${OUTPUT_DIR}/gnomad_v${GNOMAD_VERSION}_${GENOME}_${TYPE}_master_af${AF_CLEAN}.bcf"
else
  CHR_CLEAN=$(echo "${CHROMOSOMES}" | tr ' ' '_')
  OUTPUT_FILE="${OUTPUT_DIR}/gnomad_v${GNOMAD_VERSION}_${GENOME}_${TYPE}_${CHR_CLEAN}_master_af${AF_CLEAN}.bcf"
fi

# Build chromosome arguments
CHR_ARGS=""
if [ "${CHROMOSOMES}" != "all" ]; then
  CHR_ARGS="--chr ${CHROMOSOMES}"
fi

echo "==============================================================================="
echo "Generating Master BCF from gnomAD"
echo "==============================================================================="
echo "Configuration:"
echo "  Genome:           ${GENOME}"
echo "  gnomAD version:   ${GNOMAD_VERSION}"
echo "  Dataset type:     ${TYPE}"
echo "  Chromosomes:      ${CHROMOSOMES}"
echo "  Master AF:        ${MASTER_AF} (0.1%)"
echo "  Output:           ${OUTPUT_FILE}"
echo "==============================================================================="
echo ""

# Check if output already exists
if [ -f "${OUTPUT_FILE}" ]; then
  echo "⚠️  WARNING: Output file already exists: ${OUTPUT_FILE}"
  read -p "Overwrite? (y/N): " -n 1 -r
  echo
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 1
  fi
  rm -f "${OUTPUT_FILE}" "${OUTPUT_FILE}.csi"
fi

# Check if Hail is installed
if ! python3 -c "import hail" 2>/dev/null; then
  echo "❌ ERROR: Hail is not installed"
  echo "Install with: pip install hail"
  exit 1
fi

# Check if bcftools is installed
if ! command -v bcftools &>/dev/null; then
  echo "❌ ERROR: bcftools is not installed"
  echo "Install with: sudo apt-get install bcftools"
  exit 1
fi

# Run the Hail export script
echo "🚀 Starting Hail export (this may take 1-3 hours for full genome)..."
echo ""

START_TIME=$(date +%s)

python3 scripts/export_gnomad_hail.py \
  --output "${OUTPUT_FILE}" \
  --af "${MASTER_AF}" \
  --genome "${GENOME}" \
  --gnomad-version "${GNOMAD_VERSION}" \
  ${CHR_ARGS}

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
DURATION_MIN=$((DURATION / 60))

echo ""
echo "==============================================================================="
echo "✅ Master BCF generated successfully!"
echo "==============================================================================="
echo "  File:      ${OUTPUT_FILE}"
echo "  Duration:  ${DURATION_MIN} minutes"
echo ""

# Show statistics
echo "Statistics:"
bcftools stats "${OUTPUT_FILE}" | grep -E "^SN"

echo ""
echo "Next steps:"
echo "  1. Run 02-subset-af-thresholds.sh to create AF-filtered versions"
echo "  2. Run 03-build-blueprint.sh to build Docker images"
echo ""
