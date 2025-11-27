#!/usr/bin/env bash
set -euo pipefail

#===============================================================================
# Build Blueprint Docker Image Locally
#===============================================================================
# This script builds the vcfstash-blueprint (lean) Docker image from a
# pre-generated BCF file.
#
# Usage:
#   ./03-build-blueprint.sh BCF_FILE [OPTIONS]
#
# Arguments:
#   BCF_FILE                Path to gnomAD BCF file
#
# Options:
#   --genome GENOME         Genome build (default: GRCh38)
#   --type TYPE             Dataset type (default: joint)
#   --tag TAG               Docker tag (default: auto-generated)
#   --registry REGISTRY     Docker registry (default: ghcr.io/julius-muller)
#   --push                  Push to registry after build
#   --no-cache              Build without using Docker cache
#===============================================================================

# Check arguments
if [ $# -lt 1 ]; then
  grep "^#" "$0" | grep -v "#!/" | sed 's/^# \?//'
  exit 1
fi

BCF_FILE="$1"
shift

# Default configuration
GENOME="GRCh38"
TYPE="joint"
TAG=""
REGISTRY="ghcr.io/julius-muller"
PUSH=false
NO_CACHE=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome)    GENOME="$2"; shift 2;;
    --type)      TYPE="$2"; shift 2;;
    --tag)       TAG="$2"; shift 2;;
    --registry)  REGISTRY="$2"; shift 2;;
    --push)      PUSH=true; shift;;
    --no-cache)  NO_CACHE="--no-cache"; shift;;
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

# Validate BCF file
if [ ! -f "${BCF_FILE}" ]; then
  echo "❌ ERROR: BCF file not found: ${BCF_FILE}"
  exit 1
fi

if [ ! -f "${BCF_FILE}.csi" ]; then
  echo "⚠️  BCF index not found, creating..."
  bcftools index "${BCF_FILE}"
fi

# Extract AF from filename (e.g., gnomad_v4.1_GRCh38_joint_af010.bcf -> 010)
BASENAME=$(basename "${BCF_FILE}" .bcf)
if [[ $BASENAME =~ _af([0-9]+) ]]; then
  AF_CLEAN="${BASH_REMATCH[1]}"
  # Convert to decimal (010 -> 0.10)
  AF="0.${AF_CLEAN}"
else
  echo "❌ ERROR: Cannot extract AF from filename: ${BASENAME}"
  echo "Expected format: *_af010.bcf (for AF 0.10)"
  exit 1
fi

# Auto-generate tag if not provided
if [ -z "${TAG}" ]; then
  GENOME_LOWER=$(echo "${GENOME}" | tr '[:upper:]' '[:lower:]')
  TAG="gnomad-${GENOME_LOWER}-${TYPE}-af${AF_CLEAN}"
fi

# Generate cache name
GENOME_LOWER=$(echo "${GENOME}" | tr '[:upper:]' '[:lower:]')
CACHE_NAME="gnomad_${GENOME_LOWER}_${TYPE}_af${AF_CLEAN}"

# Full image name
IMAGE_NAME="${REGISTRY}/vcfstash-blueprint:${TAG}"

echo "==============================================================================="
echo "Building Blueprint Docker Image"
echo "==============================================================================="
echo "BCF file:        ${BCF_FILE}"
echo "Genome:          ${GENOME}"
echo "Type:            ${TYPE}"
echo "AF threshold:    ${AF}"
echo "Cache name:      ${CACHE_NAME}"
echo "Image tag:       ${TAG}"
echo "Full image:      ${IMAGE_NAME}"
echo "Push to GHCR:    ${PUSH}"
echo "==============================================================================="
echo ""

# Prepare build context
BUILD_DIR="./docker/build-context"
mkdir -p "${BUILD_DIR}/gnomad-data"

echo "📦 Copying BCF to build context..."
cp "${BCF_FILE}" "${BUILD_DIR}/gnomad-data/"
cp "${BCF_FILE}.csi" "${BUILD_DIR}/gnomad-data/"

BCF_BASENAME=$(basename "${BCF_FILE}")

echo "🐳 Building Docker image..."
START_TIME=$(date +%s)

docker build \
  -f docker/Dockerfile.cache-hail \
  --build-arg AF="${AF}" \
  --build-arg GENOME="${GENOME}" \
  --build-arg CACHE_NAME="${CACHE_NAME}" \
  --build-arg BCF_FILE="gnomad-data/${BCF_BASENAME}" \
  -t "${IMAGE_NAME}" \
  -t "${REGISTRY}/vcfstash-blueprint:latest" \
  ${NO_CACHE} \
  --progress=plain \
  -f docker/Dockerfile.cache-hail \
  .

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
DURATION_MIN=$((DURATION / 60))

echo ""
echo "==============================================================================="
echo "✅ Blueprint image built successfully!"
echo "==============================================================================="
echo "Image:     ${IMAGE_NAME}"
echo "Duration:  ${DURATION_MIN} minutes"
echo ""

# Get image size
IMAGE_SIZE=$(docker images "${IMAGE_NAME}" --format "{{.Size}}")
echo "Size:      ${IMAGE_SIZE}"
echo ""

# Test the image
echo "🧪 Testing image..."
docker run --rm "${IMAGE_NAME}" -v

echo ""
echo "✅ Image verification passed!"
echo ""

# Push if requested
if [ "${PUSH}" = true ]; then
  echo "📤 Pushing to registry..."
  docker push "${IMAGE_NAME}"
  docker push "${REGISTRY}/vcfstash-blueprint:latest"
  echo "✅ Pushed to ${REGISTRY}"
  echo ""
fi

# Cleanup
echo "🧹 Cleaning up build context..."
rm -rf "${BUILD_DIR}"

echo "==============================================================================="
echo "Done!"
echo "==============================================================================="
echo ""
echo "To test the blueprint image:"
echo "  docker run --rm ${IMAGE_NAME} --help"
echo ""
echo "To run tests:"
echo "  docker run --rm --entrypoint /bin/sh ${IMAGE_NAME} \\"
echo "    -c 'cd /app && export PYTHONPATH=/app/venv/lib/python3.13/site-packages && python3 -m pytest -m blueprint tests/ -v'"
echo ""
if [ "${PUSH}" = false ]; then
  echo "To push to registry:"
  echo "  docker push ${IMAGE_NAME}"
  echo ""
fi
