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
NETWORK=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome)        GENOME="$2"; shift 2;;
    --type)          TYPE="$2"; shift 2;;
    --tag)           TAG="$2"; shift 2;;
    --registry)      REGISTRY="$2"; shift 2;;
    --push)          PUSH=true; shift;;
    --no-cache)      NO_CACHE="--no-cache"; shift;;
    --host-network)  NETWORK="--network=host"; shift;;
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

# Prepare build context on large disk (/mnt/data/vcfstash_data)
BUILD_CONTEXT_DIR="/mnt/data/vcfstash_data/build-context"
rm -rf "${BUILD_CONTEXT_DIR}"
mkdir -p "${BUILD_CONTEXT_DIR}"

BCF_BASENAME=$(basename "${BCF_FILE}")

echo "📦 Preparing build context on /mnt/data..."
echo "  - Copying project files..."
# Copy project files (code, configs, etc - small)
# Exclude main data directory but keep tests/data
rsync -a --exclude='.git' --exclude='/data' --exclude='*.pyc' --exclude='__pycache__' \
  ./ "${BUILD_CONTEXT_DIR}/"

echo "  - Linking BCF files (hard links, no extra space)..."
# Create hard links to BCF (same filesystem, zero copy)
mkdir -p "${BUILD_CONTEXT_DIR}/docker/gnomad-data"
ln -f "${BCF_FILE}" "${BUILD_CONTEXT_DIR}/docker/gnomad-data/${BCF_BASENAME}"
ln -f "${BCF_FILE}.csi" "${BUILD_CONTEXT_DIR}/docker/gnomad-data/${BCF_BASENAME}.csi"

echo "🐳 Building Docker image..."
START_TIME=$(date +%s)

# Enable BuildKit if available (optional, legacy builder works fine)
if docker buildx version &>/dev/null; then
  echo "  Using BuildKit (modern builder)"
  export DOCKER_BUILDKIT=1
else
  echo "  Using legacy builder (buildx not installed, this is fine)"
  unset DOCKER_BUILDKIT
fi

# Build from the temporary context on /mnt/data
docker build \
  -f "${BUILD_CONTEXT_DIR}/docker/Dockerfile.cache-hail" \
  --build-arg AF="${AF}" \
  --build-arg GENOME="${GENOME}" \
  --build-arg CACHE_NAME="${CACHE_NAME}" \
  --build-arg BCF_FILE="docker/gnomad-data/${BCF_BASENAME}" \
  -t "${IMAGE_NAME}" \
  -t "${REGISTRY}/vcfstash-blueprint:latest" \
  ${NO_CACHE} \
  ${NETWORK} \
  "${BUILD_CONTEXT_DIR}"

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
echo "🧹 Cleaning up temporary build context..."
rm -rf "${BUILD_CONTEXT_DIR}"

echo "==============================================================================="
echo "Done!"
echo "==============================================================================="
echo ""
echo "To test the blueprint image:"
echo "  docker run --rm ${IMAGE_NAME} --help"
echo ""
echo "To run tests (auto-detects blueprint scenario):"
echo "  docker run --rm --entrypoint /bin/sh ${IMAGE_NAME} \\"
echo "    -c 'cd /app && export PYTHONPATH=/app/venv/lib/python3.13/site-packages && python3 -m pytest tests/ -v'"
echo ""
if [ "${PUSH}" = false ]; then
  echo "To push to registry:"
  echo "  docker push ${IMAGE_NAME}"
  echo ""
fi
