#!/usr/bin/env bash
set -euo pipefail

#===============================================================================
# Build Annotated Base Image (WITHOUT annotation)
#===============================================================================
# This script builds ONLY the Docker base image with VEP installed.
# It does NOT run the annotation step - that's done by 04b-annotate-and-commit.sh
#
# Usage:
#   ./04a-build-base-image.sh BCF_FILE [OPTIONS]
#
# Arguments:
#   BCF_FILE                Path to gnomAD BCF file
#
# Options:
#   --genome GENOME         Genome build (default: GRCh38)
#   --type TYPE             Dataset type (default: joint)
#   --vep-version VER       VEP version (default: 115.2)
#   --annotation-name NAME  Annotation stash name (default: vep_gnomad)
#   --tag TAG               Docker tag (default: auto-generated)
#   --registry REGISTRY     Docker registry (default: ghcr.io/julius-muller)
#   --no-cache              Build without using Docker cache
#   --host-network          Use host network for Docker build (fixes PyPI timeouts)
#   --yes|-y                Skip confirmation prompt
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
VEP_VERSION="115.2"
ANNOTATION_NAME="vep_gnomad"
TAG=""
REGISTRY="ghcr.io/julius-muller"
NO_CACHE=""
HOST_NETWORK=""
YES=false

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome)           GENOME="$2"; shift 2;;
    --type)             TYPE="$2"; shift 2;;
    --vep-version)      VEP_VERSION="$2"; shift 2;;
    --annotation-name)  ANNOTATION_NAME="$2"; shift 2;;
    --tag)              TAG="$2"; shift 2;;
    --registry)         REGISTRY="$2"; shift 2;;
    --no-cache)         NO_CACHE="--no-cache"; shift;;
    --host-network)     HOST_NETWORK="--network host"; shift;;
    --yes|-y)           YES=true; shift;;
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

# Extract AF from filename
BASENAME=$(basename "${BCF_FILE}" .bcf)
if [[ $BASENAME =~ _af([0-9]+) ]]; then
  AF_CLEAN="${BASH_REMATCH[1]}"
  # Strip leading zeros to avoid octal interpretation (010 -> 10)
  AF_NUM=$((10#${AF_CLEAN}))
  # Convert to decimal: af001 -> 0.01, af010 -> 0.10, af005 -> 0.05
  # Formula: AF = (number after af) / 100
  AF=$(awk "BEGIN {printf \"%.2f\", ${AF_NUM}/100}")
else
  echo "❌ ERROR: Cannot extract AF from filename: ${BASENAME}"
  echo "Expected format: *_af010.bcf (for AF 0.10) or *_af001.bcf (for AF 0.01)"
  exit 1
fi

# Auto-generate tag if not provided (includes AF threshold for unique tags)
if [ -z "${TAG}" ]; then
  GENOME_LOWER=$(echo "${GENOME}" | tr '[:upper:]' '[:lower:]')
  TYPE_LOWER=$(echo "${TYPE}" | tr '[:upper:]' '[:lower:]')
  VEP_MAJOR=$(echo "${VEP_VERSION}" | cut -d. -f1)
  TAG="gnomad-v41-${GENOME_LOWER}-${TYPE_LOWER}-af${AF_CLEAN}-vep${VEP_MAJOR}-py"
fi

# Generate cache name
GENOME_LOWER=$(echo "${GENOME}" | tr '[:upper:]' '[:lower:]')
CACHE_NAME="gnomad_${GENOME_LOWER}_${TYPE}_af${AF_CLEAN}"

# Full image name (base image, not final)
IMAGE_NAME="${REGISTRY}/vcfstash-annotated:${TAG}-base"

echo "==============================================================================="
echo "Building Annotated Base Image (Pure Python - Docker only, NO annotation)"
echo "==============================================================================="
echo "BCF file:           ${BCF_FILE}"
echo "Genome:             ${GENOME}"
echo "Type:               ${TYPE}"
echo "AF threshold:       ${AF}"
echo "VEP version:        ${VEP_VERSION}"
echo "Annotation name:    ${ANNOTATION_NAME}"
echo "Cache name:         ${CACHE_NAME}"
echo "Image tag:          ${TAG}-base"
echo "Full image:         ${IMAGE_NAME}"
echo "==============================================================================="
echo ""
if [ "${YES}" = false ]; then
  read -p "Continue? (y/N): " -n 1 -r
  echo
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 1
  fi
else
  echo "Auto-confirming (--yes flag provided)"
fi

# Prepare build context on large disk (/mnt/data/vcfstash_data)
BUILD_CONTEXT_DIR="/mnt/data/vcfstash_data/build-context"
rm -rf "${BUILD_CONTEXT_DIR}"
mkdir -p "${BUILD_CONTEXT_DIR}"

BCF_BASENAME=$(basename "${BCF_FILE}")

echo "📦 Preparing build context on /mnt/data..."
echo "  - Copying project files..."
# Copy project files (code, configs, etc - small)
rsync -a --exclude='.git' --exclude='/data' --exclude='*.pyc' --exclude='__pycache__' \
  ./ "${BUILD_CONTEXT_DIR}/"

echo "  - Linking BCF files (hard links, no extra space)..."
# Create hard links to BCF (same filesystem, zero copy)
mkdir -p "${BUILD_CONTEXT_DIR}/docker/gnomad-data"
ln -f "${BCF_FILE}" "${BUILD_CONTEXT_DIR}/docker/gnomad-data/${BCF_BASENAME}"
ln -f "${BCF_FILE}.csi" "${BUILD_CONTEXT_DIR}/docker/gnomad-data/${BCF_BASENAME}.csi"

echo "🐳 Building Docker base image..."
START_TIME=$(date +%s)

# Enable BuildKit if available
if docker buildx version &>/dev/null; then
  echo "  Using BuildKit (modern builder)"
  export DOCKER_BUILDKIT=1
else
  echo "  Using legacy builder (buildx not installed, this is fine)"
  unset DOCKER_BUILDKIT
fi

# Build the base image
docker build \
  -f "${BUILD_CONTEXT_DIR}/docker/Dockerfile.annotated" \
  --build-arg AF="${AF}" \
  --build-arg GENOME="${GENOME}" \
  --build-arg CACHE_NAME="${CACHE_NAME}" \
  --build-arg BCF_FILE="docker/gnomad-data/${BCF_BASENAME}" \
  --build-arg ANNOTATION_NAME="${ANNOTATION_NAME}" \
  -t "${IMAGE_NAME}" \
  ${NO_CACHE} \
  ${HOST_NETWORK} \
  "${BUILD_CONTEXT_DIR}"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
DURATION_MIN=$((DURATION / 60))

echo ""
echo "==============================================================================="
echo "✅ Base image built successfully!"
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

# Verify VEP is available
echo "🧪 Verifying VEP..."
docker run --rm --entrypoint vep "${IMAGE_NAME}" --help | head -5

echo ""
echo "✅ Base image ready!"
echo ""

# Cleanup
echo "🧹 Cleaning up temporary build context..."
rm -rf "${BUILD_CONTEXT_DIR}"

echo "==============================================================================="
echo "Done!"
echo "==============================================================================="
echo ""
echo "Next step: Run annotation with 04b-annotate-and-commit.sh"
echo ""
echo "  ./scripts/local-build/04b-annotate-and-commit.sh \\"
echo "    --base-image ${IMAGE_NAME} \\"
echo "    --vep-cache-dir /path/to/vep/cache"
echo ""
