#!/usr/bin/env bash
set -euo pipefail

#===============================================================================
# Build Annotated Docker Image Locally
#===============================================================================
# This script builds the vcfstash-annotated (fat) Docker image with VEP
# pipeline and pre-annotated cache from a pre-generated BCF file.
#
# Usage:
#   ./04-build-annotated.sh BCF_FILE [OPTIONS]
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
#   --push                  Push to registry after build
#   --no-cache              Build without using Docker cache
#   --host-network          Use host network for Docker build (fixes PyPI timeouts)
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
PUSH=false
NO_CACHE=""
HOST_NETWORK=""
VEP_CACHE_DIR=""
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
    --push)             PUSH=true; shift;;
    --no-cache)         NO_CACHE="--no-cache"; shift;;
    --host-network)     HOST_NETWORK="--network host"; shift;;
    --vep-cache-dir)    VEP_CACHE_DIR="$2"; shift 2;;
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
  TAG="gnomad-v41-${GENOME_LOWER}-${TYPE_LOWER}-af${AF_CLEAN}-vep${VEP_MAJOR}"
fi

# Generate cache name
GENOME_LOWER=$(echo "${GENOME}" | tr '[:upper:]' '[:lower:]')
CACHE_NAME="gnomad_${GENOME_LOWER}_${TYPE}_af${AF_CLEAN}"

# Full image name
IMAGE_NAME="${REGISTRY}/vcfstash-annotated:${TAG}"

echo "==============================================================================="
echo "Building Annotated Docker Image (with VEP)"
echo "==============================================================================="
echo "BCF file:           ${BCF_FILE}"
echo "Genome:             ${GENOME}"
echo "Type:               ${TYPE}"
echo "AF threshold:       ${AF}"
echo "VEP version:        ${VEP_VERSION}"
echo "Annotation name:    ${ANNOTATION_NAME}"
echo "Cache name:         ${CACHE_NAME}"
echo "Image tag:          ${TAG}"
echo "Full image:         ${IMAGE_NAME}"
echo "Push to GHCR:       ${PUSH}"
if [ -n "${VEP_CACHE_DIR}" ]; then
  echo "Host VEP cache:     ${VEP_CACHE_DIR}"
  if [ ! -d "${VEP_CACHE_DIR}" ]; then
    echo "❌ ERROR: VEP cache dir not found: ${VEP_CACHE_DIR}"
    exit 1
  fi
else
  echo "❌ ERROR: --vep-cache-dir is required (cache must be mounted, not copied)."
  exit 1
fi
echo "==============================================================================="
echo ""
echo "⚠️  WARNING: This build includes VEP annotation and may take 30-60 minutes!"
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
# Exclude main data directory but keep tests/data
rsync -a --exclude='.git' --exclude='/data' --exclude='*.pyc' --exclude='__pycache__' \
  ./ "${BUILD_CONTEXT_DIR}/"

echo "  - Linking BCF files (hard links, no extra space)..."
# Create hard links to BCF (same filesystem, zero copy)
mkdir -p "${BUILD_CONTEXT_DIR}/docker/gnomad-data"
ln -f "${BCF_FILE}" "${BUILD_CONTEXT_DIR}/docker/gnomad-data/${BCF_BASENAME}"
ln -f "${BCF_FILE}.csi" "${BUILD_CONTEXT_DIR}/docker/gnomad-data/${BCF_BASENAME}.csi"

echo "🐳 Building Docker image (this will take a while)..."
START_TIME=$(date +%s)

# Enable BuildKit if available (optional, legacy builder works fine)
if docker buildx version &>/dev/null; then
  echo "  Using BuildKit (modern builder)"
  export DOCKER_BUILDKIT=1
else
  echo "  Using legacy builder (buildx not installed, this is fine)"
  unset DOCKER_BUILDKIT
fi

# Build from the temporary context on /mnt/data (Stage 1: Base image without annotation)
TEMP_IMAGE="${IMAGE_NAME}-temp"
docker build \
  -f "${BUILD_CONTEXT_DIR}/docker/Dockerfile.annotated" \
  --build-arg AF="${AF}" \
  --build-arg GENOME="${GENOME}" \
  --build-arg CACHE_NAME="${CACHE_NAME}" \
  --build-arg BCF_FILE="docker/gnomad-data/${BCF_BASENAME}" \
  --build-arg ANNOTATION_NAME="${ANNOTATION_NAME}" \
  -t "${TEMP_IMAGE}" \
  ${NO_CACHE} \
  ${HOST_NETWORK} \
  "${BUILD_CONTEXT_DIR}"

echo ""
echo "✅ Base image built successfully!"
echo ""

# Stage 2: Run annotation with VEP cache mounted
if [ -z "${VEP_CACHE_DIR}" ]; then
  echo "❌ ERROR: --vep-cache-dir is required for annotated build"
  exit 1
fi

echo "🧬 Running VEP annotation with mounted cache..."
echo "   VEP cache: ${VEP_CACHE_DIR}"
echo "   (This may take 30-60 minutes depending on dataset size)"
echo ""

# Run container with VEP cache mounted and execute annotation script
# Run as root so we can write to /cache directory
docker run \
  --name vcfstash-annotate-temp \
  --user root \
  -v "${VEP_CACHE_DIR}:/opt/vep/.vep:ro" \
  --entrypoint /bin/bash \
  "${TEMP_IMAGE}" \
  -c "export VCFSTASH_ROOT=/app && \
      export PATH=/usr/local/bin:/opt/vep/ensembl-vep:\$PATH && \
      bash /tmp/build-annotated-cache.sh \
        --bcf-file /tmp/gnomad.bcf \
        --gnomad-af '${AF}' \
        --cache-dir '/cache' \
        --threads 8 \
        --cache-name '${CACHE_NAME}' \
        --genome '${GENOME}' \
        --params /app/recipes/docker-annotated/params.yaml \
        --annotation-config /app/recipes/docker-annotated/annotation.config \
        --annotation-name '${ANNOTATION_NAME}' \
        --vep-cache '/opt/vep/.vep'"

# Check if annotation succeeded
if [ $? -ne 0 ]; then
  echo "❌ ERROR: Annotation failed!"
  docker rm vcfstash-annotate-temp
  exit 1
fi

echo ""
echo "✅ Annotation completed successfully!"
echo ""
echo "📦 Committing container to final image..."

# Commit the container with annotation results to final image
docker commit vcfstash-annotate-temp "${IMAGE_NAME}"
docker tag "${IMAGE_NAME}" "${REGISTRY}/vcfstash-annotated:latest"

# Clean up temporary container and image
docker rm vcfstash-annotate-temp
docker rmi "${TEMP_IMAGE}"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
DURATION_MIN=$((DURATION / 60))

echo ""
echo "==============================================================================="
echo "✅ Annotated image built successfully!"
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

# Verify cache exists
echo "🧪 Verifying annotated cache..."
docker run --rm --entrypoint ls "${IMAGE_NAME}" -la /cache/db/stash/${ANNOTATION_NAME}/

echo ""
echo "✅ Image verification passed!"
echo ""

# Push if requested
if [ "${PUSH}" = true ]; then
  echo "📤 Pushing to registry..."
  docker push "${IMAGE_NAME}"
  docker push "${REGISTRY}/vcfstash-annotated:latest"
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
echo "To test the annotated image:"
echo "  docker run --rm ${IMAGE_NAME} --help"
echo ""
echo "To run tests (auto-detects annotated scenario):"
echo "  docker run --rm --entrypoint /bin/sh ${IMAGE_NAME} \\"
echo "    -c 'cd /app && export PYTHONPATH=/app/venv/lib/python3.13/site-packages && python3 -m pytest tests/ -v'"
echo ""
echo "To annotate a sample:"
echo "  docker run --rm \\"
echo "    -v /path/to/sample.vcf:/data/sample.vcf \\"
echo "    -v /path/to/output:/output \\"
echo "    ${IMAGE_NAME} \\"
echo "    annotate \\"
echo "      -a /cache/db/stash/${ANNOTATION_NAME} \\"
echo "      --vcf /data/sample.vcf \\"
echo "      --output /output \\"
echo "      -y /app/recipes/docker-annotated/params.yaml"
echo ""
if [ "${PUSH}" = false ]; then
  echo "To push to registry:"
  echo "  docker push ${IMAGE_NAME}"
  echo ""
fi
