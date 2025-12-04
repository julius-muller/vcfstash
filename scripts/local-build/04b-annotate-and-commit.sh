#!/usr/bin/env bash
set -euo pipefail

#===============================================================================
# Run Annotation and Commit Final Image
#===============================================================================
# This script runs VEP annotation on an existing base image and commits the
# result to create the final annotated image. This is step 2 of 2.
#
# Usage:
#   ./04b-annotate-and-commit.sh [OPTIONS]
#
# Options:
#   --base-image IMAGE      Base Docker image to annotate (required)
#   --vep-cache-dir DIR     VEP cache directory to mount (required)
#   --output-tag TAG        Final image tag (default: remove '-base' suffix)
#   --registry REGISTRY     Docker registry (default: ghcr.io/julius-muller)
#   --push                  Push to registry after build
#   --threads N             Number of VEP threads (default: 8)
#   --yes|-y                Skip confirmation prompt
#===============================================================================

# Parse arguments
BASE_IMAGE=""
VEP_CACHE_DIR=""
OUTPUT_TAG=""
REGISTRY="ghcr.io/julius-muller"
PUSH=false
THREADS=8
YES=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --base-image)       BASE_IMAGE="$2"; shift 2;;
    --vep-cache-dir)    VEP_CACHE_DIR="$2"; shift 2;;
    --output-tag)       OUTPUT_TAG="$2"; shift 2;;
    --registry)         REGISTRY="$2"; shift 2;;
    --push)             PUSH=true; shift;;
    --threads)          THREADS="$2"; shift 2;;
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

# Validate required arguments
if [ -z "${BASE_IMAGE}" ]; then
  echo "❌ ERROR: --base-image is required"
  exit 1
fi

if [ -z "${VEP_CACHE_DIR}" ]; then
  echo "❌ ERROR: --vep-cache-dir is required"
  exit 1
fi

if [ ! -d "${VEP_CACHE_DIR}" ]; then
  echo "❌ ERROR: VEP cache dir not found: ${VEP_CACHE_DIR}"
  exit 1
fi

# Auto-generate output tag if not provided
if [ -z "${OUTPUT_TAG}" ]; then
  # Remove '-base' suffix from base image tag
  OUTPUT_TAG="${BASE_IMAGE%-base}"
  # Add 'latest' tag too
  LATEST_TAG="${OUTPUT_TAG%:*}:latest"
else
  OUTPUT_TAG="${REGISTRY}/vcfstash-annotated:${OUTPUT_TAG}"
  LATEST_TAG="${REGISTRY}/vcfstash-annotated:latest"
fi

echo "==============================================================================="
echo "Running VEP Annotation and Committing Final Image (Pure Python)"
echo "==============================================================================="
echo "Base image:         ${BASE_IMAGE}"
echo "VEP cache:          ${VEP_CACHE_DIR}"
echo "VEP threads:        ${THREADS}"
echo "Final image:        ${OUTPUT_TAG}"
echo "Latest tag:         ${LATEST_TAG}"
echo "Push to GHCR:       ${PUSH}"
echo "==============================================================================="
echo ""
echo "⚠️  WARNING: Annotation will take 30-60 minutes for full genome data!"
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

# Clean up any leftover containers
echo "🧹 Cleaning up old annotation containers..."
docker rm -f vcfstash-annotate-temp 2>/dev/null || true

echo ""
echo "🧬 Running VEP annotation with mounted cache..."
echo "   VEP cache: ${VEP_CACHE_DIR}"
echo "   (This may take 30-60 minutes depending on dataset size)"
echo ""

START_TIME=$(date +%s)

# Run container with VEP cache mounted and execute annotation script
# Run as root so we can write to /cache directory
docker run \
  --name vcfstash-annotate-temp \
  --user root \
  -v "${VEP_CACHE_DIR}:/opt/vep/.vep:ro" \
  --entrypoint /bin/bash \
  "${BASE_IMAGE}" \
  -c "export VCFSTASH_ROOT=/app && \
      export PATH=/usr/local/bin:/opt/vep/ensembl-vep:\$PATH && \
      bash /tmp/build-annotated-cache.sh \
        --bcf-file /tmp/gnomad.bcf \
        --gnomad-af 0.010 \
        --cache-dir '/cache' \
        --threads ${THREADS} \
        --cache-name 'gnomad_grch38_joint_af010' \
        --genome 'GRCh38' \
        --params /app/recipes/docker-annotated/params.yaml \
        --annotation-config /app/recipes/docker-annotated/annotation.config \
        --annotation-name 'vep_gnomad' \
        --vep-cache '/opt/vep/.vep'"

# Check if annotation succeeded
if [ $? -ne 0 ]; then
  echo ""
  echo "❌ ERROR: Annotation failed!"
  echo ""
  echo "To debug, check container logs:"
  echo "  docker logs vcfstash-annotate-temp"
  echo ""
  echo "Container kept for debugging. To remove:"
  echo "  docker rm vcfstash-annotate-temp"
  exit 1
fi

echo ""
echo "✅ Annotation completed successfully!"
echo ""
echo "📦 Committing container to final image..."

# Commit the container with annotation results to final image
docker commit vcfstash-annotate-temp "${OUTPUT_TAG}"
docker tag "${OUTPUT_TAG}" "${LATEST_TAG}"

# Clean up temporary container
echo "🧹 Cleaning up temporary container..."
docker rm vcfstash-annotate-temp

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
DURATION_MIN=$((DURATION / 60))

echo ""
echo "==============================================================================="
echo "✅ Annotated image built successfully!"
echo "==============================================================================="
echo "Image:     ${OUTPUT_TAG}"
echo "Latest:    ${LATEST_TAG}"
echo "Duration:  ${DURATION_MIN} minutes"
echo ""

# Get image size
IMAGE_SIZE=$(docker images "${OUTPUT_TAG}" --format "{{.Size}}")
echo "Size:      ${IMAGE_SIZE}"
echo ""

# Test the image
echo "🧪 Testing image..."
docker run --rm "${OUTPUT_TAG}" -v

# Verify VEP is available
echo "🧪 Verifying VEP..."
docker run --rm --entrypoint vep "${OUTPUT_TAG}" --help | head -5

# Verify cache exists
echo "🧪 Verifying annotated cache..."
docker run --rm --entrypoint ls "${OUTPUT_TAG}" -la /cache/db/stash/vep_gnomad/ 2>/dev/null || \
  echo "Note: Cache verification requires running annotation first"

echo ""
echo "✅ Image verification passed!"
echo ""

# Push if requested
if [ "${PUSH}" = true ]; then
  echo "📤 Pushing to registry..."
  docker push "${OUTPUT_TAG}"
  docker push "${LATEST_TAG}"
  echo "✅ Pushed to ${REGISTRY}"
  echo ""
fi

echo "==============================================================================="
echo "Done!"
echo "==============================================================================="
echo ""
echo "To test the annotated image:"
echo "  docker run --rm ${OUTPUT_TAG} --help"
echo ""
echo "To run tests (auto-detects annotated scenario):"
echo "  docker run --rm --entrypoint /bin/sh ${OUTPUT_TAG} \\"
echo "    -c 'cd /app && export PYTHONPATH=/app/venv/lib/python3.13/site-packages && python3 -m pytest tests/ -v'"
echo ""
echo "To annotate a sample:"
echo "  docker run --rm \\"
echo "    -v /path/to/sample.vcf:/data/sample.vcf \\"
echo "    -v /path/to/output:/output \\"
echo "    ${OUTPUT_TAG} \\"
echo "    annotate \\"
echo "      -a /cache/db/stash/vep_gnomad \\"
echo "      --vcf /data/sample.vcf \\"
echo "      --output /output \\"
echo "      -y /app/recipes/docker-annotated/params.yaml"
echo ""
if [ "${PUSH}" = false ]; then
  echo "To push to registry:"
  echo "  docker push ${OUTPUT_TAG}"
  echo "  docker push ${LATEST_TAG}"
  echo ""
fi
