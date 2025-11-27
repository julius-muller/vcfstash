#!/usr/bin/env bash
set -euo pipefail

#===============================================================================
# Push Docker Images to Registry
#===============================================================================
# This script pushes locally-built Docker images to GitHub Container Registry.
#
# Usage:
#   ./05-push-images.sh [OPTIONS]
#
# Options:
#   --type TYPE             Image type: blueprint/annotated/both (default: both)
#   --tag TAG               Specific tag to push (default: all local tags)
#   --registry REGISTRY     Docker registry (default: ghcr.io/julius-muller)
#   --dry-run               Show what would be pushed without pushing
#===============================================================================

# Default configuration
TYPE="both"
TAG=""
REGISTRY="ghcr.io/julius-muller"
DRY_RUN=false

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --type)      TYPE="$2"; shift 2;;
    --tag)       TAG="$2"; shift 2;;
    --registry)  REGISTRY="$2"; shift 2;;
    --dry-run)   DRY_RUN=true; shift;;
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

# Validate type
if [[ ! "${TYPE}" =~ ^(blueprint|annotated|both)$ ]]; then
  echo "❌ ERROR: Invalid type: ${TYPE}"
  echo "Must be: blueprint, annotated, or both"
  exit 1
fi

echo "==============================================================================="
echo "Push Docker Images to Registry"
echo "==============================================================================="
echo "Registry:    ${REGISTRY}"
echo "Type:        ${TYPE}"
if [ -n "${TAG}" ]; then
  echo "Tag:         ${TAG}"
else
  echo "Tag:         all local tags"
fi
echo "Dry run:     ${DRY_RUN}"
echo "==============================================================================="
echo ""

# Function to push image
push_image() {
  local image=$1

  if [ "${DRY_RUN}" = true ]; then
    echo "  [DRY RUN] Would push: ${image}"
  else
    echo "  Pushing: ${image}"
    docker push "${image}"
    echo "  ✅ Pushed successfully"
  fi
}

# Function to list and push images
push_images_by_type() {
  local img_type=$1
  local img_name="${REGISTRY}/vcfstash-${img_type}"

  echo "───────────────────────────────────────────────────────────────────────────"
  echo "Processing ${img_type} images"
  echo "───────────────────────────────────────────────────────────────────────────"

  # Get all matching images
  if [ -n "${TAG}" ]; then
    # Specific tag
    IMAGES=$(docker images "${img_name}:${TAG}" --format "{{.Repository}}:{{.Tag}}" | grep -v "<none>")
  else
    # All tags
    IMAGES=$(docker images "${img_name}" --format "{{.Repository}}:{{.Tag}}" | grep -v "<none>")
  fi

  if [ -z "${IMAGES}" ]; then
    echo "⚠️  No ${img_type} images found"
    echo ""
    return
  fi

  echo "Found images:"
  echo "${IMAGES}" | sed 's/^/  - /'
  echo ""

  # Push each image
  while IFS= read -r image; do
    push_image "${image}"
  done <<< "${IMAGES}"

  echo ""
}

# Check if logged in to registry
if ! docker info 2>/dev/null | grep -q "${REGISTRY}"; then
  echo "⚠️  Not logged in to ${REGISTRY}"
  echo ""
  echo "Log in with:"
  echo "  echo \$GITHUB_TOKEN | docker login ghcr.io -u USERNAME --password-stdin"
  echo ""
  read -p "Continue anyway? (y/N): " -n 1 -r
  echo
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 1
  fi
fi

# Push images based on type
if [ "${TYPE}" = "blueprint" ] || [ "${TYPE}" = "both" ]; then
  push_images_by_type "blueprint"
fi

if [ "${TYPE}" = "annotated" ] || [ "${TYPE}" = "both" ]; then
  push_images_by_type "annotated"
fi

echo "==============================================================================="
echo "Done!"
echo "==============================================================================="
echo ""

if [ "${DRY_RUN}" = true ]; then
  echo "This was a dry run. To actually push, remove --dry-run flag."
  echo ""
fi

echo "View your images at:"
echo "  https://github.com/users/julius-muller/packages?repo_name=vcfstash"
echo ""
