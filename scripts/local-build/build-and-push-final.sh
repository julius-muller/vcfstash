#!/usr/bin/env bash
set -euo pipefail

# Build and (optionally) push the two code-only images:
#  - ghcr.io/julius-muller/vcfcache:latest          (lean CLI + bcftools)
#  - ghcr.io/julius-muller/vcfcache:vep115.2_basic  (CLI + bcftools + VEP + fixed recipe)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

SKIP_TESTS=false
SKIP_PUSH=false
FORCE_BUILD=false

while [[ $# -gt 0 ]]; do
  case $1 in
    --skip-tests) SKIP_TESTS=true; shift ;;
    --skip-push)  SKIP_PUSH=true; shift ;;
    --force)      FORCE_BUILD=true; shift ;;
    -h|--help)
      echo "Usage: $0 [--skip-tests] [--skip-push] [--force]"; exit 0 ;;
    *) echo "Unknown arg $1"; exit 1 ;;
  esac
done

log() { echo "[$(date +%H:%M:%S)] $*"; }

image_exists() { docker image inspect "$1" >/dev/null 2>&1; }

build_image() {
  local dockerfile=$1 tag=$2
  if image_exists "$tag" && ! $FORCE_BUILD; then
    log "Image $tag exists, skip (use --force to rebuild)"; return
  fi
  log "Building $tag from $dockerfile"
  (cd "$PROJECT_ROOT" && docker build -f "$dockerfile" -t "$tag" .)
}

test_image() {
  local tag=$1
  log "Testing $tag"
  docker run --rm "$tag" python -m pytest tests/test_archive.py tests/test_manifest.py tests/test_cli_alias_and_pull.py tests/test_integration_annotation.py -q
}

smoke_image() {
  local tag=$1
  log "Smoke testing $tag"
  docker run --rm "$tag" vcfcache --version
}

push_image() { log "Pushing $1"; docker push "$1"; }

# Build lean image
LEAN_TAG="ghcr.io/julius-muller/vcfcache:latest"
build_image docker/Dockerfile.vcfcache "$LEAN_TAG"
if ! $SKIP_TESTS; then test_image "$LEAN_TAG"; fi
if ! $SKIP_PUSH;  then push_image "$LEAN_TAG"; fi

# Build preset
PRESET_TAG="ghcr.io/julius-muller/vcfcache:vep115.2_basic"
build_image docker/Dockerfile.vep115.2_basic "$PRESET_TAG"
if ! $SKIP_TESTS; then smoke_image "$PRESET_TAG"; fi
if ! $SKIP_PUSH;  then push_image "$PRESET_TAG"; fi

log "Done. Built: $LEAN_TAG and $PRESET_TAG"
