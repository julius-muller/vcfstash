#!/usr/bin/env bash
set -euo pipefail

# Build and (optionally) push the runtime image:
#  - ghcr.io/julius-muller/vcfcache:latest (CLI + bcftools)

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
  (cd "$PROJECT_ROOT" && docker build -f "$dockerfile" -t "$tag" --network=host .)
}

test_image() {
  local dockerfile=$1 tag=$2
  log "Testing $tag (using test stage)"
  # Build and run the test stage which has dev dependencies and runs pytest
  docker build -f "$dockerfile" --target test -t "$tag-test" --network=host "$PROJECT_ROOT"
  docker run --rm "$tag-test"
}

smoke_image() {
  local tag=$1
  log "Smoke testing $tag"
  docker run --rm "$tag" --version
}

push_image() { log "Pushing $1"; docker push "$1"; }

LEAN_TAG="ghcr.io/julius-muller/vcfcache:latest"
build_image docker/Dockerfile.vcfcache "$LEAN_TAG"
smoke_image "$LEAN_TAG"
if ! $SKIP_TESTS; then test_image docker/Dockerfile.vcfcache "$LEAN_TAG"; fi
if ! $SKIP_PUSH;  then push_image "$LEAN_TAG"; fi

log "Done. Built: $LEAN_TAG"
