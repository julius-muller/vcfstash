#!/usr/bin/env bash
#
# Build and push final production Docker images for VCFstash
#
# This script:
# 1. Builds blueprint images for AF 0.10 and AF 0.01
# 2. Builds base annotated images
# 3. Annotates with VEP and commits
# 4. Runs tests on each image
# 5. Pushes to registry if tests pass
#
# Usage: ./build-and-push-final.sh [--skip-tests] [--skip-push] [--af010-only] [--af001-only]

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Data paths
GNOMAD_DIR="/mnt/data/vcfstash_data/gnomad"
VEP_CACHE_DIR="/mnt/data/apps/ensembl-vep/115/cachedir"

# gnomAD files
GNOMAD_AF010="${GNOMAD_DIR}/gnomad_v4.1_GRCh38_joint_af010.bcf"
GNOMAD_AF001="${GNOMAD_DIR}/gnomad_v4.1_GRCh38_joint_af001.bcf"

# Build flags
SKIP_TESTS=false
SKIP_PUSH=false
BUILD_AF010=true
BUILD_AF001=true
FORCE_BUILD=false

# ============================================================================
# Parse arguments
# ============================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-tests)
            SKIP_TESTS=true
            shift
            ;;
        --skip-push)
            SKIP_PUSH=true
            shift
            ;;
        --force)
            FORCE_BUILD=true
            shift
            ;;
        --af010-only)
            BUILD_AF001=false
            shift
            ;;
        --af001-only)
            BUILD_AF010=false
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --skip-tests    Skip running tests after build"
            echo "  --skip-push     Skip pushing images to registry"
            echo "  --force         Force rebuild even if images exist"
            echo "  --af010-only    Build only AF ≥ 0.10 images"
            echo "  --af001-only    Build only AF ≥ 0.01 images"
            echo "  -h, --help      Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# ============================================================================
# Helper functions
# ============================================================================

log_info() {
    echo ""
    echo "═══════════════════════════════════════════════════════════════════════"
    echo "  $1"
    echo "═══════════════════════════════════════════════════════════════════════"
    echo ""
}

log_step() {
    echo ""
    echo "───────────────────────────────────────────────────────────────────────"
    echo "  ▸ $1"
    echo "───────────────────────────────────────────────────────────────────────"
}

log_success() {
    echo ""
    echo "✓ SUCCESS: $1"
    echo ""
}

log_error() {
    echo ""
    echo "✗ ERROR: $1" >&2
    echo ""
}

run_tests() {
    local image=$1
    local test_name=$2

    log_step "Running tests for $test_name"

    # Check if VEP cache is needed
    if [[ "$image" == *"annotated"* ]]; then
        log_info "Testing annotated image (requires VEP cache)"
        docker run --rm \
            -v "${VEP_CACHE_DIR}:/opt/vep/.vep:ro" \
            --entrypoint /bin/bash \
            "$image" \
            -c "cd /app && python3.13 -m pytest tests/ -v --color=yes"
    else
        log_info "Testing blueprint image (no VEP cache needed)"
        docker run --rm \
            --entrypoint /bin/bash \
            "$image" \
            -c "cd /app && python3.13 -m pytest tests/ -v --color=yes"
    fi

    if [[ $? -eq 0 ]]; then
        log_success "Tests passed for $test_name"
        return 0
    else
        log_error "Tests failed for $test_name"
        return 1
    fi
}

image_exists() {
    local image=$1
    docker image inspect "$image" &>/dev/null
}

push_image() {
    local image=$1
    local description=$2

    log_step "Pushing $description: $image"
    docker push "$image"

    if [[ $? -eq 0 ]]; then
        log_success "Pushed $description"
        return 0
    else
        log_error "Failed to push $description"
        return 1
    fi
}

# ============================================================================
# Build AF 0.10 (10%) images
# ============================================================================

build_af010() {
    log_info "Building AF ≥ 0.10 image set"

    # Define expected image names (match what the scripts actually generate)
    BLUEPRINT_IMAGE="ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af010"
    BASE_IMAGE="ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af010-vep115-base"
    ANNOTATED_IMAGE="ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af010-vep115"

    # Check if images already exist
    if [[ "$FORCE_BUILD" == "false" ]]; then
        if image_exists "$ANNOTATED_IMAGE"; then
            log_info "Annotated image already exists: $ANNOTATED_IMAGE"
            log_info "Skipping build. Use --force to rebuild."
            return 0
        fi
    fi

    # 1. Build blueprint
    log_step "Step 1/3: Building blueprint image"
    if [[ "$FORCE_BUILD" == "false" ]] && image_exists "$BLUEPRINT_IMAGE"; then
        log_info "Blueprint image already exists, skipping build"
    else
        if ! "$SCRIPT_DIR/03-build-blueprint.sh" \
            "$GNOMAD_AF010" \
            --host-network \
            --genome GRCh38 \
            --type joint; then
            log_error "Failed to build blueprint image"
            return 1
        fi
    fi

    # Test blueprint
    if [[ "$SKIP_TESTS" == "false" ]]; then
        run_tests "$BLUEPRINT_IMAGE" "AF 0.10 blueprint" || return 1
    fi

    # Push blueprint
    if [[ "$SKIP_PUSH" == "false" ]]; then
        push_image "$BLUEPRINT_IMAGE" "AF 0.10 blueprint" || return 1
    fi

    # 2. Build base image
    log_step "Step 2/3: Building base annotated image"
    if [[ "$FORCE_BUILD" == "false" ]] && image_exists "$BASE_IMAGE"; then
        log_info "Base image already exists, skipping build"
    else
        if ! "$SCRIPT_DIR/04a-build-base-image.sh" \
            "$GNOMAD_AF010" \
            --host-network \
            -y; then
            log_error "Failed to build base annotated image"
            return 1
        fi
    fi

    # 3. Annotate with VEP and commit
    log_step "Step 3/3: Annotating with VEP"
    if [[ "$FORCE_BUILD" == "false" ]] && image_exists "$ANNOTATED_IMAGE"; then
        log_info "Annotated image already exists, skipping annotation"
    else
        if ! "$SCRIPT_DIR/04b-annotate-and-commit.sh" \
            --base-image "$BASE_IMAGE" \
            --vep-cache-dir "$VEP_CACHE_DIR" \
            -y; then
            log_error "Failed to annotate with VEP"
            return 1
        fi
    fi

    # Test annotated image
    if [[ "$SKIP_TESTS" == "false" ]]; then
        run_tests "$ANNOTATED_IMAGE" "AF 0.10 annotated" || return 1
    fi

    # Push annotated image
    if [[ "$SKIP_PUSH" == "false" ]]; then
        push_image "$ANNOTATED_IMAGE" "AF 0.10 annotated" || return 1

        # Tag as latest
        log_step "Tagging AF 0.10 as latest"
        docker tag "$ANNOTATED_IMAGE" "ghcr.io/julius-muller/vcfstash-annotated:latest"
        push_image "ghcr.io/julius-muller/vcfstash-annotated:latest" "latest tag"
    fi

    log_success "AF 0.10 image set complete"
}

# ============================================================================
# Build AF 0.01 (1%) images
# ============================================================================

build_af001() {
    log_info "Building AF ≥ 0.01 image set"

    # Define expected image names (match what the scripts actually generate)
    BLUEPRINT_IMAGE="ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af001"
    BASE_IMAGE="ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af001-vep115-base"
    ANNOTATED_IMAGE="ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af001-vep115"

    # Check if images already exist
    if [[ "$FORCE_BUILD" == "false" ]]; then
        if image_exists "$ANNOTATED_IMAGE"; then
            log_info "Annotated image already exists: $ANNOTATED_IMAGE"
            log_info "Skipping build. Use --force to rebuild."
            return 0
        fi
    fi

    # 1. Build blueprint
    log_step "Step 1/3: Building blueprint image"
    if [[ "$FORCE_BUILD" == "false" ]] && image_exists "$BLUEPRINT_IMAGE"; then
        log_info "Blueprint image already exists, skipping build"
    else
        if ! "$SCRIPT_DIR/03-build-blueprint.sh" \
            "$GNOMAD_AF001" \
            --host-network \
            --genome GRCh38 \
            --type joint; then
            log_error "Failed to build blueprint image"
            return 1
        fi
    fi

    # Test blueprint
    if [[ "$SKIP_TESTS" == "false" ]]; then
        run_tests "$BLUEPRINT_IMAGE" "AF 0.01 blueprint" || return 1
    fi

    # Push blueprint
    if [[ "$SKIP_PUSH" == "false" ]]; then
        push_image "$BLUEPRINT_IMAGE" "AF 0.01 blueprint" || return 1
    fi

    # 2. Build base image
    log_step "Step 2/3: Building base annotated image"
    if [[ "$FORCE_BUILD" == "false" ]] && image_exists "$BASE_IMAGE"; then
        log_info "Base image already exists, skipping build"
    else
        if ! "$SCRIPT_DIR/04a-build-base-image.sh" \
            "$GNOMAD_AF001" \
            --host-network \
            -y; then
            log_error "Failed to build base annotated image"
            return 1
        fi
    fi

    # 3. Annotate with VEP and commit
    log_step "Step 3/3: Annotating with VEP"
    if [[ "$FORCE_BUILD" == "false" ]] && image_exists "$ANNOTATED_IMAGE"; then
        log_info "Annotated image already exists, skipping annotation"
    else
        if ! "$SCRIPT_DIR/04b-annotate-and-commit.sh" \
            --base-image "$BASE_IMAGE" \
            --vep-cache-dir "$VEP_CACHE_DIR" \
            -y; then
            log_error "Failed to annotate with VEP"
            return 1
        fi
    fi

    # Test annotated image
    if [[ "$SKIP_TESTS" == "false" ]]; then
        run_tests "$ANNOTATED_IMAGE" "AF 0.01 annotated" || return 1
    fi

    # Push annotated image
    if [[ "$SKIP_PUSH" == "false" ]]; then
        push_image "$ANNOTATED_IMAGE" "AF 0.01 annotated" || return 1
    fi

    log_success "AF 0.01 image set complete"
}

# ============================================================================
# Main execution
# ============================================================================

main() {
    log_info "VCFstash Final Image Build & Push"

    echo "Configuration:"
    echo "  - Skip tests: $SKIP_TESTS"
    echo "  - Skip push: $SKIP_PUSH"
    echo "  - Force build: $FORCE_BUILD"
    echo "  - Build AF 0.10: $BUILD_AF010"
    echo "  - Build AF 0.01: $BUILD_AF001"
    echo ""

    # Verify required files exist
    if [[ "$BUILD_AF010" == "true" ]] && [[ ! -f "$GNOMAD_AF010" ]]; then
        log_error "gnomAD AF 0.10 file not found: $GNOMAD_AF010"
        exit 1
    fi

    if [[ "$BUILD_AF001" == "true" ]] && [[ ! -f "$GNOMAD_AF001" ]]; then
        log_error "gnomAD AF 0.01 file not found: $GNOMAD_AF001"
        exit 1
    fi

    if [[ ! -d "$VEP_CACHE_DIR" ]]; then
        log_error "VEP cache directory not found: $VEP_CACHE_DIR"
        exit 1
    fi

    # Build images
    if [[ "$BUILD_AF010" == "true" ]]; then
        build_af010 || {
            log_error "AF 0.10 build failed"
            exit 1
        }
    fi

    if [[ "$BUILD_AF001" == "true" ]]; then
        build_af001 || {
            log_error "AF 0.01 build failed"
            exit 1
        }
    fi

    # Summary
    log_info "Build Complete!"
    echo "Images built successfully:"
    echo ""

    if [[ "$BUILD_AF010" == "true" ]]; then
        echo "AF ≥ 0.10 (10%):"
        echo "  Blueprint:  ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af010"
        echo "  Annotated:  ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af010-vep115"
        echo "  Latest:     ghcr.io/julius-muller/vcfstash-annotated:latest"
        echo ""
    fi

    if [[ "$BUILD_AF001" == "true" ]]; then
        echo "AF ≥ 0.01 (1%):"
        echo "  Blueprint:  ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af001"
        echo "  Annotated:  ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af001-vep115"
        echo ""
    fi

    if [[ "$SKIP_PUSH" == "true" ]]; then
        echo "NOTE: Images were not pushed (--skip-push)"
        echo "To push manually:"
        echo "  docker push ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af010"
        echo "  docker push ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af010-vep115"
        echo ""
    fi
}

# Run main
main
