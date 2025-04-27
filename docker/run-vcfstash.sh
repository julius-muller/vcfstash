#!/bin/bash
# run-vcfstash.sh - A wrapper script for running VCFstash with Docker Compose

# Set default directories or use environment variables if already set
export REFERENCE_DIR=${REFERENCE_DIR:-"$(pwd)/reference"}
export DATA_DIR=${DATA_DIR:-"$(pwd)/data"}
export CACHE_DIR=${CACHE_DIR:-"$(pwd)/cache"}

# Ensure directories exist
mkdir -p "$REFERENCE_DIR" "$DATA_DIR" "$CACHE_DIR"

# Print configuration
echo "VCFstash Docker configuration:"
echo "  Reference directory: $REFERENCE_DIR"
echo "  Data directory: $DATA_DIR"
echo "  Cache directory: $CACHE_DIR"
echo ""

# Run docker-compose with all arguments passed to this script
cd "$(dirname "$0")" # Change to the script directory
docker-compose run --rm vcfstash "$@"