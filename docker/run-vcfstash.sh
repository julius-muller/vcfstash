#!/bin/bash
set -e

# Default directories
DEFAULT_DATA_DIR="./data"
DEFAULT_CACHE_DIR="./cache"

# Use environment variables if set, otherwise use defaults
DATA_DIR=${DATA_DIR:-$DEFAULT_DATA_DIR}
CACHE_DIR=${CACHE_DIR:-$DEFAULT_CACHE_DIR}

# Export variables for docker-compose
export DATA_DIR
export CACHE_DIR

# Run docker with the provided arguments
docker run --rm \
  -v ${DATA_DIR}:/data \
  -v ${CACHE_DIR}:/cache \
  vcfcache:latest "$@"
