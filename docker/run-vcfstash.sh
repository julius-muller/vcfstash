#!/bin/bash
set -e

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Default directories
DEFAULT_REFERENCE_DIR="./reference"
DEFAULT_DATA_DIR="./data"
DEFAULT_CACHE_DIR="./cache"

# Use environment variables if set, otherwise use defaults
REFERENCE_DIR=${REFERENCE_DIR:-$DEFAULT_REFERENCE_DIR}
DATA_DIR=${DATA_DIR:-$DEFAULT_DATA_DIR}
CACHE_DIR=${CACHE_DIR:-$DEFAULT_CACHE_DIR}

# Export variables for docker-compose
export REFERENCE_DIR
export DATA_DIR
export CACHE_DIR

# Run docker-compose with the provided arguments
cd "$SCRIPT_DIR/.."
docker-compose -f docker/docker-compose.yml run --rm vcfstash "$@"