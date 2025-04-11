# scripts/entrypoint.sh
#!/bin/bash
set -e

# Ensure directories exist
mkdir -p /cache /data

# Check if command is provided
if [ "$1" = "" ]; then
    echo "Usage: vcfstash <command> [options]"
    exit 1
fi

# Set VCF environment variables
export VCF_PATH=/opt/vcf
export VCF_DATA=/opt/vcf/.vcf

# Execute command
exec "$@"