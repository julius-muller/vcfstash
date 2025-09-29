#!/usr/bin/env bash
set -euo pipefail

# --- argument parser -------------------------------------------------------
AF="0.10"
CACHE_DIR="/home/micromamba/cache"
THREADS="8"
TOOL_VER="115.1"
CNAME="vep_gnomad"
GENOME="GRCh38"
GNOMAD_URL="${GNOMAD_URL:-}"   # will come from Docker ARG/ENV, may be overridden here

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gnomad-af)    AF="$2"; shift 2;;
    --cache-dir)    CACHE_DIR="$2"; shift 2;;
    --threads)      THREADS="$2"; shift 2;;
    --tool-version) TOOL_VER="$2"; shift 2;;
    --cache-name)   CNAME="$2"; shift 2;;
    --genome)       GENOME="$2"; shift 2;;
    --params)       PARAMS_FILE="$2"; shift 2;;
    --config)       CONFIG_FILE="$2"; shift 2;;
    --url|-u)       GNOMAD_URL="$2"; shift 2;;
    *) echo "Unknown arg $1"; exit 1;;
  esac
done

# --------------------------------------------------------------------------
# 1. choose source VCF  (we're inside the argument-parser block)

# BGZF-compressed file that actually exists (≈27 kB)
DEFAULT_URL="https://raw.githubusercontent.com/samtools/htslib/master/test/CEU-exon-realigned.vcf.bgz"

if [[ -z "${GNOMAD_URL}" ]]; then
  GNOMAD_URL="${DEFAULT_URL}"
fi
# --------------------------------------------------------------------------
# 2. fetch or synthesise
G_SRC="/tmp/gnomad.${GENOME}.vcf.bgz"

echo "Attempting download: ${GNOMAD_URL}"
if curl -fL "${GNOMAD_URL}" -o "${G_SRC}"; then
  echo "✓ downloaded test VCF"
else
  echo "× download failed – creating tiny inline VCF"
  cat > /tmp/toy.vcf <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=248956422>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM  POS     ID  REF ALT QUAL FILTER INFO
1       10000   .   G   A   .    PASS  AF=0.15
1       10500   .   C   T   .    PASS  AF=0.20
EOF
  bgzip -c /tmp/toy.vcf > "${G_SRC}"
fi

# ensure index
tabix -p vcf "${G_SRC}"

# --------------------------------------------------------------------------
# 3. Build blueprint & annotate
vcfstash stash-init \
        --vcf /tmp/gnomad_af.bcf \
        --output "${CACHE_DIR}" \
        -y "${PARAMS_FILE}"

vcfstash stash-annotate \
        --db    "${CACHE_DIR}" \
        --name  "${CNAME}" \
        -a      "${CONFIG_FILE}"

echo "Cache created in ${CACHE_DIR}"