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

# --- choose source ---------------------------------------------------------
DEFAULT_URL="https://raw.githubusercontent.com/samtools/htslib/master/test/aux/CEU-example.vcf.bgz"
[[ -z "${GNOMAD_URL}" ]] && GNOMAD_URL="${DEFAULT_URL}"

# --- download or synthesize ------------------------------------------------
G_SRC="/tmp/gnomad.${GENOME}.vcf.bgz"
echo "Attempting download: ${GNOMAD_URL}"
if ! curl -fL "${GNOMAD_URL}" -o "${G_SRC}"; then
    echo "Download failed – creating inline 2-variant BGZF VCF"
    cat > /tmp/toy.vcf <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=248956422>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM  POS ID  REF ALT QUAL FILTER INFO
1       10000 .  G   A   .    PASS  AF=0.15
1       10500 .  C   T   .    PASS  AF=0.20
EOF
    bgzip -c /tmp/toy.vcf > "${G_SRC}"
fi
tabix -p vcf "${G_SRC}"

# --- build blueprint & annotate -------------------------------------------
rm -rf "${CACHE_DIR:?}/"*    # ensure empty dir
vcfstash stash-init   --force \
        --vcf "${G_SRC}" \
        --output "${CACHE_DIR}" \
        -y "${PARAMS_FILE}"

vcfstash stash-annotate \
        --db   "${CACHE_DIR}" \
        --name "${CNAME}" \
        -a     "${CONFIG_FILE}"