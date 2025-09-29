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
# 1. create minimal reference genome if needed
REF_PATH="${CACHE_DIR}/reference.fa"
if [[ ! -f "${REF_PATH}" ]]; then
    echo "Creating minimal reference genome for testing..."
    # Create a minimal chromosome 1 with enough bases for our test variants
    echo ">1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF" > "${REF_PATH}"
    # Generate 20000 bases of 'A' to cover our test variant positions (10001-15000)
    python3 -c "print('A' * 20000)" >> "${REF_PATH}"
fi

# Update params.yaml to use our temporary reference
TEMP_PARAMS="${CACHE_DIR}/params_temp.yaml"
sed "s|reference: \"/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa\"|reference: \"${REF_PATH}\"|" "${PARAMS_FILE}" > "${TEMP_PARAMS}"
PARAMS_FILE="${TEMP_PARAMS}"

# --------------------------------------------------------------------------
# 2. obtain BGZF-indexed VCF  (download or synthesize)

WORK_TMP="${CACHE_DIR}/_tmp"
mkdir -p "${WORK_TMP}"
G_SRC="${WORK_TMP}/source.vcf.bgz"

echo "Attempting download: ${GNOMAD_URL:-<none>}"
if [[ -n "${GNOMAD_URL}" ]] && curl -fL "${GNOMAD_URL}" -o "${G_SRC}"; then
    echo "✓ downloaded test VCF"
else
    echo "× download failed – generating inline 2-variant BGZF VCF"
    cat > /tmp/toy.vcf <<'HEADER'
##fileformat=VCFv4.2
##contig=<ID=1,length=248956422>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
HEADER

    # create 5 000 fake single-alt SNPs; AF alternates between 0.11 and 0.12
    for i in $(seq 1 5000); do
        pos=$((10000 + i))
        af=$(printf "0.%02d" $((10 + i % 2)))   # 0.11 / 0.12
        echo -e "1\t${pos}\t.\tA\tG\t.\tPASS\tAF=${af}" >> /tmp/toy.vcf
    done

    bgzip -c /tmp/toy.vcf > "${G_SRC}"
fi

tabix -f -p vcf "${G_SRC}"

# --- build blueprint & annotate -------------------------------------------
DB_DIR="${CACHE_DIR}/db"
rm -rf "${DB_DIR}"

vcfstash stash-init   --force \
        --vcf "${G_SRC}" \
        --output "${DB_DIR}" \
        -y "${PARAMS_FILE}"

vcfstash stash-annotate \
        --db   "${DB_DIR}" \
        --name "${CNAME}" \
        -a     "${CONFIG_FILE}"