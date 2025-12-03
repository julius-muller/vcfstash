#!/usr/bin/env bash
set -euo pipefail

# Benchmark cached vs uncached annotation on down-sampled subsets of a large BCF.
#
# Scales:
#   PANEL  - 5k variants (approx. targeted panel)
#   WES    - 100k variants (approx. exome)
#   WGS    - all variants from source (full WGS)
#
# Requirements on host:
#   - bcftools, tabix, shuf
#   - docker with access to:
#       ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af010-vep115
#       ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af001-vep115
#   - VEP cache directory (set VEP_CACHE_DIR)
#
# Usage:
#   ./tests/run_benchmarks_subset.sh [SOURCE_BCF]
#   SOURCE_BCF defaults to /mnt/data/samples/test_mgm/mgm_WGS_32.gatkWGS_norm_hg38.bcf
#
# Output:
#   Logs are written to ./tests/benchmarks/ in TSV format:
#     timestamp\timage\tmode\tscale\tvariants\tseconds\tstatus\toutput_bcf

SOURCE_BCF=${1:-/mnt/data/samples/test_mgm/mgm_WGS_32.gatkWGS_norm_hg38.bcf}
VEP_CACHE_DIR=${VEP_CACHE_DIR:-/mnt/data/apps/ensembl-vep/115/cachedir}

if [[ ! -f "$SOURCE_BCF" ]]; then
  echo "Source BCF not found: $SOURCE_BCF" >&2
  exit 1
fi
if [[ ! -f "${SOURCE_BCF}.csi" ]]; then
  echo "Index not found for $SOURCE_BCF (.csi required)" >&2
  exit 1
fi
if [[ ! -d "$VEP_CACHE_DIR" ]]; then
  echo "VEP cache dir not found: $VEP_CACHE_DIR" >&2
  exit 1
fi

SCALES=("PANEL:5000" "WES:100000" "WGS:FULL")
IMAGES=(
  "ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af010-vep115"
  "ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af001-vep115"
)

LOG_DIR="$(pwd)/tests/benchmarks"
mkdir -p "$LOG_DIR"

tsv_log() {
  local line="$1"
  echo -e "$line" | tee -a "$LOG_FILE"
}

mk_subset() {
  local scale_name=$1
  local n=$2
  local dest_dir=$3

  # WGS/FULL: use source directly
  if [[ "$n" == "FULL" ]]; then
    echo "$SOURCE_BCF"
    return
  fi
  mkdir -p "$dest_dir"
  local list="$dest_dir/${scale_name}.sites"
  local out="$dest_dir/${scale_name}.bcf"

  # Random subset of sites
  bcftools view -H "$SOURCE_BCF" | shuf -n "$n" | cut -f1-2 > "$list"
  bcftools view -R "$list" -Ob -o "$out" "$SOURCE_BCF"
  bcftools index "$out"
  echo "$out"
}

run_bench() {
  local image=$1
  local mode=$2          # cached or uncached
  local scale=$3
  local bcf=$4
  local outdir=$5

  mkdir -p "$outdir"
  local bname
  bname="$(basename "$bcf")"
  local out_name="${bname%.bcf}_vst.bcf"
  local run_name="run_${mode:-cached}_${scale}"
  local run_dir_host="${outdir}/${run_name}"
  local run_dir_cont="/out/${run_name}"
  local outfile="${run_dir_host}/${out_name}"

  # Check if this benchmark was already completed (log file has at least 2 lines)
  local log_entry_count=0
  if [ -f "$LOG_FILE" ]; then
    log_entry_count=$(grep -c "^[0-9]" "$LOG_FILE" || true)
  fi

  # Check if this specific run is already in the log
  if [ -f "$LOG_FILE" ] && grep -q "${mode}.*${scale}" "$LOG_FILE" 2>/dev/null; then
    if [ "$log_entry_count" -ge 1 ]; then
      echo "Skipping $mode $scale - already completed"
      return 0
    else
      # Log exists but incomplete, delete and re-run
      rm -f "$LOG_FILE"
    fi
  fi

  # ensure fresh run dir on host (vcfstash will create it)
  rm -rf "${run_dir_host}"
  local start=$(date -u +%s)
  set +e

  # Get the project root directory (parent of tests/)
  local script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  local project_root="$(dirname "$script_dir")"

  docker run --rm \
    -v "$bcf":/work/input.bcf:ro \
    -v "${bcf}.csi":/work/input.bcf.csi:ro \
    -v "$VEP_CACHE_DIR":/opt/vep/.vep:ro \
    -v "$outdir":/out \
    -v /tmp:/tmp \
    -v "${project_root}/vcfstash":/app/venv/lib/python3.13/site-packages/vcfstash:ro \
    -w /app \
    --user "$(id -u):$(id -g)" \
    --entrypoint /bin/bash \
    "$image" \
    -lc "NXF_HOME=${run_dir_cont}/.nxf NXF_WORK=${run_dir_cont}/work \
         VCFSTASH_LOGLEVEL=ERROR VCFSTASH_FILE_LOGLEVEL=ERROR NXF_ANSI_LOG=false \
         vcfstash annotate ${mode} \
         --force \
         -a /cache/db/stash/vep_gnomad \
         --vcf /work/input.bcf \
         --output ${run_dir_cont} \
         -y /app/recipes/docker-annotated/params.yaml"
  status=$?
  set -e
  local end=$(date -u +%s)
  local elapsed=$((end - start))
  local outfile="${run_dir_host}/${out_name}"

  # Log the result
  tsv_log "$(date -Iseconds)\t${image}\t${mode}\t${scale}\t$(bcftools index -n "$bcf")\t${elapsed}\t${status}\t${outfile}"

  # Clean up intermediate files on success, keep only logs
  if [ "$status" -eq 0 ]; then
    rm -rf "${run_dir_host}/work" "${run_dir_host}/.nxf" 2>/dev/null || true
    find "${run_dir_host}" -name "*.bcf" -o -name "*.bcf.csi" -o -name "*.html" | xargs rm -f 2>/dev/null || true
  fi
}

main() {
  # Pre-flight summary
  echo "=== VCFstash Benchmark Plan ==="
  echo "Source: $SOURCE_BCF"
  echo ""
  echo "Scales to test:"
  for scale_def in "${SCALES[@]}"; do
    IFS=":" read -r scale_name nvars <<<"$scale_def"
    if [[ "$nvars" == "FULL" ]]; then
      echo "  - $scale_name (full source file)"
    else
      echo "  - $scale_name ($nvars variants)"
    fi
  done
  echo ""
  echo "Docker images:"
  for image in "${IMAGES[@]}"; do
    echo "  - ${image##*:}"
  done
  echo ""
  echo "Modes: cached, uncached"
  echo "Total benchmarks: $((${#SCALES[@]} * ${#IMAGES[@]} * 2))"
  echo "==============================="
  echo ""

  for scale_def in "${SCALES[@]}"; do
    IFS=":" read -r scale_name nvars <<<"$scale_def"
    bench_dir="$LOG_DIR/${scale_name,,}"
    mkdir -p "$bench_dir"

    # Use full source file for WGS (nvars='FULL'), create/reuse subset for others
    if [[ "$nvars" == "FULL" ]]; then
      subset_bcf="$SOURCE_BCF"
    else
      if [ -f "$bench_dir/subset/${scale_name}.bcf" ]; then
        subset_bcf="$bench_dir/subset/${scale_name}.bcf"
      else
        subset_bcf=$(mk_subset "$scale_name" "$nvars" "$bench_dir/subset")
      fi
    fi

    for image in "${IMAGES[@]}"; do
      LOG_FILE="$bench_dir/${image##*:}_${scale_name}.log"
      # Only write header if log file doesn't exist
      if [ ! -f "$LOG_FILE" ]; then
        tsv_log "timestamp\timage\tmode\tscale\tvariants\tseconds\tstatus\toutput_bcf"
      fi
      run_bench "$image" "--uncached" "$scale_name" "$subset_bcf" "$bench_dir/out_uncached"
      run_bench "$image" ""           "$scale_name" "$subset_bcf" "$bench_dir/out_cached"
    done
  done
}

main "$@"
