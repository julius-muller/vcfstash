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
#       ghcr.io/julius-muller/vcfcache-annotated:gnomad-grch38-joint-af010-vep115
#       ghcr.io/julius-muller/vcfcache-annotated:gnomad-grch38-joint-af001-vep115
#   - VEP cache directory (set VEP_CACHE_DIR)
#
# Usage:
#   ./tests/run_benchmarks_subset.sh [OPTIONS] [SOURCE_BCF]
#   Options:
#     -a    Append mode: always run benchmarks even if they exist in logs
#     -f    Force mode: delete all existing logs and re-run everything
#
#   SOURCE_BCF defaults to /mnt/data/samples/test_mgm/mgm_WGS_32.gatkWGS_norm_hg38.bcf
#
# Environment Variables:
#   VEP_CACHE_DIR   VEP cache directory (default: /mnt/data/apps/ensembl-vep/115/cachedir)
#   TEMP_DIR        Temporary directory for annotation files (default: /mnt/data/tmp)
#
# Output:
#   Logs are written to ./tests/benchmarks/ in TSV format:
#     timestamp\timage\tmode\tscale\tvariants\tseconds\tstatus\toutput_bcf

# Parse command-line flags
APPEND_MODE=false
FORCE_MODE=false
while getopts "af" opt; do
  case $opt in
    a) APPEND_MODE=true ;;
    f) FORCE_MODE=true ;;
    *) echo "Usage: $0 [-a] [-f] [SOURCE_BCF]" >&2; exit 1 ;;
  esac
done
shift $((OPTIND-1))

SOURCE_BCF=${1:-/mnt/data/samples/test_mgm/mgm_WGS_32.gatkWGS_norm_hg38.bcf}
VEP_CACHE_DIR=${VEP_CACHE_DIR:-/mnt/data/apps/ensembl-vep/115/cachedir}
TEMP_DIR=${TEMP_DIR:-/mnt/data/tmp}

# Create temp directory if it doesn't exist
mkdir -p "$TEMP_DIR"

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
  "ghcr.io/julius-muller/vcfcache-annotated:gnomad-grch38-joint-af010-vep115"
  "ghcr.io/julius-muller/vcfcache-annotated:gnomad-grch38-joint-af001-vep115"
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

  # Format mode for display (remove leading dashes)
  local mode_display="${mode#--}"
  mode_display="${mode_display:-cached}"

  # Print what we're about to benchmark
  echo ">>> Running: ${image##*:} | $mode_display | $scale"

  mkdir -p "$outdir"
  local bname
  bname="$(basename "$bcf")"
  local out_name="${bname%.bcf}_vc.bcf"
  local run_name="run_${mode:-cached}_${scale}"
  local run_dir_host="${outdir}/${run_name}"
  local run_dir_cont="/out/${run_name}"
  local outfile="${run_dir_host}/${out_name}"

  # Check log file status
  local total_lines=0
  if [ -f "$LOG_FILE" ]; then
    total_lines=$(wc -l < "$LOG_FILE")
  fi

  # If log has less than 2 lines (header only or corrupted), delete it
  if [ -f "$LOG_FILE" ] && [ "$total_lines" -lt 2 ]; then
    echo "    Log incomplete, deleting and re-running"
    rm -f "$LOG_FILE"
    total_lines=0
  fi

  # Skip logic (only applies if NOT in append mode and NOT in force mode)
  if [ "$APPEND_MODE" = false ] && [ "$total_lines" -ge 2 ]; then
    # Create a precise grep pattern that matches tab-separated fields
    # Log format: timestamp\timage\tmode\tscale\tvariants\tseconds\tstatus\toutput
    # We need to match: \t{mode}\t{scale}\t
    local grep_pattern=$'\t'"${mode}"$'\t'"${scale}"$'\t'
    # Check if this specific benchmark already exists in the log
    if grep -qF "${grep_pattern}" "$LOG_FILE" 2>/dev/null; then
      echo "    Skipped - already completed (use -a to append)"
      return 0
    fi
  fi

  # If we get here, we're running the benchmark
  # ensure fresh run dir on host (vcfcache will create it)
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
    -v "$TEMP_DIR":/tmp \
    -v "${project_root}/vcfcache":/app/venv/lib/python3.13/site-packages/vcfcache:ro \
    -w /app \
    --user "$(id -u):$(id -g)" \
    --entrypoint /bin/bash \
    "$image" \
    -lc "NXF_HOME=${run_dir_cont}/.nxf NXF_WORK=${run_dir_cont}/work \
         VCFCACHE_LOGLEVEL=ERROR VCFCACHE_FILE_LOGLEVEL=ERROR NXF_ANSI_LOG=false \
         vcfcache annotate ${mode} \
         --force \
         -a /cache/db/cache/vep_gnomad \
         --vcf /work/input.bcf \
         --output ${run_dir_cont}/${out_name} \
         --stats-dir ${run_dir_cont} \
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
    find "${run_dir_host}" -name "*.bcf" -o -name "*.bcf.csi" -o -name "*.html" -print0 2>/dev/null | xargs -0 rm -f 2>/dev/null || true
  fi
}

main() {
  # Pre-flight summary
  echo "=== VCFcache Benchmark Plan ==="
  echo "Source: $SOURCE_BCF"
  echo "Temp directory: $TEMP_DIR"
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
  echo "Docker images (cached mode):"
  for image in "${IMAGES[@]}"; do
    echo "  - ${image##*:}"
  done
  echo ""
  echo "Modes:"
  echo "  - Uncached: ${#SCALES[@]} benchmarks (once per scale)"
  echo "  - Cached: $((${#SCALES[@]} * ${#IMAGES[@]})) benchmarks (all scales × all cache AFs)"
  echo "Total benchmarks: $((${#SCALES[@]} + ${#SCALES[@]} * ${#IMAGES[@]}))"
  echo ""
  if [ "$FORCE_MODE" = true ]; then
    echo "FORCE MODE: Will delete all existing logs and re-run"
  elif [ "$APPEND_MODE" = true ]; then
    echo "APPEND MODE: Will run all benchmarks and append to logs"
  else
    echo "DEFAULT MODE: Will skip already completed benchmarks"
  fi
  echo "==============================="
  echo ""

  # Force mode: delete all existing logs
  if [ "$FORCE_MODE" = true ]; then
    echo "Deleting all existing log files..."
    find "$LOG_DIR" -name "*.log" -type f -delete 2>/dev/null || true
  fi

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

    # Run UNCACHED once per scale (cache AF doesn't matter for uncached)
    # Use first image for consistency
    first_image="${IMAGES[0]}"
    LOG_FILE="$bench_dir/uncached_${scale_name}.log"
    if [ ! -f "$LOG_FILE" ]; then
      tsv_log "timestamp\timage\tmode\tscale\tvariants\tseconds\tstatus\toutput_bcf"
    fi
    run_bench "$first_image" "--uncached" "$scale_name" "$subset_bcf" "$bench_dir/out_uncached"

    # Run CACHED for each cache image (AF010, AF001, etc.)
    for image in "${IMAGES[@]}"; do
      LOG_FILE="$bench_dir/${image##*:}_${scale_name}.log"
      if [ ! -f "$LOG_FILE" ]; then
        tsv_log "timestamp\timage\tmode\tscale\tvariants\tseconds\tstatus\toutput_bcf"
      fi
      run_bench "$image" "" "$scale_name" "$subset_bcf" "$bench_dir/out_cached"
    done
  done
}

main "$@"
