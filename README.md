[![DOI](https://zenodo.org/badge/947952659.svg)](https://zenodo.org/badge/latestdoi/947952659)
[![CI](https://github.com/julius-muller/vcfcache/actions/workflows/ci.yml/badge.svg)](https://github.com/julius-muller/vcfcache/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/julius-muller/vcfcache)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/vcfcache)](https://pypi.org/project/vcfcache/)
[![Cite](https://img.shields.io/badge/Cite-CITATION.cff-blue)](CITATION.cff)
[![codecov](https://codecov.io/github/julius-muller/vcfcache/graph/badge.svg?token=ELV3PZ6PNL)](https://codecov.io/github/julius-muller/vcfcache)


# VCFcache – cache once, annotate fast

VCFcache builds a normalized blueprint of common variants, annotates it once, and reuses those annotations so only novel variants are processed at runtime. It is genome‑agnostic and tool‑agnostic (VEP, SnpEff, ANNOVAR, custom scripts).

Important: to use a cache, you must have the same annotation tool (and compatible version) installed locally.

See [WIKI.md](WIKI.md) for full documentation, performance notes, and cache distribution via Zenodo.

---

## Quick Start (containers: Docker or Apptainer)

Containers include a modern `bcftools`, which avoids OS‑level version issues.

```bash
docker pull ghcr.io/julius-muller/vcfcache:latest

# List available public caches
docker run --rm ghcr.io/julius-muller/vcfcache:latest list caches

# Use a public cache from Zenodo
docker run --rm -v $(pwd):/work ghcr.io/julius-muller/vcfcache:latest \
  annotate \
    -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic \
    --vcf /work/sample.vcf.gz \
    --output /work/sample_vc.bcf \
    --stats-dir /work
```

Apptainer/Singularity (example):
```bash
apptainer exec docker://ghcr.io/julius-muller/vcfcache:latest \
  vcfcache list caches
```

---

## Quick Start (pip)

Requires: Python >= 3.11 and `bcftools >= 1.20`.

```bash
uv pip install vcfcache
vcfcache demo -q
vcfcache --help
```

Install `bcftools` separately:
- Ubuntu/Debian: `sudo apt-get install bcftools`
- macOS: `brew install bcftools`
- Conda: `conda install -c bioconda bcftools`

---

## Quick Start (from source with test suite)

```bash
git clone https://github.com/julius-muller/vcfcache.git
cd vcfcache
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
vcfcache --help
```

---

## Build your own cache

1. **Create blueprint** (normalize/deduplicate variants):
```bash
vcfcache blueprint-init --vcf gnomad.bcf --output ./cache -y params.yaml
```

2. **Annotate blueprint** (create cache):
```bash
vcfcache cache-build --name vep_cache --db ./cache -a annotation.yaml -y params.yaml
```

3. **Use cache** on samples:
```bash
vcfcache annotate --requirements -a ./cache/cache/vep_cache
vcfcache annotate -a ./cache/cache/vep_cache --vcf sample.vcf.gz --output ./sample_vc.bcf --stats-dir ./results
```

If `--stats-dir` is provided, stats are written to `<stats_dir>/<input_basename>_vcstats`. Otherwise, they default to `<cwd>/<input_basename>_vcstats`.
Use `--no-stats` to skip writing stats/logs (disables `vcfcache compare`).

---

## Public caches/blueprints (Zenodo)

List and download:
```bash
vcfcache list caches
vcfcache list blueprints
vcfcache cache-build --doi <DOI>                 # download cache
vcfcache blueprint-init --doi <DOI> -o ./cache   # download blueprint
```
Use a downloaded cache:
```bash
vcfcache annotate -a ~/.cache/vcfcache/caches/<cache_name> --vcf sample.vcf.gz --output sample_vc.bcf
```

---

## params.yaml (runtime settings)

`params.yaml` defines tool paths and runtime settings (e.g., `bcftools_cmd`, `annotation_tool_cmd`, `threads`, `temp_dir`, `genome_build`).  
Pass it with `-y/--yaml` for `cache-build` and `annotate`. If omitted for `annotate`, the cache’s `params.snapshot.yaml` is used.

## annotation.yaml (annotation recipe)

`annotation.yaml` defines the annotation command, required tool version, and output tag (`must_contain_info_tag`).  
It is required for `cache-build` (`-a/--anno-config`) and is stored in the cache as `annotation.snapshot.yaml`.

The key field is `annotation_cmd`. It is a shell command string that must read from `$INPUT_BCF` and write to `$OUTPUT_BCF`.  
You can include `$AUXILIARY_DIR` for tool side‑outputs. This is typically a direct translation of your annotation pipeline.

To see requirements for a downloaded cache, run:
```bash
vcfcache annotate --requirements -a <cache_dir>
```

For publishing to Zenodo, prepare a metadata YAML and use `vcfcache push` (details in the wiki).

Minimal example (bcftools annotate):
```yaml
annotation_cmd: "bcftools annotate -a /path/to/anno.bcf -c INFO -o $OUTPUT_BCF -Ob -W --threads ${params.threads} $INPUT_BCF"
must_contain_info_tag: "CSQ"
required_tool_version: "1.0"
genome_build: "GRCh38"
```

---

## Stats directory contents

`<input_basename>_vcstats` contains:
- `annotation.log` and `workflow.log` — run logs and timing breakdown.
- `workflow/params.snapshot.yaml` and `workflow/annotation.snapshot.yaml` — exact configs used for the run.
- `auxiliary/` — tool side‑outputs (e.g., VEP stats).
- `compare_stats.yaml` and `.vcfcache_complete` — summary metrics and completion metadata.

---

## Configuration (common)

Override system bcftools (if needed):
```bash
export VCFCACHE_BCFTOOLS=/path/to/bcftools-1.22
```

Change where downloaded caches/blueprints are stored (default: `~/.cache/vcfcache`):
```bash
export VCFCACHE_DIR=/path/to/vcfcache_cache_dir
```

Or in `params.yaml`:
```yaml
bcftools_cmd: "/path/to/bcftools"
```

See [WIKI.md](WIKI.md) for detailed configuration, cache distribution via Zenodo, and troubleshooting.

---

## Links

- **Documentation**: [WIKI.md](WIKI.md)
- **Source**: https://github.com/julius-muller/vcfcache
- **Issues**: https://github.com/julius-muller/vcfcache/issues
- **Docker**: ghcr.io/julius-muller/vcfcache
