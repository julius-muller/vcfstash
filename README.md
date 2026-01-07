[![DOI](https://zenodo.org/badge/947952659.svg)](https://zenodo.org/badge/latestdoi/947952659)
[![CI](https://github.com/julius-muller/vcfcache/actions/workflows/ci.yml/badge.svg)](https://github.com/julius-muller/vcfcache/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/julius-muller/vcfcache)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/vcfcache)](https://pypi.org/project/vcfcache/)
[![Cite](https://img.shields.io/badge/Cite-CITATION.cff-blue)](CITATION.cff)
[![codecov](https://codecov.io/github/julius-muller/vcfcache/graph/badge.svg?token=ELV3PZ6PNL)](https://codecov.io/github/julius-muller/vcfcache)


# VCFcache – Cache once, annotate fast

Cache common variants once, reuse them for every sample. VCFcache builds a normalized blueprint, annotates it once, and reuses those results so only rare/novel variants are annotated at runtime.

**Performance**: With 60-90% cache hit rates on typical samples, VCFcache achieves 2-10× speed-ups compared to standard annotation pipelines. Cache lookups are constant-time operations regardless of cache size, making the tool highly scalable. See [WIKI.md](WIKI.md#performance-model) for the detailed runtime efficiency model.

**Contig naming**: VCFcache requires matching contig naming between cache and samples (e.g., both use `chr1` or both use `1`). At runtime, vcfcache reports the contig overlap and fails fast if there is no overlap.

Works with any genome/build (human, mouse, plants, model organisms) as long as your inputs and annotation pipeline use compatible reference builds and contig naming. The genome build must be explicitly set in both `params.yaml` and `annotation.yaml`.

---

## Quick Start - pip install

Requires: Python >= 3.11 (earlier versions untested), bcftools >= 1.20

```bash
pip install vcfcache
vcfcache demo --smoke-test  # Run comprehensive demo
vcfcache --help
```

Install bcftools separately:
- Ubuntu/Debian: `sudo apt-get install bcftools`
- macOS: `brew install bcftools`
- Conda: `conda install -c bioconda bcftools`

---

## Quick Start - Docker

**Docker includes bcftools** - no separate installation needed.

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

---

## Quick Start - from source

```bash
git clone https://github.com/julius-muller/vcfcache.git
cd vcfcache
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
vcfcache --help
```

---

## Build Your Own Cache

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
vcfcache annotate -a ./cache/cache/vep_cache --vcf sample.vcf.gz --output ./sample_vc.bcf --stats-dir ./results
```

If `--stats-dir` is omitted, stats are written to `<cwd>/<input_basename>_vcstats`.
Use `--no-stats` to skip writing stats/logs (disables `vcfcache compare`).

---

## Configuration

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
