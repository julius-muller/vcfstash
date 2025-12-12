# VCFcache – Cache once, annotate fast

Cache common variants once, reuse them for every sample. VCFcache builds a normalized blueprint, annotates it once, and reuses those results so only rare/novel variants are annotated at runtime.

---

## Quick Start - pip install

Requires: Python >= 3.11, bcftools >= 1.20

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

# Use a public cache from Zenodo
docker run --rm -v $(pwd):/work ghcr.io/julius-muller/vcfcache:latest \
  annotate \
    -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic \
    --vcf /work/sample.vcf.gz \
    --output /work/out \
    --force

# List available public caches
docker run --rm ghcr.io/julius-muller/vcfcache:latest list caches
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
vcfcache annotate -a ./cache/cache/vep_cache --vcf sample.vcf.gz --output ./results
```

---

## Configuration

Override system bcftools (if needed):
```bash
export VCFCACHE_BCFTOOLS=/path/to/bcftools-1.22
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
