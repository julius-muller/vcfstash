# VCFcache – Accelerate Variant Annotation by 70%+

**Cache common variants once, annotate samples instantly.**

VCFcache builds a local cache of pre-annotated variants and intelligently applies them to your samples, reducing annotation time by 70%+ for typical workflows. Works with any annotation tool (VEP, ANNOVAR, SnpEff, custom scripts).

---

## What is VCFcache?

VCFcache solves a simple problem: **you're annotating the same common variants over and over**. Instead:

1. **Build a cache** from representative variants (e.g., gnomAD)
2. **Annotate the cache once** with your preferred tool (VEP, etc.)
3. **Reuse cached annotations** for all your samples

Only novel/rare variants get annotated in real-time. Common variants (60-90% of most samples) come from cache instantly.

---

## Three Deployment Modes

### 🔧 **Mode 1: From Source (Development)**
Full flexibility, requires manual setup.

**Requirements**: Python 3.13+, uv (package manager)

```bash
git clone https://github.com/julius-muller/vcfcache.git
cd vcfcache
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
vcfcache --help
```

**Note**: bcftools 1.22 is bundled with the package, no separate installation needed.

**Use when**: Developing vcfcache, customizing workflows, or running tests.

**Tests**: `python -m pytest tests/ -v`

---

### 🐳 **Mode 2: Blueprint Docker (Build Your Own Cache)**
Lightweight image for creating custom caches.

**Requirements**: Docker, gnomAD BCF file, ~10GB disk space

```bash
# Pull blueprint image
docker pull ghcr.io/julius-muller/vcfcache-blueprint:latest

# Create cache from your gnomAD subset
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/cache:/cache \
  ghcr.io/julius-muller/vcfcache-blueprint:latest \
  blueprint-init \
    --vcf /data/gnomad_subset.bcf \
    --output /cache \
    -y /app/recipes/docker-cache/params.yaml

# Annotate cache with bcftools (mock annotation for testing)
docker run --rm \
  -v $(pwd)/cache:/cache \
  ghcr.io/julius-muller/vcfcache-blueprint:latest \
  cache-build \
    --name my_annotation \
    --db /cache \
    -a /app/recipes/docker-cache/annotation.yaml \
    -y /app/recipes/docker-cache/params.yaml
```

**Use when**: Building custom caches with specific variant sets or annotation tools.

**Tests**: Mock bcftools annotation (no VEP required)

---

### 🚀 **Mode 3: Annotated Docker (Production Ready)**
Pre-built cache with VEP annotations, ready to use.

**Requirements**: Docker, VEP cache (for annotation updates), ~5GB disk space

```bash
# Pull pre-built image with gnomAD cache + VEP annotations
docker pull ghcr.io/julius-muller/vcfcache-annotated:gnomad-grch38-joint-af010-vep115

# Annotate your sample (VEP cache hit rate: 60-90% typical)
docker run --rm \
  -v $(pwd)/samples:/data \
  -v $(pwd)/results:/output \
  ghcr.io/julius-muller/vcfcache-annotated:gnomad-grch38-joint-af010-vep115 \
  annotate \
    -a /cache/db/cache/vep_gnomad \
    --vcf /data/sample.vcf.gz \
    --output /output \
    -y /app/recipes/docker-annotated/params.yaml

# Output: /output/sample_vst.bcf with VEP annotations
```

**Use when**: Production annotation pipelines, reproducible research, instant setup.

**Cache variants**: gnomAD v4.1 GRCh38, AF ≥ 1% (~10 million variants)

**Tests**: Run with VEP cache mounted
```bash
docker run --rm \
  -v /path/to/vep/cache:/opt/vep/.vep:ro \
  --entrypoint /bin/sh \
  ghcr.io/julius-muller/vcfcache-annotated:latest \
  -c 'cd /app && python3.13 -m pytest tests/ -v'
```

---

## Quick Start Example

### Vanilla (from source)

```bash
# 1. Create cache from gnomAD
vcfcache blueprint-init \
  --vcf gnomad_chr1_af0.01.bcf \
  --output ./cache \
  -y params.yaml

# 2. Annotate cache with VEP
vcfcache cache-build \
  --name vep_gnomad \
  --db ./cache \
  -a annotation.yaml \
  -y params.yaml

# 3. Annotate your samples
vcfcache annotate \
  -a ./cache/cache/vep_gnomad \
  --vcf sample1.vcf.gz \
  --output ./results \
  -y params.yaml
```

### Docker (annotated mode)

```bash
# One command - instant annotation
docker run --rm \
  -v $(pwd):/work \
  ghcr.io/julius-muller/vcfcache-annotated:latest \
  annotate \
    -a /cache/db/cache/vep_gnomad \
    --vcf /work/sample.vcf.gz \
    --output /work/results
```

---

## Configuration Files (YAML)

VCFcache uses two YAML files for configuration:

### 1. `params.yaml` - Resource Paths

Defines where tools and resources are located:

```yaml
# Tool paths
bcftools_cmd: "bcftools"
annotation_tool_cmd: "vep"
tool_version_command: "vep --version"

# Resources
chr_add: "${VCFCACHE_ROOT}/resources/chr_add.txt"
temp_dir: "/tmp"

# VEP settings (example)
vep_buffer: 500000
vep_forks: 8
vep_cache: "/path/to/vep/cache"

# Optional validation
optional_checks:
  genome_build: "GRCh38"
  vep_cache_version: "115"
```

### 2. `annotation.yaml` - Annotation Command

Defines how to annotate variants:

```yaml
# Command to run (supports variable substitution)
annotation_cmd: |
  ${params.bcftools_cmd} view ${INPUT_BCF} | \
  ${params.annotation_tool_cmd} \
    --cache \
    --dir_cache ${params.vep_cache} \
    --format vcf \
    --vcf \
    -i STDIN \
    -o STDOUT | \
  ${params.bcftools_cmd} view -o ${OUTPUT_BCF} -Ob -W

# Validation
must_contain_info_tag: CSQ
required_tool_version: "115.2"

optional_checks:
  vep_cache_version: "115"
```

**Variables available**:
- `${params.*}` - Any key from params.yaml
- `${INPUT_BCF}` - Input file path
- `${OUTPUT_BCF}` - Output file path
- `${AUXILIARY_DIR}` - Directory for temporary files

---

## Performance

Typical speedups with gnomAD AF ≥ 1% cache:

| Sample Type | Cache Hit Rate | Speedup |
|-------------|----------------|---------|
| Population study (WGS) | 85-95% | 8-15x |
| Clinical exome | 60-80% | 3-5x |
| Rare disease | 40-60% | 2-3x |

**Why**: Most samples share common variants (dbSNP, gnomAD). VCFcache caches these once and reuses them across all samples.

---

## Documentation

- **[WIKI.md](WIKI.md)**: Comprehensive guide - cache structure, configuration, best practices, troubleshooting
- **[CLAUDE.md](CLAUDE.md)**: Developer guide for working with the codebase
- **[CHANGELOG.md](CHANGELOG.md)**: Version history and updates

---

## Requirements by Mode

| Mode | Python | Docker | VEP | Disk Space |
|------|--------|--------|-----|------------|
| **Vanilla** | 3.13+ | - | Optional | 5-10GB |
| **Blueprint** | - | ✓ | - | 2GB |
| **Annotated** | - | ✓ | Cache for updates | 5GB |

**Note**: bcftools 1.22 is bundled with all modes, no separate installation needed.

---

## Architecture

VCFcache uses a **pure Python workflow system** (no Java/Nextflow dependency):

- **Workflow Manager**: Orchestrates bcftools commands via subprocess
- **Database**: SQLite-backed blueprint and cache management
- **Annotation**: 4-step caching strategy (add cache → filter misses → annotate misses → merge back)

---

## Credits

- **[bcftools](https://github.com/samtools/bcftools)**: Efficient VCF/BCF manipulation
- **[uv](https://github.com/astral-sh/uv)**: Fast Python package management
- **[VEP](https://ensembl.org/vep)**: Variant annotation

---

## License

MIT License - see [LICENSE](LICENSE) file

---

## Citation

If you use VCFcache in your research, please cite:

```
VCFcache: Fast Variant Annotation Caching
https://github.com/julius-muller/vcfcache
```
