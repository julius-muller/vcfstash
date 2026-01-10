[![DOI](https://zenodo.org/badge/947952659.svg)](https://zenodo.org/badge/latestdoi/947952659)
[![CI](https://github.com/julius-muller/vcfcache/actions/workflows/ci.yml/badge.svg)](https://github.com/julius-muller/vcfcache/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/julius-muller/vcfcache)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/vcfcache)](https://pypi.org/project/vcfcache/)
[![Cite](https://img.shields.io/badge/Cite-CITATION.cff-blue)](CITATION.cff)
[![codecov](https://codecov.io/github/julius-muller/vcfcache/graph/badge.svg?token=ELV3PZ6PNL)](https://codecov.io/github/julius-muller/vcfcache)

# VCFcache – cache once, annotate fast

VCFcache builds a normalized blueprint of common variants, annotates it once, and reuses those annotations so only novel variants are processed at runtime.

## When VCFcache helps

VCFcache is useful when you either (a) repeatedly annotate many samples with a stable pipeline, or (b) want to quickly apply common annotations (e.g., VEP --everything) to a large VCF/BCF. Speed increase depends on cache hit rate of the input sample and per-variant annotation speed of the original pipeline.

## Key properties

* **Drop-in integration:** keep your existing annotation command; place it into a simple `annotation.yaml` and run `vcfcache annotate`.
* **Cache reuse with automatic fallback:** cache hits are reused; cache misses are annotated with your configured command and merged into one output.
* **Genome- and tool-agnostic:** works with arbitrary reference builds and organisms, and with any annotator or pipeline that can be expressed as a command (stdin/stdout or file-based).
* **Pre-built caches available:** published caches built from public aggregation resources and annotation tools can be downloaded and used immediately (currently hg19/hg38 annotated with VEP --everything, more on request); the same tooling can generate highly efficient custom caches for your specific options and datasets.
* **BCF-native I/O:** VCFcache reads and writes **BCF** for performance and indexing; use `bcftools view` to convert VCF/VCF.gz at the boundaries.

**Important**: to use a prebuilt cache, you must have the same annotation tool (and compatible version) installed locally.

See [WIKI.md](WIKI.md) for full documentation, performance notes, and cache distribution via Zenodo.

---

## Quick Start

### Installation

```bash
# Via pip (requires Python 3.11+ and bcftools >= 1.20)
uv pip install vcfcache

# Via Docker (includes bcftools)
docker pull ghcr.io/julius-muller/vcfcache:latest

# Via Apptainer
apptainer exec docker://ghcr.io/julius-muller/vcfcache:latest vcfcache --help
```

See [WIKI.md - Section 2 (Quick Start)](WIKI.md#2-quick-start) for development installation and troubleshooting.

---

## Your First Annotation (3-minute tutorial)

This minimal example shows the complete workflow using a public cache.

### 1. List available caches

```bash
vcfcache list caches
```

### 2. Pull a specific cache

```bash
# Cache auto-downloads on first use, or download explicitly:
vcfcache cache-build --doi 10.5281/zenodo.18189447
```

### 3. Check cache requirements

```bash
vcfcache annotate --requirements -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic
```

This shows:
- Required annotation tool version (e.g., VEP 115.2)
- Required params (reference cache paths, etc.)
- The exact annotation command that will run

**Critical**: Install the exact tool version shown (e.g., `vep --version` must match).

### 4. Annotate your sample

```bash
vcfcache annotate \
  -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic \
  --vcf sample.bcf \
  --output sample_annotated.bcf \
  --stats-dir ./results
```

**Input format**: VCFcache operates on BCF. Convert VCF/VCF.gz at boundaries:
```bash
bcftools view -Ob sample.vcf.gz | vcfcache annotate -a <cache-alias> -i - -o - | ...
```

See [WIKI.md - Section 7 (Using a cache to annotate samples)](WIKI.md#7-using-a-cache-to-annotate-samples) for all annotation options.

---

## Building Your Own Cache

If you need different annotation settings (plugins, flags, tool version):

### From a public cache's blueprint
```bash
# Downloaded caches include the blueprint! Use it to build a cache variant:
# 1. Download cache (includes blueprint)
vcfcache cache-build --doi 10.5281/zenodo.18189447

# 2. Create annotation.yaml with your pipeline (see below)

# 3. Build cache from the downloaded blueprint
vcfcache cache-build \
  --db ~/.cache/vcfcache/caches/<cache-alias> \
  --name my_vep_cache \
  -a annotation.yaml \
  -y params.yaml
```

### From a public blueprint only
```bash
# 1. Download blueprint (smaller, no annotations)
vcfcache blueprint-init --doi <blueprint_DOI> -o ./cache_root

# 2. Create annotation.yaml with your pipeline (see below)

# 3. Build cache
vcfcache cache-build \
  --db ./cache_root \
  --name my_vep_cache \
  -a annotation.yaml \
  -y params.yaml
```

### From your own variants
```bash
# Use your cohort's common variants for maximum cache hit rate
vcfcache blueprint-init --vcf cohort_common.bcf --output ./cache_root -y params.yaml
vcfcache cache-build --db ./cache_root --name my_cache -a annotation.yaml -y params.yaml
```

See [WIKI.md - Section 6 (Building your own cache)](WIKI.md#6-building-your-own-cache-end-to-end) for the complete workflow including sharing via Zenodo.

---

## Setting Up annotation.yaml

The annotation.yaml defines your annotation pipeline **at cache-build time**. It is **immutable**: once baked into a cache, it cannot be changed without rebuilding.

**What goes in annotation.yaml vs params.yaml:**
- **annotation.yaml** (immutable): annotation logic, tool flags, plugins — anything that affects the annotation semantics
- **params.yaml** (runtime): machine-specific paths, threads, tool locations — settings that differ between environments but don't change the annotation results

**When you annotate samples**, you can only provide params.yaml (via `-y`). The annotation.yaml is fixed in the cache.

**Minimal VEP example**:
```yaml
annotation_cmd: |
  ${params.bcftools_cmd} view -Ov ${INPUT_BCF} | \
  ${params.annotation_tool_cmd} \
    --offline --cache --vcf \
    --dir_cache ${params.vep_cache} \
    --assembly ${params.genome_build} \
    --stats_file ${AUXILIARY_DIR}/vep_stats.html \
    -i stdin -o stdout | \
  ${params.bcftools_cmd} view -Ob -o ${OUTPUT_BCF} -W

must_contain_info_tag: CSQ
required_tool_version: "115.2"
genome_build: "GRCh38"
```

**Variable substitution**:
- `${INPUT_BCF}`, `${OUTPUT_BCF}`, `${AUXILIARY_DIR}`: provided by VCFcache automatically
- `${params.*}`: substituted from your params.yaml at runtime (allows different paths per machine)
- Use `${params.*}` for paths/settings that vary between machines (e.g., `/path/to/vep/cache`)
- Hardcode values that define the annotation semantics (e.g., `--offline --cache`)

**params.yaml example**:
```yaml
genome_build: "GRCh38"
bcftools_cmd: "bcftools"
annotation_tool_cmd: "vep"
vep_cache: "/path/to/vep/cache"  # Machine-specific path
threads: 8                        # Machine-specific setting
```

See [WIKI.md - Section 8 (Configuration reference)](WIKI.md#8-configuration-reference-paramsyaml--annotationyaml) for full configuration details and advanced examples.

---

## Links

- **Full Documentation**: [WIKI.md](WIKI.md)
- **Performance Model**: [WIKI.md - Section 11](WIKI.md#11-performance-model-runtime-efficiency)
- **CLI Reference**: [WIKI.md - Section 10](WIKI.md#10-cli-reference-all-commands--flags)
- **Source**: https://github.com/julius-muller/vcfcache
- **Issues**: https://github.com/julius-muller/vcfcache/issues
- **Docker**: ghcr.io/julius-muller/vcfcache
