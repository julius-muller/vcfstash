# VCFcache Wiki

Comprehensive guide to VCFcache configuration, deployment, and best practices.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Cache Structure](#cache-structure)
3. [Configuration Files](#configuration-files)
4. [Deployment Modes](#deployment-modes)
5. [Building Caches](#building-caches)
6. [Annotation Tools](#annotation-tools)
7. [Performance Optimization](#performance-optimization)
8. [Testing & Validation](#testing--validation)
9. [Troubleshooting](#troubleshooting)
10. [Docker Best Practices](#docker-best-practices)

---

## Quick Start

### 30-Second Test (Pre-built Cache)

```bash
# Pull image and annotate
docker pull ghcr.io/julius-muller/vcfcache-annotated:latest
docker run --rm \
  -v $(pwd)/samples:/data \
  -v $(pwd)/results:/output \
  ghcr.io/julius-muller/vcfcache-annotated:latest \
  annotate \
    -a /cache/db/cache/vep_gnomad \
    --vcf /data/sample.vcf.gz \
    --output /output
```

### 15-Minute Local Test

```bash
# Clone and install
git clone https://github.com/julius-muller/vcfcache.git
cd vcfcache
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"

# Run tests
python -m pytest tests/ -v
```

### Production Setup (2-3 hours)

See [Building Caches](#building-caches) section below.

---

## Cache Structure

VCFcache organizes caches in a standardized directory structure:

```
cache_dir/
├── blueprint/                    # Normalized variant blueprint
│   ├── vcfcache.bcf             # All variants from input VCFs (normalized, deduplicated)
│   ├── vcfcache.bcf.csi         # Index
│   └── sources.info             # Tracking file (MD5 hashes of input VCFs)
├── cache/                        # Annotation caches
│   └── <annotation_name>/       # Named annotation (e.g., "vep_gnomad")
│       ├── vcfcache_annotated.bcf     # Annotated blueprint (THE CACHE)
│       ├── vcfcache_annotated.bcf.csi # Index
│       ├── annotation.yaml            # Annotation command used
│       ├── params.snapshot.yaml       # Params at annotation time
│       └── auxiliary/                 # Optional tool outputs (logs, etc.)
├── db/                           # Database metadata (SQLite)
│   ├── cache.db                 # Cache metadata, version tracking
│   └── workflow/                # Workflow files (copied during init)
└── work/                         # Temporary workflow files (can be deleted)
```

### Key Concepts

- **Blueprint**: The normalized, deduplicated variant set used as the basis for all annotations
- **Cache**: A named annotation of the blueprint (can have multiple cachees per blueprint)
- **Sources tracking**: MD5 hashes prevent re-adding duplicate input VCFs

---

## Configuration Files

VCFcache uses **two YAML files** for configuration (no Groovy, no Java required):

### 1. `params.yaml` - Resource Paths and Settings

**Purpose**: Define where tools and resources are located. This file can be shared across multiple annotations.

**Structure**:

```yaml
## Required Resources (keys are fixed - don't rename)
bcftools_cmd: "bcftools"                    # Path to bcftools binary
annotation_tool_cmd: "vep"                  # Path to annotation tool
tool_version_command: "vep --version"       # Command to get tool version
chr_add: "${VCFCACHE_ROOT}/resources/chr_add.txt"  # Chromosome renaming file
temp_dir: "/tmp"                            # Temporary directory

## Optional Resources (can be referenced in annotation.yaml as ${params.KEY})
vep_buffer: 500000                          # VEP buffer size
vep_forks: 8                                # VEP parallel processes
vep_cache: "/path/to/vep/cache"            # VEP cache directory

## Optional Checks (validated during annotation)
optional_checks:
  genome_build: "GRCh38"                   # Expected genome build
  vep_cache_version: "115"                 # Expected VEP cache version
```

**Variable Substitution**:
- `${VCFCACHE_ROOT}`: Package installation directory
- Environment variables: `${HOME}`, `${USER}`, etc.

### 2. `annotation.yaml` - Annotation Command

**Purpose**: Define the exact command used to annotate variants. This file is **locked** after `cache-build` and stored with the cache to ensure reproducibility.

**Structure**:

```yaml
# Annotation command (uses pipes for streaming)
annotation_cmd: |
  ${params.bcftools_cmd} view ${INPUT_BCF} | \
  ${params.annotation_tool_cmd} \
    --offline \
    --cache \
    --dir_cache ${params.vep_cache} \
    --buffer_size ${params.vep_buffer} \
    --fork ${params.vep_forks} \
    --format vcf \
    --vcf \
    -i STDIN \
    -o STDOUT | \
  ${params.bcftools_cmd} view -o ${OUTPUT_BCF} -Ob -W

# Validation tag (must appear in output)
must_contain_info_tag: CSQ

# Tool version check
required_tool_version: "115.2"

# Additional checks (locked with cache)
optional_checks:
  vep_cache_version: "115"
  genome_build: "GRCh38"
```

**Special Variables** (automatically replaced during execution):
- `${INPUT_BCF}`: Path to input BCF file
- `${OUTPUT_BCF}`: Path to output BCF file
- `${AUXILIARY_DIR}`: Directory for temporary/auxiliary files
- `${params.*}`: Any key from params.yaml (e.g., `${params.vep_buffer}`)

---

## Deployment Modes

### Mode 1: Vanilla (From Source)

**Best for**: Development, testing, custom workflows

**Setup**:
```bash
git clone https://github.com/julius-muller/vcfcache.git
cd vcfcache
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
```

**Requirements**:
- Python 3.13+
- bcftools 1.20+
- tabix
- Annotation tool (VEP, etc.) - optional

**Tests**:
```bash
python -m pytest tests/ -v
```

---

### Mode 2: Blueprint Docker

**Best for**: Building custom caches, testing workflows

**Pull Image**:
```bash
docker pull ghcr.io/julius-muller/vcfcache-blueprint:latest
```

**Create Cache**:
```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/cache:/cache \
  ghcr.io/julius-muller/vcfcache-blueprint:latest \
  blueprint-init \
    --vcf /data/gnomad_subset.bcf \
    --output /cache \
    -y /app/recipes/docker-cache/params.yaml
```

**Annotate Cache** (mock bcftools annotation):
```bash
docker run --rm \
  -v $(pwd)/cache:/cache \
  ghcr.io/julius-muller/vcfcache-blueprint:latest \
  cache-build \
    --name my_annotation \
    --db /cache \
    -a /app/recipes/docker-cache/annotation.yaml \
    -y /app/recipes/docker-cache/params.yaml
```

**Tests** (no VEP required):
```bash
docker run --rm \
  --entrypoint /bin/sh \
  ghcr.io/julius-muller/vcfcache-blueprint:latest \
  -c 'cd /app && python3.13 -m pytest tests/ -v'
```

---

### Mode 3: Annotated Docker (Production)

**Best for**: Production annotation, reproducible research

**Pull Image**:
```bash
docker pull ghcr.io/julius-muller/vcfcache-annotated:gnomad-grch38-joint-af010-vep115
```

**Annotate Samples**:
```bash
docker run --rm \
  -v $(pwd)/samples:/data \
  -v $(pwd)/results:/output \
  ghcr.io/julius-muller/vcfcache-annotated:gnomad-grch38-joint-af010-vep115 \
  annotate \
    -a /cache/db/cache/vep_gnomad \
    --vcf /data/sample.vcf.gz \
    --output /output \
    -y /app/recipes/docker-annotated/params.yaml
```

**Cache Details**:
- **gnomAD v4.1** GRCh38
- **Joint exomes + genomes**
- **AF ≥ 1%** (~10 million variants)
- **VEP 115** with HGVS, SIFT, PolyPhen, transcript info

**Tests** (requires VEP cache mounted):
```bash
docker run --rm \
  -v /path/to/vep/cache:/opt/vep/.vep:ro \
  --entrypoint /bin/sh \
  ghcr.io/julius-muller/vcfcache-annotated:latest \
  -c 'cd /app && python3.13 -m pytest tests/ -v'
```

---

## Building Caches

### Preparing gnomAD Data

```bash
# Download gnomAD v4.1 (example: chr1 only)
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr1.vcf.bgz

# Filter by allele frequency (optional but recommended)
bcftools view -i 'AF>=0.01' \
  gnomad.genomes.v4.1.sites.chr1.vcf.bgz \
  -Ob -o gnomad_chr1_af0.01.bcf -W

# The index (.csi) is created automatically with -W flag
```

### Creating a Cache (Vanilla)

```bash
# 1. Initialize blueprint
vcfcache blueprint-init \
  --vcf gnomad_chr1_af0.01.bcf \
  --output ./my_cache \
  --normalize \
  -y params.yaml

# 2. Annotate blueprint
vcfcache cache-build \
  --name vep_gnomad \
  --db ./my_cache \
  -a annotation.yaml \
  -y params.yaml

# 3. Verify cache
ls -lh ./my_cache/cache/vep_gnomad/
bcftools stats ./my_cache/cache/vep_gnomad/vcfcache_annotated.bcf | grep "number of records"
```

### Creating a Cache (Docker)

```bash
# Build blueprint Docker image with your gnomAD data
./scripts/local-build/03-build-blueprint.sh gnomad_af0.01.bcf --push

# Build annotated Docker image with VEP
./scripts/local-build/04a-build-base-image.sh gnomad_af0.01.bcf --yes
./scripts/local-build/04b-annotate-and-commit.sh \
  --base-image <image-from-04a> \
  --vep-cache-dir /path/to/vep/cache \
  --yes
```

---

## Annotation Tools

### VEP (Ensembl Variant Effect Predictor)

**params.yaml**:
```yaml
annotation_tool_cmd: "vep"
tool_version_command: "vep --version | head -1"
vep_buffer: 500000
vep_forks: 8
vep_cache: "/path/to/vep/cache"
```

**annotation.yaml**:
```yaml
annotation_cmd: |
  ${params.bcftools_cmd} view ${INPUT_BCF} | \
  ${params.annotation_tool_cmd} \
    --offline \
    --cache \
    --dir_cache ${params.vep_cache} \
    --buffer_size ${params.vep_buffer} \
    --fork ${params.vep_forks} \
    --format vcf \
    --vcf \
    -i STDIN \
    -o STDOUT \
    --hgvs \
    --symbol \
    --canonical \
    --sift b \
    --polyphen b | \
  ${params.bcftools_cmd} view -o ${OUTPUT_BCF} -Ob -W

must_contain_info_tag: CSQ
required_tool_version: "115"
```

### SnpEff

**params.yaml**:
```yaml
annotation_tool_cmd: "snpEff"
tool_version_command: "snpEff -version 2>&1 | head -1"
snpeff_db: "GRCh38.99"
```

**annotation.yaml**:
```yaml
annotation_cmd: |
  ${params.bcftools_cmd} view ${INPUT_BCF} | \
  ${params.annotation_tool_cmd} \
    -noStats \
    ${params.snpeff_db} \
    /dev/stdin | \
  ${params.bcftools_cmd} view -o ${OUTPUT_BCF} -Ob -W

must_contain_info_tag: ANN
```

### ANNOVAR

**params.yaml**:
```yaml
annotation_tool_cmd: "table_annovar.pl"
annovar_db: "/path/to/annovar/humandb"
annovar_protocol: "refGene,gnomad41_genome"
annovar_operation: "g,f"
```

**annotation.yaml**:
```yaml
annotation_cmd: |
  ${params.bcftools_cmd} view ${INPUT_BCF} > ${AUXILIARY_DIR}/input.vcf && \
  ${params.annotation_tool_cmd} \
    ${AUXILIARY_DIR}/input.vcf \
    ${params.annovar_db} \
    -buildver hg38 \
    -out ${AUXILIARY_DIR}/output \
    -remove \
    -protocol ${params.annovar_protocol} \
    -operation ${params.annovar_operation} \
    -vcfinput && \
  ${params.bcftools_cmd} view -o ${OUTPUT_BCF} -Ob -W \
    ${AUXILIARY_DIR}/output.hg38_multianno.vcf

must_contain_info_tag: Func.refGene
```

---

## Performance Optimization

### Cache Hit Rates

Typical hit rates with gnomAD AF ≥ 1% cache:

| Sample Source | Cache Hits | Speedup |
|---------------|------------|---------|
| PopGen (WGS) | 85-95% | 8-15x |
| Clinical Exome | 60-80% | 3-5x |
| Rare Disease | 40-60% | 2-3x |
| De novo calls | 10-30% | 1.2-1.5x |

### Optimizing Cache Size vs. Hit Rate

AF threshold trade-offs:

| gnomAD AF | Cache Size | Variants | Typical Hit Rate |
|-----------|------------|----------|------------------|
| ≥ 10% | ~500MB | ~1M | 50-70% |
| ≥ 5% | ~1.5GB | ~3M | 65-80% |
| ≥ 1% | ~5GB | ~10M | 75-90% |
| ≥ 0.1% | ~20GB | ~50M | 85-95% |

**Recommendation**: Start with AF ≥ 1% (good balance of size vs. performance).

### Parallelization

**VEP forks**:
```yaml
vep_forks: 8  # Use number of CPU cores available
```

**bcftools threads** (not yet supported in all commands):
```yaml
bcftools_threads: 4
```

### Disk I/O Optimization

- Store cache on **fast SSD** for best performance
- Use **local cache** rather than network-mounted storage
- Enable **bcftools index caching** with `-W` flag (writes index alongside BCF)

---

## Testing & Validation

### Unit Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Run specific test file
python -m pytest tests/test_annotate.py -xvs

# Run with coverage
python -m pytest tests/ --cov=vcfcache --cov-report=html
```

### Docker Tests

**Blueprint image**:
```bash
docker run --rm \
  --entrypoint /bin/sh \
  ghcr.io/julius-muller/vcfcache-blueprint:latest \
  -c 'cd /app && python3.13 -m pytest tests/ -v'
```

**Annotated image** (requires VEP cache):
```bash
docker run --rm \
  -v /path/to/vep/cache:/opt/vep/.vep:ro \
  --entrypoint /bin/sh \
  ghcr.io/julius-muller/vcfcache-annotated:latest \
  -c 'cd /app && python3.13 -m pytest tests/ -v'
```

### Validating Annotations

```bash
# Check that annotations were applied
bcftools view -h annotated_sample_vst.bcf | grep "##INFO=<ID=CSQ"

# Count annotated variants
bcftools query -f '%INFO/CSQ\n' annotated_sample_vst.bcf | grep -v "^\.$" | wc -l

# Compare with direct annotation (should be identical)
diff <(bcftools view -H direct_annotated.vcf | cut -f1-8) \
     <(bcftools view -H cached_annotated_vst.vcf | cut -f1-8)
```

---

## Troubleshooting

### Common Issues

#### 1. "Failed to check bcftools version"

**Cause**: bcftools not in PATH or GLIBC incompatibility

**Solution**:
```bash
# Check bcftools version
bcftools --version

# For Docker images, ensure you're using compatible base
# Annotated images use bcftools 1.22 compiled with GLIBC 2.35 compat
```

#### 2. "Cache directory /opt/vep/.vep/homo_sapiens not found"

**Cause**: VEP cache not mounted or not at expected location

**Solution**:
```bash
# Mount VEP cache when running annotated Docker images
docker run --rm \
  -v /path/to/vep/cache:/opt/vep/.vep:ro \
  ...
```

#### 3. "must_contain_info_tag 'CSQ' not found in output"

**Cause**: Annotation command failed or didn't add expected tag

**Solution**:
- Check annotation command in `annotation.yaml`
- Verify annotation tool is working: `vep --help`
- Check `work/` directory logs for error messages
- Test annotation command manually

#### 4. "MD5 mismatch: input file has changed"

**Cause**: Input VCF was modified after being added to cache

**Solution**:
- Use `--force` flag to override (not recommended)
- Regenerate cache from scratch if source data changed

#### 5. "No variants in cache overlap with sample"

**Cause**: Different chromosome naming (chr1 vs. 1) or different genome build

**Solution**:
- Use `--normalize` flag during `blueprint-init` to standardize chromosome names
- Verify genome build matches between cache and sample

---

## Docker Best Practices

### Volume Mounts

```bash
# Best practice: use absolute paths
docker run --rm \
  -v /absolute/path/to/data:/data \
  -v /absolute/path/to/results:/output \
  vcfcache-image ...

# Avoid relative paths (can be ambiguous)
```

### Resource Limits

```bash
# Limit memory and CPUs
docker run --rm \
  --memory=16g \
  --cpus=8 \
  vcfcache-image ...
```

### Image Tags

```bash
# Use specific tags for reproducibility
docker pull ghcr.io/julius-muller/vcfcache-annotated:gnomad-grch38-joint-af010-vep115

# Avoid :latest in production
```

### Cleanup

```bash
# Remove stopped containers
docker container prune

# Remove unused images
docker image prune

# Remove all unused data
docker system prune -a
```

---

## Advanced Topics

### Multi-Sample Annotation

```bash
# Annotate multiple samples in parallel
for sample in samples/*.vcf.gz; do
  docker run --rm \
    -v $(pwd):/work \
    vcfcache-annotated:latest \
    annotate \
      -a /cache/db/cache/vep_gnomad \
      --vcf /work/$sample \
      --output /work/results &
done
wait
```

### Custom Chromosome Sets

**chr_add.txt** (for renaming chromosomes):
```
1 chr1
2 chr2
...
X chrX
Y chrY
MT chrM
```

### Cache Versioning

Store `optional_checks` in params.yaml to track cache versions:

```yaml
optional_checks:
  genome_build: "GRCh38"
  vep_cache_version: "115"
  gnomad_version: "4.1"
  cache_created: "2025-12-04"
```

These are validated during annotation to prevent version mismatches.

---

## Getting Help

- **GitHub Issues**: https://github.com/julius-muller/vcfcache/issues
- **Documentation**: See [README.md](README.md) and [CLAUDE.md](CLAUDE.md)
- **Tests**: Run `python -m pytest tests/ -v` to see example workflows

---

**Last Updated**: 2025-12-04
