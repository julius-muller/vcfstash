# Local Docker Image Build Pipeline

This directory contains scripts for building VCFstash Docker images locally. This is the **recommended approach** for production builds, as it provides better control, faster iteration, and doesn't consume GitHub Actions minutes.

## Why Build Locally?

**Problems with GitHub Actions for large genomics data:**
- ⏱️ **Time limits:** 6-hour max per job (full genome build takes 3+ hours)
- 💰 **Cost:** Limited free minutes (2000/month = ~11 full builds)
- 🔧 **Resources:** Limited CPU/RAM for large-scale data processing
- 🐛 **Debugging:** Difficult to debug 3-hour builds
- 📦 **Storage:** Artifact retention limits and size constraints

**Benefits of local builds:**
- ✅ **Full control** over hardware and resources
- ✅ **No time/quota limits**
- ✅ **Easy debugging** and iteration
- ✅ **Faster** with better hardware
- ✅ **Master BCF strategy** - generate once, reuse forever
- ✅ **Reproducible** with version-controlled scripts

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────┐
│  1. Generate Master BCF (AF ≥ 0.001 = 0.1%)                    │
│     • One-time Hail export from gnomAD                          │
│     • Takes 1-3 hours for full genome                           │
│     • Stores all variants with AF ≥ 0.001                       │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  2. Subset to AF Thresholds                                     │
│     • Fast bcftools filtering (minutes, not hours!)             │
│     • Create AF ≥ 0.01, 0.05, 0.10 versions                     │
│     • Reuse master BCF for different thresholds                 │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  3. Build Blueprint Images (lean)                               │
│     • Data-only Docker images                                   │
│     • No VEP, no annotation tools                               │
│     • Fast builds (~5-10 minutes)                               │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  4. Build Annotated Images (fat) - Optional                     │
│     • Full VEP pipeline + pre-annotated cache                   │
│     • Slower builds (~30-60 minutes)                            │
│     • Ready-to-use for production                               │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  5. Push to GitHub Container Registry                           │
│     • Upload finished images to GHCR                            │
│     • Share with team / deploy to production                    │
└─────────────────────────────────────────────────────────────────┘
```

## Quick Start

### Prerequisites

```bash
# Install system dependencies
sudo apt-get install bcftools tabix samtools

# Install Python dependencies
pip install hail

# Install Docker
# See: https://docs.docker.com/engine/install/

# Login to GitHub Container Registry
echo $GITHUB_TOKEN | docker login ghcr.io -u USERNAME --password-stdin
```

### Full Pipeline

```bash
# Step 1: Generate master BCF (one-time, 1-3 hours)
./01-generate-master-bcf.sh \
  --genome GRCh38 \
  --type joint \
  --chromosomes all \
  --master-af 0.001 \
  --output-dir ./data

# Step 2: Subset to different AF thresholds (fast!)
./02-subset-af-thresholds.sh ./data/gnomad_v4.1_GRCh38_joint_master_af0001.bcf \
  --af-thresholds "0.01 0.05 0.10" \
  --output-dir ./data

# Step 3: Build blueprint images
./03-build-blueprint.sh ./data/gnomad_v4.1_GRCh38_joint_af010.bcf \
  --genome GRCh38 \
  --type joint \
  --push

# Step 4: Build annotated image (optional)
./04-build-annotated.sh ./data/gnomad_v4.1_GRCh38_joint_af010.bcf \
  --genome GRCh38 \
  --type joint \
  --vep-version 115.2 \
  --push

# Or push all at once
./05-push-images.sh --type both
```

### Test First with chrY

```bash
# Step 1: Generate test BCF (chrY only, ~5 minutes)
./01-generate-master-bcf.sh \
  --chromosomes chrY \
  --master-af 0.001 \
  --output-dir ./data

# Step 2: Subset
./02-subset-af-thresholds.sh ./data/gnomad_v4.1_GRCh38_joint_chrY_master_af0001.bcf \
  --af-thresholds "0.10" \
  --output-dir ./data

# Step 3: Build and test
./03-build-blueprint.sh ./data/gnomad_v4.1_GRCh38_joint_chrY_af010.bcf

# Verify it works!
docker run --rm ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af010 -v
```

## Script Reference

### 01-generate-master-bcf.sh

Generates the master BCF file with minimal AF filtering (0.1%).

**Usage:**
```bash
./01-generate-master-bcf.sh [OPTIONS]
```

**Options:**
- `--genome GENOME` - Genome build (default: GRCh38)
- `--gnomad-version VER` - gnomAD version (default: 4.1)
- `--type TYPE` - Dataset type: joint/genomes/exomes (default: joint)
- `--chromosomes CHR` - Chromosomes: all/chrY/chr22 (default: all)
- `--output-dir DIR` - Output directory (default: ./data)
- `--master-af AF` - Master AF threshold (default: 0.001)

**Output:**
- `gnomad_v4.1_GRCh38_joint_master_af0001.bcf` (and .csi index)

**Time:** 1-3 hours for full genome, 5-10 minutes for chrY

### 02-subset-af-thresholds.sh

Subsets the master BCF to different AF thresholds.

**Usage:**
```bash
./02-subset-af-thresholds.sh MASTER_BCF [OPTIONS]
```

**Options:**
- `--af-thresholds "X Y Z"` - Space-separated AF values (default: "0.01 0.05 0.10")
- `--output-dir DIR` - Output directory (default: ./data)
- `--threads N` - Number of threads (default: 4)

**Output:**
- `gnomad_v4.1_GRCh38_joint_af010.bcf`
- `gnomad_v4.1_GRCh38_joint_af005.bcf`
- `gnomad_v4.1_GRCh38_joint_af001.bcf`

**Time:** Minutes (very fast!)

### 03-build-blueprint.sh

Builds the blueprint (lean) Docker image.

**Usage:**
```bash
./03-build-blueprint.sh BCF_FILE [OPTIONS]
```

**Options:**
- `--genome GENOME` - Genome build (default: GRCh38)
- `--type TYPE` - Dataset type (default: joint)
- `--tag TAG` - Docker tag (default: auto-generated)
- `--registry REGISTRY` - Docker registry (default: ghcr.io/julius-muller)
- `--push` - Push to registry after build
- `--no-cache` - Build without Docker cache

**Output:**
- Docker image: `ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af010`

**Time:** 5-10 minutes

### 04-build-annotated.sh

Builds the annotated (fat) Docker image with VEP.

**Usage:**
```bash
./04-build-annotated.sh BCF_FILE [OPTIONS]
```

**Options:**
- `--genome GENOME` - Genome build (default: GRCh38)
- `--type TYPE` - Dataset type (default: joint)
- `--vep-version VER` - VEP version (default: 115.2)
- `--annotation-name NAME` - Annotation stash name (default: vep_gnomad)
- `--tag TAG` - Docker tag (default: auto-generated)
- `--registry REGISTRY` - Docker registry (default: ghcr.io/julius-muller)
- `--push` - Push to registry after build
- `--no-cache` - Build without Docker cache

**Output:**
- Docker image: `ghcr.io/julius-muller/vcfstash-annotated:gnomad-grch38-vep115`

**Time:** 30-60 minutes (includes VEP annotation)

### 05-push-images.sh

Pushes Docker images to GitHub Container Registry.

**Usage:**
```bash
./05-push-images.sh [OPTIONS]
```

**Options:**
- `--type TYPE` - Image type: blueprint/annotated/both (default: both)
- `--tag TAG` - Specific tag to push (default: all local tags)
- `--registry REGISTRY` - Docker registry (default: ghcr.io/julius-muller)
- `--dry-run` - Show what would be pushed without pushing

**Examples:**
```bash
# Push all images
./05-push-images.sh --type both

# Push only blueprint images
./05-push-images.sh --type blueprint

# Dry run (test without pushing)
./05-push-images.sh --dry-run
```

## Master BCF Strategy

The key insight is to generate the BCF **once** with minimal filtering, then subset as needed:

```
Master BCF (AF ≥ 0.001 = 0.1%)
    ↓ bcftools view -i 'AF>=0.10'  [fast!]
    ├─→ AF ≥ 0.10 (10%)     → Blueprint docker
    ↓ bcftools view -i 'AF>=0.05'  [fast!]
    ├─→ AF ≥ 0.05 (5%)      → Blueprint docker
    ↓ bcftools view -i 'AF>=0.01'  [fast!]
    └─→ AF ≥ 0.01 (1%)      → Blueprint docker
```

**Why 0.001 (0.1%) as master threshold?**
- Balances file size vs. flexibility
- Captures rare but potentially important variants
- Still much smaller than no filtering
- Easy to subset to common thresholds (1%, 5%, 10%)

## File Size Estimates

| Filter | Variants | BCF Size | Use Case |
|--------|----------|----------|----------|
| AF ≥ 0.001 (0.1%) | ~XX M | ~XX GB | Master file |
| AF ≥ 0.01 (1%) | ~XX M | ~XX GB | Research |
| AF ≥ 0.05 (5%) | ~XX M | ~XX GB | Clinical |
| AF ≥ 0.10 (10%) | ~XX M | ~XX GB | Common variants |

*(Actual sizes depend on gnomAD version and genome build)*

## Reproducibility

All scripts are version-controlled in this repository. To reproduce a build:

1. **Document your build:**
   ```bash
   # Save build info
   cat > build-info.txt <<EOF
   Date: $(date)
   Master BCF: gnomad_v4.1_GRCh38_joint_master_af0001.bcf
   Master BCF MD5: $(md5sum data/*.bcf | head -1)
   vcfstash version: $(git rev-parse HEAD)
   Scripts version: $(git rev-parse HEAD)
   Docker version: $(docker --version)
   bcftools version: $(bcftools --version | head -1)
   Hail version: $(python3 -c 'import hail; print(hail.__version__)')
   EOF
   ```

2. **Commit to your repo:**
   ```bash
   git add build-info.txt
   git commit -m "Document Docker build $(date +%Y-%m-%d)"
   ```

3. **Anyone can reproduce:**
   ```bash
   git checkout <commit-hash>
   ./01-generate-master-bcf.sh [same params]
   ```

## GitHub Actions vs Local

| Aspect | GitHub Actions | Local Build |
|--------|---------------|-------------|
| **Best for** | Testing (chrY) | Production (full genome) |
| **Time limit** | 6 hours max | Unlimited |
| **Cost** | CI minutes | Hardware only |
| **Control** | Limited | Full |
| **Debugging** | Difficult | Easy |
| **Iteration** | Slow (queue) | Fast (immediate) |
| **Hardware** | Fixed | Configurable |

**Recommendation:**
- Use **GitHub Actions** for testing changes with chrY
- Use **local builds** for production full-genome images

## Troubleshooting

### Hail Issues

**Problem:** Hail runs out of memory

**Solution:**
```bash
# Increase Spark memory
export SPARK_DRIVER_MEMORY=16g
export SPARK_EXECUTOR_MEMORY=16g
./01-generate-master-bcf.sh ...
```

### Docker Build Issues

**Problem:** Docker build runs out of space

**Solution:**
```bash
# Clean up Docker
docker system prune -a

# Increase Docker disk space
# Edit Docker Desktop settings → Resources → Disk image size
```

**Problem:** Build fails during annotation

**Solution:**
```bash
# Check VEP is working
docker run --rm ensemblorg/ensembl-vep:release_115.2 vep --help

# Build without cache to force fresh download
./04-build-annotated.sh --no-cache ...
```

### bcftools Issues

**Problem:** bcftools view is slow

**Solution:**
```bash
# Use more threads
./02-subset-af-thresholds.sh --threads 8 ...

# Or install newer bcftools
conda install -c bioconda bcftools=1.20
```

## Support

For issues or questions:
- GitHub Issues: https://github.com/julius-muller/vcfstash/issues
- Check existing issues first
- Include build-info.txt when reporting problems
