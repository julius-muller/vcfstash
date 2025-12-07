# Docker Build Scripts

This directory contains scripts for building and managing VCFcache Docker images.

## Quick Start - Build Final Images

To build and push all production images:

```bash
./build-and-push-final.sh
```

This will:
1. Build blueprint images for AF ≥ 0.10 and AF ≥ 0.01
2. Build base annotated images
3. Annotate with VEP
4. Run tests on each image
5. Push to registry if tests pass

## Options

```bash
./build-and-push-final.sh [OPTIONS]

Options:
  --skip-tests    Skip running tests after build
  --skip-push     Skip pushing images to registry
  --af010-only    Build only AF ≥ 0.10 images
  --af001-only    Build only AF ≥ 0.01 images
  -h, --help      Show help message
```

## Examples

### Build and test locally (don't push):
```bash
./build-and-push-final.sh --skip-push
```

### Build only AF 0.10 images:
```bash
./build-and-push-final.sh --af010-only
```

### Build quickly without tests:
```bash
./build-and-push-final.sh --skip-tests
```

## Final Images

The script produces these images:

### AF ≥ 0.10 (10% - recommended for most users)
- **Blueprint**: `ghcr.io/julius-muller/vcfcache-blueprint:gnomad-grch38-joint-af010`
- **Annotated**: `ghcr.io/julius-muller/vcfcache-annotated:gnomad-grch38-joint-af010-vep115`
- **Latest**: `ghcr.io/julius-muller/vcfcache-annotated:latest` (tagged to AF 0.10)

### AF ≥ 0.01 (1% - larger cache, higher hit rate)
- **Blueprint**: `ghcr.io/julius-muller/vcfcache-blueprint:gnomad-grch38-joint-af001`
- **Annotated**: `ghcr.io/julius-muller/vcfcache-annotated:gnomad-grch38-joint-af001-vep115`

## Requirements

- Docker installed and running
- Access to gnomAD BCF files:
  - `/mnt/data/vcfcache_data/gnomad/gnomad_v4.1_GRCh38_joint_af010.bcf`
  - `/mnt/data/vcfcache_data/gnomad/gnomad_v4.1_GRCh38_joint_af001.bcf`
- VEP cache directory: `/mnt/data/apps/ensembl-vep/115/cachedir`
- Docker registry credentials configured (for push)

## Individual Build Scripts

If you need more control, you can run the individual scripts:

### 1. Blueprint Image
```bash
./03-build-blueprint.sh /path/to/gnomad.bcf --host-network --genome GRCh38 --type joint
```

### 2. Base Annotated Image
```bash
./04a-build-base-image.sh /path/to/gnomad.bcf --host-network -y
```

### 3. Annotate with VEP
```bash
./04b-annotate-and-commit.sh \
  --base-image ghcr.io/julius-muller/vcfcache-annotated:TAG-base \
  --vep-cache-dir /path/to/vep/cache \
  -y
```

## GitHub Actions Alternative

For automated builds on every release, consider setting up GitHub Actions. However, note:

**Pros**:
- Automated on every release/tag
- Reproducible builds
- No local resources needed

**Cons**:
- Large gnomAD files (~10GB+) difficult to cache
- Long build times may hit GitHub Actions limits
- VEP cache (~20GB) too large for GitHub

**Recommendation**: Use local builds for production images, GitHub Actions for development/testing with smaller datasets.
