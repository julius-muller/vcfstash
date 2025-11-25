# gnomAD BCF Generation with Hail

This document explains the smart caching system for gnomAD BCF files used by VCFstash.

## Overview

VCFstash uses gnomAD data to pre-populate variant caches. Instead of downloading huge VCF files every time, we:

1. **Query gnomAD directly** using Hail from the public bucket (no auth required!)
2. **Pre-compute BCF files** for common configurations
3. **Cache them as GitHub Releases** for fast reuse
4. **Docker builds** download pre-computed files when available

## Architecture

### Two-Tier Caching System

```
┌─────────────────────────────────────────┐
│ Tier 1: Pre-computed BCF Files          │
│ (GitHub Releases)                        │
│ - Fast download (~30 seconds)            │
│ - Common configurations (AF 0.10, chr1)  │
│ - Updated weekly                         │
└─────────────────────────────────────────┘
                  ↓ (if available)
┌─────────────────────────────────────────┐
│ Tier 2: On-demand Hail Query             │
│ (GitHub Actions)                         │
│ - Runs only if no pre-computed file      │
│ - Custom configurations                  │
│ - Takes 10-30 minutes                    │
└─────────────────────────────────────────┘
```

## GitHub Actions Workflows

### 1. Precompute gnomAD BCF (`.github/workflows/precompute_gnomad_bcf.yml`)

**Purpose:** Pre-compute BCF files for common configurations

**Triggers:**
- Manual: Configure AF threshold and chromosomes
- Weekly: Auto-generates AF=0.10, all chromosomes

**Output:** Creates GitHub Release with BCF files as assets

**Example Usage:**
```bash
# Manual trigger via GitHub UI
Actions → Precompute gnomAD BCF → Run workflow
  AF: 0.10
  Chromosomes: all
  gnomAD version: 4.1
```

### 2. Build Cache Image with Hail (`.github/workflows/build_cache_hail.yml`)

**Purpose:** Build Docker images with gnomAD cache

**Smart Caching:**
1. Checks if pre-computed BCF exists in releases
2. If yes → Downloads (~30 seconds)
3. If no → Runs Hail query (~10-30 minutes)
4. Builds Docker image with VCFstash + gnomAD cache

**Triggers:**
- Manual: Custom configurations
- Weekly: Builds latest cache

## Local Usage

### Running Hail Export Locally

```bash
# Install dependencies
pip install hail
sudo apt-get install bcftools

# Export gnomAD data (no authentication needed!)
python scripts/export_gnomad_hail.py \
  --output gnomad_chr1_af0.10.bcf \
  --af 0.10 \
  --chr chr1

# Or all chromosomes
python scripts/export_gnomad_hail.py \
  --output gnomad_af0.10.bcf \
  --af 0.10

# Or multiple specific chromosomes
python scripts/export_gnomad_hail.py \
  --output gnomad_chr1-2_af0.05.bcf \
  --af 0.05 \
  --chr chr1 chr2
```

### Using Pre-computed BCF Files

```bash
# Download from releases
wget https://github.com/YOUR_USER/vcfstash/releases/download/gnomad-bcf-v4.1-af0.10-genome/gnomad_v4.1_GRCh38_af0.10.bcf
wget https://github.com/YOUR_USER/vcfstash/releases/download/gnomad-bcf-v4.1-af0.10-genome/gnomad_v4.1_GRCh38_af0.10.bcf.csi

# Use with VCFstash
vcfstash stash-init \
  --vcf gnomad_v4.1_GRCh38_af0.10.bcf \
  --output /path/to/cache \
  -y params.yaml
```

## Configuration Options

### Allele Frequency Thresholds

| AF Threshold | Variant Count (chr1) | File Size | Use Case |
|--------------|---------------------|-----------|----------|
| 0.10 (10%)   | ~50K                | ~10 MB    | Very common variants |
| 0.05 (5%)    | ~200K               | ~40 MB    | Common variants |
| 0.01 (1%)    | ~1M                 | ~200 MB   | Moderately common |

### Chromosome Options

- **Single**: `chr1` - Fast testing
- **Multiple**: `chr1 chr2` - Specific regions
- **All**: `all` - Genome-wide (recommended for production)

## Technical Details

### How Anonymous Access Works

The scripts use Google Cloud Storage (GCS) connector with anonymous access:

```python
# Disable authentication for public buckets
os.environ['GCE_METADATA_HOST'] = 'metadata.google.internal.invalid'

spark_conf = {
    'spark.hadoop.google.cloud.auth.service.account.enable': 'false',
    'spark.hadoop.google.cloud.auth.null.enable': 'true',
}
```

This allows reading from `gs://gcp-public-data--gnomad/` without any GCP credentials!

### BCF File Naming Convention

```
gnomad_v{VERSION}_{GENOME}_{CHROMOSOMES}_af{AF}.bcf
```

Examples:
- `gnomad_v4.1_GRCh38_af0.10.bcf` - All chromosomes, AF ≥ 10%
- `gnomad_v4.1_GRCh38_chr1_af0.05.bcf` - Chromosome 1 only, AF ≥ 5%
- `gnomad_v4.1_GRCh38_chr1_chr2_af0.01.bcf` - Chr 1 & 2, AF ≥ 1%

### Release Tag Convention

```
gnomad-bcf-v{VERSION}-af{AF}-{CHROMOSOMES}
```

Examples:
- `gnomad-bcf-v4.1-af0.10-genome` - All chromosomes
- `gnomad-bcf-v4.1-af0.05-chr1` - Chromosome 1 only

## Performance Comparison

### Docker Build Times

| Scenario | Time | Network | Notes |
|----------|------|---------|-------|
| Pre-computed (chr1) | ~2 min | 10 MB | Downloads from releases |
| Pre-computed (genome) | ~5 min | 240 MB | Downloads from releases |
| Hail query (chr1) | ~15 min | 50 MB | Queries gnomAD directly |
| Hail query (genome) | ~30 min | 200 MB | Queries gnomAD directly |

### Cache Hit Rates

Expected cache hit rates for different populations:

| Population | chr1 (10% AF) | Genome (10% AF) | Genome (1% AF) |
|------------|---------------|-----------------|----------------|
| European | ~75% | ~80% | ~90% |
| African | ~60% | ~70% | ~85% |
| East Asian | ~70% | ~75% | ~88% |

## Maintenance

### Weekly Updates

The `schedule` trigger runs weekly to:
1. Generate latest BCF files with current gnomAD data
2. Update GitHub Releases
3. Rebuild Docker images

This ensures users always have access to the latest gnomAD data.

### Storage Costs

All BCF files are stored as GitHub Release assets:
- **Free** for public repositories
- **Counted towards storage quota** for private repositories
- Typical usage: ~500 MB per configuration

## Troubleshooting

### Hail Java Version Warning

```
Hail was built with Java 11. You are using Java 21
```

**Solution:** This is a warning, not an error. Hail works fine with newer Java versions for basic operations.

### GCS Connection Timeout

If you see timeouts connecting to GCS:
1. Check internet connection
2. Verify the public bucket path is correct
3. Ensure firewall isn't blocking GCS access

### BCF File Too Large

If the BCF file is larger than expected:
1. Increase AF threshold (e.g., 0.05 → 0.10)
2. Reduce chromosome scope
3. Consider splitting by chromosome

## Future Enhancements

Potential improvements:
1. **Multi-population filtering**: Filter by specific populations
2. **Gene-based filtering**: Include only coding regions
3. **Compressed storage**: Use additional compression
4. **CDN distribution**: Faster downloads via CDN
5. **Incremental updates**: Only download changes
