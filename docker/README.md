# VCFstash Docker Images

VCFstash provides two types of Docker images to suit different workflows:

## Image Types

### 🔷 Blueprint (Lean) - `vcfstash-blueprint`
**Data-only image for power users**

Contains:
- Pre-built variant cache (blueprint) from gnomAD
- vcfstash CLI tool
- Minimal dependencies (no VEP, no annotation tools)

Use when:
- You want to use your own annotation tools
- You need maximum flexibility
- You want the smallest image size
- You'll run annotation separately

**Image naming:**
```
ghcr.io/USER/vcfstash-blueprint:<source>-<genome>-<type>-<chr>-af<threshold>

Examples:
ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af010           # Full genome
ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-chry-af010     # chrY test
```

### 🔶 Annotated (Fat) - `vcfstash-annotated`
**All-in-one image with VEP pipeline**

Contains:
- Pre-built AND pre-annotated cache with VEP
- Full VEP installation (ensembl-vep:release_115.2)
- VEP cache (Homo sapiens GRCh38)
- vcfstash CLI tool
- Ready to annotate samples immediately

Use when:
- You want convenience over flexibility
- You're using VEP for annotation
- You need production-ready setup
- You want to skip annotation pipeline setup

**Image naming:**
```
ghcr.io/USER/vcfstash-annotated:<source>-<genome>-vep<version>

Examples:
ghcr.io/julius-muller/vcfstash-annotated:gnomad-grch38-vep115              # Full genome
ghcr.io/julius-muller/vcfstash-annotated:gnomad-grch38-chry-vep115         # chrY test
```

---

## Building Images

### Blueprint Image

```bash
# Trigger via GitHub Actions
gh workflow run build_cache_hail.yml \
  -f af=0.10 \
  -f chromosomes=all \
  -f type=joint \
  -f genome=GRCh38 \
  -f gnomad_version=4.1

# Or build locally
docker build \
  -f docker/Dockerfile.cache-hail \
  --build-arg AF=0.10 \
  --build-arg GENOME=GRCh38 \
  --build-arg CACHE_NAME=gnomad_grch38_joint_af010 \
  --build-arg BCF_FILE=path/to/gnomad.bcf \
  -t vcfstash-blueprint:custom \
  .
```

### Annotated Image

```bash
# Trigger via GitHub Actions
gh workflow run build_annotated_hail.yml \
  -f af=0.10 \
  -f chromosomes=all \
  -f type=joint \
  -f genome=GRCh38 \
  -f gnomad_version=4.1 \
  -f vep_version=115.2 \
  -f annotation_name=vep_gnomad

# Or build locally (requires VEP base image)
docker build \
  -f docker/Dockerfile.annotated \
  --build-arg AF=0.10 \
  --build-arg GENOME=GRCh38 \
  --build-arg CACHE_NAME=gnomad_grch38_joint_af010 \
  --build-arg BCF_FILE=path/to/gnomad.bcf \
  --build-arg ANNOTATION_NAME=vep_gnomad \
  -t vcfstash-annotated:custom \
  .
```

---

## Using Images

### Blueprint Image Usage

```bash
# Pull the image
docker pull ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af010

# Run tests
docker run --rm --entrypoint /bin/sh \
  ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af010 \
  -c 'cd /app && export PYTHONPATH=/app/venv/lib/python3.13/site-packages:$PYTHONPATH && python3 -m pytest -m blueprint tests/ -v'

# Inspect the cache
docker run --rm \
  ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af010 \
  stash-init --help

# Use with your own annotation pipeline
docker run --rm \
  -v /path/to/sample.vcf:/data/sample.vcf \
  -v /path/to/output:/output \
  ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af010 \
  annotate \
    -a /cache/db/stash/YOUR_ANNOTATION \
    --vcf /data/sample.vcf \
    --output /output \
    -y /app/recipes/docker-cache/params.yaml
```

### Annotated Image Usage

```bash
# Pull the image
docker pull ghcr.io/julius-muller/vcfstash-annotated:gnomad-grch38-vep115

# Check version
docker run --rm ghcr.io/julius-muller/vcfstash-annotated:gnomad-grch38-vep115 -v

# Annotate a sample VCF (uses pre-annotated cache)
docker run --rm \
  -v /path/to/sample.vcf:/data/sample.vcf \
  -v /path/to/output:/output \
  ghcr.io/julius-muller/vcfstash-annotated:gnomad-grch38-vep115 \
  annotate \
    -a /cache/db/stash/vep_gnomad \
    --vcf /data/sample.vcf \
    --output /output \
    -y /app/recipes/docker-annotated/params.yaml

# Result: Annotated VCF in /path/to/output with 70%+ speedup!
```

---

## Image Comparison

| Feature | Blueprint (Lean) | Annotated (Fat) |
|---------|-----------------|-----------------|
| **Size** | ~2-3 GB | ~15-20 GB |
| **Build Time** | ~5-10 min | ~30-60 min |
| **Contains VEP** | ❌ No | ✅ Yes |
| **Contains VEP Cache** | ❌ No | ✅ Yes |
| **Pre-annotated** | ❌ No | ✅ Yes |
| **Flexibility** | ✅ High | ⚠️ Medium |
| **Ready to Use** | ⚠️ Setup needed | ✅ Immediate |
| **Use Case** | Custom pipelines | Production VEP |

---

## CI/CD Integration

Both image types are built via GitHub Actions workflows:

- **Blueprint:** `.github/workflows/build_cache_hail.yml`
- **Annotated:** `.github/workflows/build_annotated_hail.yml`

Workflows automatically:
1. Export gnomAD data with Hail (or download from releases)
2. Build Docker image with cache
3. Push to GitHub Container Registry
4. Cache Docker layers for faster rebuilds

---

## Development

### File Structure

```
docker/
├── Dockerfile.cache-hail       # Blueprint image (lean)
├── Dockerfile.annotated        # Annotated image (fat)
├── build-cache-with-bcf.sh     # Blueprint build script
├── build-annotated-cache.sh    # Annotated build script
└── README.md                   # This file

recipes/
├── docker-cache/
│   ├── params.yaml            # Blueprint params
│   └── annotation.config      # Blueprint annotation (not used)
└── docker-annotated/
    ├── params.yaml            # Annotated params (VEP configured)
    └── annotation.config      # VEP annotation config
```

### Testing Blueprint Image

```bash
pytest -m blueprint tests/ -v
```

Only runs tests compatible with blueprint image (no Nextflow/VEP).

### Testing Annotated Image

```bash
pytest tests/ -v
```

Runs all tests including full integration tests with VEP.

---

## Troubleshooting

### Blueprint Image

**Issue:** "Nextflow not found"
- **Solution:** Blueprint is data-only. Use annotated image for full pipeline.

**Issue:** "Cannot run annotation"
- **Solution:** Blueprint requires you to run `stash-annotate` separately with your tools.

### Annotated Image

**Issue:** "VEP cache not found"
- **Solution:** Ensure VEP cache is properly copied in Dockerfile at `/opt/vep/.vep`

**Issue:** "Build takes too long"
- **Solution:** Annotation is slow (30-60 min for full genome). Use chrY for testing first.

**Issue:** "Out of memory during build"
- **Solution:** Increase Docker memory limit or use smaller dataset (chrY, chr22).

---

## Performance

**Blueprint build times:**
- chrY: ~3-5 minutes
- chr22: ~5-10 minutes
- Full genome: ~10-15 minutes

**Annotated build times:**
- chrY: ~10-15 minutes
- chr22: ~20-30 minutes
- Full genome: ~30-60 minutes

**Annotation speedup with cache:**
- 70-90% faster than annotating from scratch
- Varies based on sample overlap with gnomAD
