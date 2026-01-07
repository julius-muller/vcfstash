# VCFcache Docker

VCFcache ships a single “runtime” container image intended to run the CLI consistently across machines.

## Published image

- `ghcr.io/julius-muller/vcfcache:latest` (runtime)
- `ghcr.io/julius-muller/vcfcache:vX.Y.Z` (versioned release tags)

The image includes `bcftools` and `tabix` and runs `vcfcache` as the entrypoint.

## Typical usage

Annotate a sample with a public cache alias (auto-download from Zenodo):
```bash
docker run --rm -v $(pwd):/work ghcr.io/julius-muller/vcfcache:latest \
  annotate \
    -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic \
    --vcf /work/sample.bcf \
    --output /work/sample_vc.bcf \
    --stats-dir /work/out \
    --force
```

List caches:
```bash
docker run --rm ghcr.io/julius-muller/vcfcache:latest list caches
```

## Building locally

Multi-stage Dockerfile: `docker/Dockerfile.vcfcache`

Build runtime image:
```bash
docker build --target final -f docker/Dockerfile.vcfcache -t vcfcache:local .
```

Build and run tests in the container (test stage):
```bash
docker build --target test -f docker/Dockerfile.vcfcache -t vcfcache:test .
docker run --rm vcfcache:test
```
- chr22: ~20-30 minutes
- Full genome: ~30-60 minutes

**Annotation speedup with cache:**
- 70-90% faster than annotating from scratch
- Varies based on sample overlap with gnomAD
