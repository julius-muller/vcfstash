# VEP 115.2 + GRCh38 Recipe

This recipe provides a complete VEP-based annotation pipeline for GRCh38 using Ensembl VEP 115.2.

## Prerequisites

1. **Docker**: Required for running VEP
2. **Reference genome**: GRCh38 primary assembly
3. **VEP cache**: Ensembl VEP cache version 115 for GRCh38

## Quick Setup

### 1. Download VEP Docker Image
```bash
docker pull ensemblorg/ensembl-vep:release_115.2
```

### 2. Download Reference Genome
```bash
mkdir -p /data/references
cd /data/references
wget http://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

### 3. Download VEP Cache
```bash
mkdir -p /data/vep_cache/115
docker run --rm -v /data:/data ensemblorg/ensembl-vep:release_115.2 \
  INSTALL.pl -a cf -s homo_sapiens -y GRCh38 -c /data/vep_cache/115
```

### 4. Use with VCFcache
```bash
# Initialize cache with gnomAD data
vcfcache blueprint-init \
  --vcf gnomad.exomes.vcf.gz \
  --output /data/vcfcache_cache \
  -y params.yaml

# Annotate the cache
vcfcache cache-build \
  --name vep_gnomad \
  --db /data/vcfcache_cache \
  -a annotation.config

# Annotate your samples
vcfcache annotate \
  -a /data/vcfcache_cache/cache/vep_gnomad \
  --vcf your_sample.vcf.gz \
  --output results \
  -y params.yaml
```

## What This Recipe Provides

- **Comprehensive VEP annotations**: Including SIFT, PolyPhen, canonical transcripts
- **Docker-based**: No complex VEP installation required
- **Production-ready**: Uses official Ensembl VEP image
- **Optimized**: High buffer size and parallel processing for speed
- **Quality checks**: Validates VEP version and reference genome

## Customization

You can modify `annotation.config` to:
- Add more VEP plugins
- Change annotation flags
- Customize output fields
- Add custom annotations

See the [VEP documentation](https://www.ensembl.org/info/docs/tools/vep/index.html) for all available options.