# VCFcache Quick Start with Pre-built Cache

Get started with VCFcache in minutes using our pre-built cache images.

## Option 1: Download Pre-built Cache Image

### Pull the Cache Image
```bash
# Pull the latest GRCh38 cache with VEP 115.2 (1% AF threshold)
docker pull ghcr.io/julius-muller/vcfcache-cache:GRCh38-af0.01-vep115.2

# Or pull a smaller cache with higher AF threshold (10% AF)
docker pull ghcr.io/julius-muller/vcfcache-cache:GRCh38-af0.10-vep115.2
```

### Annotate Your Samples Immediately
```bash
# Create results directory
mkdir -p results

# Run annotation on your VCF
docker run --rm \
  -v $(pwd)/your_sample.vcf.gz:/data/input.vcf.gz \
  -v $(pwd)/results:/results \
  ghcr.io/julius-muller/vcfcache-cache:GRCh38-af0.01-vep115.2 \
  annotate \
  -a /cache/cache/vep_gnomad \
  --vcf /data/input.vcf.gz \
  --output /results \
  -y /cache/params.yaml
```

Your annotated VCF will be in `results/input_vst.vcf.gz`!

## Option 2: Build Your Own Cache

### Setup Data Directory
```bash
mkdir -p data/{references,vep_cache,vcfs}
cd data
```

### Download Reference Genome
```bash
cd references
wget http://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
cd ..
```

### Download VEP Cache
```bash
docker run --rm -v $(pwd):/data ensemblorg/ensembl-vep:release_115.2 \
  INSTALL.pl -a cf -s homo_sapiens -y GRCh38 -c /data/vep_cache/115
```

### Get gnomAD Data (Optional)
```bash
cd vcfs
# Download gnomAD exomes (example - adjust URL for latest version)
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr1.vcf.bgz
cd ..
```

### Build Cache
```bash
# Clone VCFcache
git clone https://github.com/julius-muller/vcfcache.git
cd vcfcache

# Use the production recipe
vcfcache blueprint-init \
  --vcf ../data/vcfs/gnomad.exomes.v4.1.sites.chr1.vcf.bgz \
  --output ../data/cache \
  -y recipes/hg38_vep115_complete/params.yaml

vcfcache cache-build \
  --name vep_gnomad \
  --db ../data/cache \
  -a recipes/hg38_vep115_complete/annotation.config
```

## Available Cache Images

| Image Tag | Description | Size | AF Threshold |
|-----------|-------------|------|--------------|
| `GRCh38-af0.10-vep115.2` | Small cache, common variants only | ~500MB | 10% |
| `GRCh38-af0.05-vep115.2` | Medium cache, good coverage | ~1GB | 5% |
| `GRCh38-af0.01-vep115.2` | Large cache, comprehensive | ~3GB | 1% |

## What's Included

Each pre-built cache includes:
- ✅ **VEP 115.2 annotations** (SIFT, PolyPhen, canonical transcripts)
- ✅ **gnomAD population frequencies**
- ✅ **GRCh38 reference genome**
- ✅ **Ready-to-use configuration**
- ✅ **Complete VCFcache CLI**

## Performance

With pre-built caches, annotation performance varies based on cache overlap:
- **Cached variants**: Retrieved instantly from cache (no VEP needed)  
- **Novel variants**: Annotated with full VEP pipeline
- **Cache contains**: gnomAD v4.1 variants from chr1 at specified AF thresholds
- **Typical performance**: 2-5x speedup when 30-80% of variants are cached
- **Best results**: Population studies with common variants

## Next Steps

- See [recipes/](../recipes/) for more annotation configurations
- Read [WIKI.md](../WIKI.md) for advanced usage
- Check [examples/](../examples/) for real-world use cases