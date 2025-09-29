# **VCFstash – Fast, Flexible, and Reliable Variant Annotation Caching**

Stop re-annotating the same common variants over and over.  
`VCFstash` builds a **local, shareable cache** of already-annotated alleles and
lets your preferred tool (VEP, ANNOVAR, SnpEff, …) skip straight to the novel
ones. **Save more than 70% of your annotation time** with minimal changes to your existing pipeline.

---

## ✨ Features

- **Speed**: Typically reduces annotation run time by more than 70% by caching common variants
- **Flexibility**: Works with any annotation tool (VEP, SnpEff, ANNOVAR, custom scripts)
- **Simplicity**: Easy integration with existing pipelines
- **Efficiency**: Automatic variant normalization and deduplication
- **Portability**: Shareable caches for easy collaboration

---

## 📋 Basic Commands

VCFstash operates on indexed VCF/BCF files and provides four main commands:

1. `stash-init`: Create a new cache from a representative VCF
2. `stash-add`: Add more variants to an existing cache
3. `stash-annotate`: Run your annotation tool once on the cache
4. `annotate`: Annotate your samples using the cache

Get help on any command:
```bash
vcfstash -h
vcfstash <command> -h
```

## 🚀 Quick Start

### Installation Options

#### 🐳 Option 1: Pre-built Cache Images (Fastest)

Get started immediately with our pre-built cache images:

```bash
# Pull a ready-to-use cache with VEP annotations
docker pull ghcr.io/julius-muller/vcfstash-cache:GRCh38-af0.01-vep115.2

# Annotate your VCF immediately (70-90% faster!)
docker run --rm \
  -v $(pwd):/data \
  ghcr.io/julius-muller/vcfstash-cache:GRCh38-af0.01-vep115.2 \
  annotate -a /cache/stash/vep_gnomad --vcf your_sample.vcf.gz --output results
```

✨ **Available cache variants:** 10% AF (small), 5% AF (medium), 1% AF (comprehensive)

#### 🐳 Option 2: Using Docker

Build your own cache with Docker:

```bash
docker run --rm -v $(pwd):/data vcfstash --help
```

#### Option 3: Local Installation

If you prefer a local installation, you'll need:

- **Python 3.11+** 
- **Java 17+** (for Nextflow)
- **tabix** for VCF/BCF handling

```bash
# Clone the repository
git clone https://github.com/julius-muller/vcfstash.git
cd vcfstash

# Create a virtual environment and install vcfstash using uv
uv venv .venv
source .venv/bin/activate
uv pip install -e .
```

### Converting Your Existing Workflow

Here's how to convert your existing annotation workflow to use VCFstash:

![Workflow Conversion](resources/conv_uawf.png)

### 3 Simple Steps to Accelerate Your Annotations

1. **Create a blueprint** from a representative VCF (e.g., gnomAD)
   ```bash
   vcfstash stash-init -i input_db.vcf.gz -o /my/cache_dir -y params.yaml
   ```

2. **Run your annotation tool once** on the blueprint
   ```bash
   vcfstash stash-annotate -d /my/cache_dir -n vep_gnomad -a annotation.config
   ```

3. **Annotate your samples** using the cache
   ```bash
   vcfstash annotate -a /my/cache_dir/stash/vep_gnomad --vcf sample1.vcf --output results
   ```

Your annotations will now run significantly faster and your output vcf at results/sample1_vst.vcf
will be identical to running your annotation workflow specified annotation.config directly on sample.vcf.

## 📚 Documentation

For detailed information on configuration, advanced usage, and best practices, please refer to the [Wiki](WIKI.md).

The Wiki covers:
- Detailed configuration instructions
- Cache structure and management
- Performance optimization
- Docker usage
- Testing and validation
- Annotation tool examples
- Troubleshooting

## 📜 Changelog

See [CHANGELOG.md](CHANGELOG.md) for a detailed list of changes and improvements.

## 🙏 Acknowledgements

- [uv](https://github.com/astral-sh/uv) for fast package management
- [Nextflow](https://www.nextflow.io/) for workflow management
- [bcftools](https://github.com/samtools/bcftools) for efficient variant processing

## License

VCFstash is available under the MIT License - see the LICENSE file for details.
