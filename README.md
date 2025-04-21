# VCFstash

A tool to accelerate VCF annotations of large VCF files by maintaining a cache of frequently shared variants across human WGS samples.

## Overview

VCFstash manages a variant cache database and runs VCF annotations only on novel variants not present in the cache. This significantly reduces annotation time for WGS samples since many variants are commonly shared between individuals.

## Key Features

- **Speed**: Typically 2-5x faster than raw VCF annotation by caching common variants
- **Flexibility**: Works with any annotation pipeline (VEP, SnpEff, ANNOVAR, custom scripts)
- **Simplicity**: Easy integration with existing pipelines
- **Efficiency**: Automatic variant normalization and deduplication
- **Scalability**: Multi-processing enabled workflows

## Requirements

- Python 3.11+
- bcftools (included in `tools/` directory)
- Nextflow (included in `workflow/.nextflow/framework/24.10.5/`)

## Installation

```bash
# Clone the repository
git clone https://github.com/julius-muller/vcfstash.git
cd vcfstash

# Create and activate virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Optional: For Parquet support
pip install -r requirements_parquet.txt
```

## Quick Start

VCFstash workflow consists of four main commands:

1. **Initialize a cache database** with common variants
2. **Add more variants** to the cache (optional)
3. **Annotate the cache** with your annotation pipeline
4. **Annotate your samples** using the cache

```bash
# 1. Initialize cache with common variants (e.g., gnomAD)
./vcfstash.py stash-init --vcf gnomad.bcf --output my_cache -y params.yaml

# 2. Add more variants (optional)
./vcfstash.py stash-add --db my_cache -i more_variants.bcf

# 3. Create an annotated cache
./vcfstash.py stash-annotate --name my_annotation --db my_cache -a annotation.config -y params.yaml

# 4. Annotate your samples using the cache
./vcfstash.py annotate -a my_cache/stash/my_annotation --vcf sample.bcf --output results -y params.yaml
```

## Integrating with Your Existing Pipeline

VCFstash makes it easy to integrate with your existing annotation pipeline by splitting it into two parts:

1. **annotation.config**: Contains your annotation commands
2. **params.yaml**: Contains configurable parameters (paths, resources)

### Example: Converting an Existing Pipeline

Let's say you have an existing VEP annotation pipeline:

```bash
vep --offline --buffer_size 500000 --fork 4 --cache \
    --dir_cache /path/to/vep_cache --fasta /path/to/reference.fasta \
    -i sample.vcf -o annotated.vcf --format vcf \
    --transcript_version --symbol --canonical
```

To use this with VCFstash, you need to split it into two files:

#### Step 1: Create annotation.config (contains the command structure)

```javascript
// annotation.config - Contains the annotation command structure
params {
    annotation_cmd = """
      ${params.bcftools_cmd} view \${INPUT_BCF} |
      ${params.annotation_tool_cmd} \
        --offline \
        --buffer_size ${params.vep_buffer} \
        --fork ${params.vep_forks} \
        --cache \
        --dir_cache ${params.vep_cache} \
        --fasta ${params.reference} \
        --format vcf \
        -i stdin \
        -o stdout \
        --transcript_version \
        --symbol \
        --canonical \
        | ${params.bcftools_cmd} view -o \${OUTPUT_BCF} -Ob --write-index
    """

    must_contain_info_tags = [ 'CSQ' ]
}
```

#### Step 2: Create params.yaml (contains the configurable values)

```yaml
# params.yaml - Contains the configurable values

# Tool paths and commands
annotation_tool_cmd: "vep"
bcftools_cmd: "${VCFSTASH_ROOT}/tools/bcftools"

# Reference data
reference: "/path/to/reference.fasta"
reference_md5sum: "28a3d9f0162be1d5db2011aa30458129"

# Resources
vep_cache: "/path/to/vep_cache"
vep_buffer: 500000
vep_forks: 4
```

This separation allows you to:
1. Keep your annotation logic fixed in annotation.config
2. Easily change paths and resources in params.yaml when running in different environments

## Database Structure

```
my_cache/                      # Main cache directory
├── blueprint/                 # Blueprint database
│   ├── vcfstash.bcf           # Normalized variant cache
│   ├── vcfstash.bcf.csi       # Index file
│   └── sources.info           # Database info and logs
├── stash/                     # Annotation instances
│   └── my_annotation/         # Named annotation instance
│       ├── annotation.config  # Frozen annotation config
│       └── ...                # Annotation results
└── workflow/                  # Workflow files
```

## Command Reference

### stash-init

Initialize a new cache database with common variants.

```bash
./vcfstash.py stash-init --vcf <input.bcf> --output <cache_dir> -y <params.yaml> [-f]
```

### stash-add

Add more variants to an existing cache.

```bash
./vcfstash.py stash-add --db <cache_dir> -i <input.bcf>
```

### stash-annotate

Create an annotated cache using your annotation pipeline.

```bash
./vcfstash.py stash-annotate --name <annotation_name> --db <cache_dir> -a <annotation.config> -y <params.yaml> [-f]
```

### annotate

Annotate a sample using the cache.

```bash
./vcfstash.py annotate -a <cache_dir>/stash/<annotation_name> --vcf <sample.bcf> --output <results_dir> -y <params.yaml> [-f]
```

## Performance Tips

- Use a large, comprehensive variant database (like gnomAD) for initialization
- Adjust `vep_forks` and `vep_buffer` in params.yaml based on your system resources
- For large cohorts, the speedup increases with each additional sample

## License

MIT
