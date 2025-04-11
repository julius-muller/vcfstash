Here's a README that explains the project's purpose and usage:

# VCFstash

A tool to accelerate VCF annotations of large VCF files by maintaining a cache of frequently shared variants across human WGS samples.

## Overview

VCFstash manages a variant cache database and runs VCF annotations only on novel variants not present in the cache. This significantly reduces annotation time for WGS samples since many variants are commonly shared between individuals.

## Features

- Initialize bare bone 'blueprint' variant cache databases from user provided VCF/BCF files. Typically from large human genome allele frequency databases like gnomAD
- Extend exisiting blueprint caches by adding new variants e.g. add variants based on WES to WGS
- Instanciate a cache by annotation of a blueprint cache.
- Run VCF annotations on arbitrary user provided vcf files, leveraging overlap of cached variants resulting in a large performance boost
- Automatic variant normalization and deduplication
- Resource management through nextflow config files
- Multi-processing enabled workflows
- Detailed logging and workflow tracking
- Support for workflow parameter customization

## Requirements

- bcftools (latest version in path, also located at tools/ )
- Ensembl VCF (configured via user provided nextflow config)
- Python 3.11+

Docker is not yet implemented.
Nextflow installation is optional, as the jar is shipped in:
workflow/.nextflow/framework/24.10.5/nextflow-24.10.5-one.jar and could be invoked with:
java -jar framework/24.10.5/nextflow-24.10.5-one.jar run

## Installation

```bash
git clone https://github.com/julius-muller/vcfstash.git
cd vcfstash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

pytest tests/test_cache_init.py 
```

# *** *THE BELOW IS OUTDATED* ***

## Usage 

### Initialize Cache Database

Create a new cache database from a VCF file (e.g., gnomAD):

```bash
docker compose run vcfstash cache init \
  -n vcf_db_gnomad_41 \
  -i gnomad.vcf.gz \
  -f reference.fasta \
  -t 8
```

### Add Variants

Add variants from another VCF to the cache:

```bash
docker compose run vcfstash cache add \
  -d vcf_db_gnomad_41 \
  -i dbsnp.vcf.gz \
  -f reference.fasta \
  -t 8
```

### Run VCF Annotation

Annotate a sample using the cache:

```bash
docker compose run vcfstash annotate \
  -d vcf_db_gnomad_41 \
  -w workflow \
  -i sample.vcf.gz
```

## Database Structure

```
vcf_db_name/
├── vcf_db.bcf          # Normalized variant cache
├── vcf_db.bcf.csi      # Index file
├── vcf_db.bcf.info     # Database info and logs
└── stash/        # VCF annotation runs
    └── YYYYMMDD_HHMMSS_[hash]/
        ├── workflow/   # Workflow files
        ├── results/    # VCF output
        └── logs/       # Run logs
```

## Workflow Parameters

All workflow parameters are configured in `config/nextflow.yml`:

### VCF Configuration

```yaml
# VCF settings
vcf_assembly: 'GRCh37'  # Genome assembly version
vep_options:            # VCF annotation options
  - '--transcript_version'
  - '--total_length'
  # ... more options

# Performance settings
vep_max_chr_parallel: 2  # Maximum chromosomes to process in parallel
vep_max_forks: 2         # VCF forks per chromosome
```

## Performance

- Uses bcftools for efficient variant handling
- Multi-threaded processing
- Caches common variants to avoid redundant annotations
- Typical speedup: 2-5x compared to raw VCF annotation

## Contributing

Pull requests welcome! Please follow our contribution guidelines.

## License

MIT
```

This README provides a clear overview of the project's purpose, installation steps, usage examples, and key features while keeping a technical focus appropriate for developers.