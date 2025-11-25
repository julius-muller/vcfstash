# VCFstash Scripts

This directory contains utility scripts for working with VCFstash.

## export_gnomad_hail.py

Export gnomAD variants directly from the Hail table to BCF format. This is much more efficient than downloading full VCF files.

### Prerequisites

```bash
# Install Hail
pip install hail

# Install bcftools (for BCF conversion)
# Ubuntu/Debian:
apt-get install bcftools

# macOS:
brew install bcftools
```

**Note**: The first time you run the script, it will automatically download the Google Cloud Storage connector JAR (~30MB). This is a one-time setup that takes a minute or two.

### Usage Examples

**No authentication required!** The script uses the public gnomAD bucket by default.

#### Export common variants (AF >= 10%)

```bash
python scripts/export_gnomad_hail.py \
  --output gnomad_af0.10.bcf \
  --af 0.10
```

#### Export only chromosome 1 with AF >= 5%

```bash
python scripts/export_gnomad_hail.py \
  --output gnomad_chr1_af0.05.bcf \
  --af 0.05 \
  --chr chr1
```

#### Export multiple chromosomes

```bash
python scripts/export_gnomad_hail.py \
  --output gnomad_chr1-2_af0.01.bcf \
  --af 0.01 \
  --chr chr1 chr2
```

#### Export all chromosomes (genome-wide)

```bash
python scripts/export_gnomad_hail.py \
  --output gnomad_all_af0.10.bcf \
  --af 0.10
```

### Advanced: Using Requester-Pays Bucket

If you need to use the requester-pays bucket (e.g., for newer versions not yet in the public bucket):

```bash
# Authenticate with gcloud
gcloud auth application-default login

# Export with requester-pays
python scripts/export_gnomad_hail.py \
  --output gnomad_af0.10.bcf \
  --af 0.10 \
  --use-requester-pays \
  --gcs-project your-project-id
```

### Performance Tips

1. **Start with single chromosome**: Test with `--chr chr1` first before running genome-wide
2. **Use higher AF thresholds**: AF >= 0.10 results in much smaller files than AF >= 0.01
3. **Monitor costs**: Hail queries on GCS are billed as requester-pays

### Expected Output Sizes

| AF Threshold | Chromosomes | Approx. Variants | Approx. Size |
|--------------|-------------|------------------|--------------|
| 0.10         | chr1        | ~50,000         | ~10 MB       |
| 0.05         | chr1        | ~200,000        | ~40 MB       |
| 0.01         | chr1        | ~1,000,000      | ~200 MB      |
| 0.10         | all         | ~1,200,000      | ~240 MB      |
| 0.01         | all         | ~20,000,000     | ~4 GB        |

### Using with VCFstash

Once you have the BCF file, use it to build a cache:

```bash
# Initialize cache
vcfstash stash-init \
  --vcf gnomad_af0.10.bcf \
  --output /path/to/cache \
  -y params.yaml

# Annotate the cache
vcfstash stash-annotate \
  --name vep_gnomad \
  --db /path/to/cache \
  -a annotation.config

# Use the cache to annotate samples
vcfstash annotate \
  -a /path/to/cache/stash/vep_gnomad \
  --vcf sample.vcf.gz \
  --output results \
  -y params.yaml
```

### Troubleshooting

#### "Permission denied" errors

The default public bucket should work without authentication. If you see permission errors:
1. Make sure you're not using `--use-requester-pays` flag
2. Check your internet connection (GCS buckets are publicly accessible)
3. Try the requester-pays bucket with authentication (see Advanced section above)

#### Out of memory errors

Hail can be memory-intensive. If you encounter OOM errors:
1. Reduce the scope (fewer chromosomes or higher AF threshold)
2. Increase memory allocation: `export PYSPARK_SUBMIT_ARGS='--driver-memory 8g pyspark-shell'`
3. Use a machine with more RAM

#### Slow export performance

- First run will be slower (Hail needs to read the table schema)
- Subsequent runs benefit from caching
- Consider using Google Compute Engine for faster access to GCS data
