# VCFstash Wiki

This wiki provides detailed information about advanced usage scenarios, configuration options, and best practices for VCFstash.

## Table of Contents

1. [Setting Up a gnomAD-based Cache](#setting-up-a-gnomad-based-cache)
2. [Optional Checks](#optional-checks)
3. [Advanced Configuration Options](#advanced-configuration-options)
4. [Using Docker](#using-docker)
5. [Performance Optimization](#performance-optimization)
6. [Troubleshooting](#troubleshooting)
7. [Annotation Tool Examples](#annotation-tool-examples)

## Setting Up a gnomAD-based Cache

gnomAD (Genome Aggregation Database) is an excellent source of common variants for initializing your VCFstash cache. This section provides a step-by-step guide to setting up a comprehensive gnomAD-based cache.

### Prerequisites

- gnomAD VCF files (exomes and/or genomes)
- Sufficient disk space (varies based on filtering criteria)
- VCFstash installed and configured

### Basic Workflow

Here's a basic workflow for setting up a gnomAD-based cache:

1. **Initialize the cache with gnomAD exomes data**:
   ```bash
   vcfstash stash-init \
     --vcf /path/to/gnomad.exomes.vcf.gz \
     --output /path/to/cache_dir \
     -y params.yaml
   ```

2. **Annotate the exomes data**:
   ```bash
   vcfstash stash-annotate \
     --name gnomad_ex \
     -a annotation.config \
     --db /path/to/cache_dir
   ```

3. **Add gnomAD genomes data to the cache**:
   ```bash
   vcfstash stash-add \
     --db /path/to/cache_dir \
     -i /path/to/gnomad.genomes.vcf.gz
   ```

4. **Re-annotate the combined data**:
   ```bash
   vcfstash stash-annotate \
     --name gnomad_genex \
     -a annotation.config \
     --db /path/to/cache_dir
   ```

5. **Optionally add other sources (e.g., dbSNP)**:
   ```bash
   vcfstash stash-add \
     --db /path/to/cache_dir \
     -i /path/to/dbsnp.vcf.gz
   ```

6. **Re-annotate the final combined data**:
   ```bash
   vcfstash stash-annotate \
     --name gnomad_genex_dbsnp \
     -a annotation.config \
     --db /path/to/cache_dir
   ```

### Filtering by Allele Frequency

For optimal performance, consider filtering gnomAD data by allele frequency. This reduces the cache size while still capturing the most common variants:

1. **Filter gnomAD files by allele frequency**:
   ```bash
   bcftools view -i 'AF>=0.01' /path/to/gnomad.vcf.gz -Ob -o gnomad_af0.01.bcf
   bcftools index gnomad_af0.01.bcf
   ```

2. **Initialize cache with filtered data**:
   ```bash
   vcfstash stash-init \
     --vcf gnomad_af0.01.bcf \
     --output /path/to/cache_dir \
     -y params.yaml
   ```

Common allele frequency thresholds:
- 0.1 (10%): Very common variants, smallest cache size
- 0.05 (5%): Good balance of coverage and size
- 0.01 (1%): More comprehensive coverage, larger cache
- 0.001 (0.1%): Very comprehensive, but much larger cache

### Complete Example Script

Here's a complete script for setting up a gnomAD-based cache with different allele frequency thresholds:

```bash
#!/bin/bash
set -euo pipefail

# Define allele frequency threshold
AF="0.01"

# Define the base stash directory
STASH_DIR="/path/to/gnomad_${AF}"

# Step 1: Stash Initialization with gnomAD exomes
vcfstash stash-init -v \
  --vcf /path/to/gnomad.exomes.filtered_${AF}.bcf \
  --output ${STASH_DIR} \
  -y params.yaml

# Step 2: First Annotation
vcfstash stash-annotate -v \
  --name gnomad_ex_${AF} \
  -a annotation.config \
  --db ${STASH_DIR}

# Step 3: Add gnomAD genomes data
vcfstash stash-add -v \
  --db ${STASH_DIR} \
  -i /path/to/gnomad.genomes.filtered_${AF}.bcf

# Step 4: Second Annotation
vcfstash stash-annotate -v \
  --name gnomad_genex_${AF} \
  -a annotation.config \
  --db ${STASH_DIR}

# Step 5: Add dbSNP data
vcfstash stash-add -v \
  --db ${STASH_DIR} \
  -i /path/to/dbsnp.bcf

# Step 6: Final Annotation
vcfstash stash-annotate -v \
  --name gnomad_complete_${AF} \
  -a annotation.config \
  --db ${STASH_DIR}

echo "Cache setup complete: ${STASH_DIR}"
```

## Optional Checks

VCFstash provides a mechanism for optional checks to ensure consistency between cache creation and annotation. These checks help prevent issues that might arise from environment changes or inconsistent configurations.

### Available Optional Checks

#### 1. Reference Genome MD5 Checksum

The most common check is verifying the reference genome's MD5 checksum:

```yaml
optional_checks:
  reference_md5sum: "28a3d9f0162be1d5db2011aa30458129"
```

This ensures that the same reference genome is used for both cache creation and annotation, preventing mismatches in variant coordinates or reference alleles.

#### 2. Tool Version Checks

You can verify that the annotation tool version matches the expected version:

```yaml
# In annotation.config
required_tool_version = '113.0'

# In params.yaml
tool_version_command: "vep | grep -oP \"ensembl-vep\\s+:\\s+\\K\\d+\\.\\d+\""
```

This ensures that the same version of the annotation tool is used, preventing inconsistencies in annotation results.

#### 3. Custom Checks

You can add any custom checks needed for your specific workflow:

```yaml
optional_checks:
  reference_md5sum: "28a3d9f0162be1d5db2011aa30458129"
  vep_cache_version: "113"
  genome_build: "GRCh37"
```

These values must match exactly between the params.yaml file and the annotation.config file.

### How Optional Checks Work

1. During `stash-annotate`, VCFstash reads the optional_checks from params.yaml and stores them in the cache
2. During `annotate`, VCFstash compares the current params.yaml values with the stored values
3. If any values don't match, VCFstash raises an error, preventing inconsistent annotations

### Best Practices for Optional Checks

1. **Always include reference_md5sum**: This is the most critical check to ensure consistent variant coordinates
2. **Include tool version checks**: Especially important for tools like VEP where different versions can produce different annotations
3. **Document your checks**: Add comments explaining what each check verifies and why it's important
4. **Be specific**: Use precise version numbers and checksums rather than ranges or patterns

## Advanced Configuration Options

VCFstash offers several advanced configuration options for customizing its behavior to suit your specific needs.

### Environment Variables

VCFstash uses the `VCFSTASH_ROOT` environment variable to locate resources. This is automatically set by the application, but you can override it if needed:

```bash
export VCFSTASH_ROOT=/custom/path/to/vcfstash
```

### Nextflow Resource Configuration

You can customize CPU, memory, and executor settings in a Nextflow configuration file:

```groovy
// Process configuration
process {
    executor = 'slurm'  // Use SLURM for job submission
    cpus = 8            // Default CPUs for all processes
    memory = '16 GB'    // Default memory for all processes

    // Process-specific settings
    withName: 'RenameAndNormalizeVCF' {
        memory = '20 GB'   // More memory for this specific process
    }

    withName: 'RunAnnotation' {
        cpus = 16          // More CPUs for annotation
        memory = '32 GB'    
    }
}
```

### Docker/Singularity Integration

VCFstash works well with containerized annotation tools. Example using Docker with VEP:

```yaml
# In params.yaml
annotation_tool_cmd: "docker run --user $(id -u):$(id -g) -i -v /mnt/data:/mnt/data --rm ensemblorg/ensembl-vep:release_113.0 vep"
tool_version_command: "docker run --user $(id -u):$(id -g) -i -v /mnt/data:/mnt/data --rm ensemblorg/ensembl-vep:release_113.0 vep | grep -oP \"ensembl-vep\\s+:\\s+\\K\\d+\\.\\d+\""
```

## Using Docker

VCFstash provides Docker support for easy deployment and consistent execution across different environments. This section covers various ways to use VCFstash with Docker.

### Basic Docker Usage

#### Building the Docker Image

You can build the Docker image from the repository:

```bash
# Clone the repository
git clone https://github.com/julius-muller/vcfstash.git
cd vcfstash

# Build the Docker image
docker build -f docker/Dockerfile -t vcfstash .
```

#### Running VCFstash with Docker

Once the image is built, you can run VCFstash commands:

```bash
# Run VCFstash with the help command
docker run --rm vcfstash --help

# Mount volumes for data access
docker run --rm \
  -v /path/to/reference:/reference:ro \
  -v /path/to/data:/data \
  -v /path/to/cache:/cache \
  vcfstash <command> <options>
```

### Using Docker Compose

VCFstash includes a Docker Compose configuration for easier management of volumes and environment variables.

#### Running with Docker Compose

```bash
# Navigate to the repository
cd vcfstash

# Set environment variables (optional)
export REFERENCE_DIR=/path/to/reference
export DATA_DIR=/path/to/data
export CACHE_DIR=/path/to/cache

# Run VCFstash with Docker Compose
docker-compose -f docker/docker-compose.yml run --rm vcfstash <command> <options>
```

#### Using the Helper Script

VCFstash provides a helper script (`run-vcfstash.sh`) that simplifies running with Docker Compose:

```bash
# Make the script executable
chmod +x docker/run-vcfstash.sh

# Run VCFstash using the helper script
./docker/run-vcfstash.sh <command> <options>

# With custom directories
REFERENCE_DIR=/path/to/reference DATA_DIR=/path/to/data CACHE_DIR=/path/to/cache \
  ./docker/run-vcfstash.sh <command> <options>
```

### Volume Mounting

When using Docker with VCFstash, you need to mount volumes for:

1. **Reference data** (read-only): Contains reference genomes and other reference files
2. **Input/output data**: Contains input VCF files and output directories
3. **Cache directory**: Stores the VCFstash cache

#### Example with Explicit Volume Mounts

```bash
docker run --rm \
  -v /path/to/reference:/reference:ro \
  -v /path/to/data:/data \
  -v /path/to/cache:/cache \
  vcfstash stash-init \
  --vcf /data/gnomad.vcf.gz \
  --output /cache \
  -y /data/params.yaml
```

#### Adjusting params.yaml for Docker

When using Docker, update your `params.yaml` to use paths inside the container:

```yaml
# params.yaml for Docker
reference: "/reference/GRCh38.fa"
temp_dir: "/data/temp"
```

### Complete Examples

#### Example 1: Initialize a Cache with Docker

```bash
# Create directories
mkdir -p reference data cache

# Copy files
cp /path/to/reference.fasta reference/
cp /path/to/gnomad.vcf.gz data/
cp /path/to/params.yaml data/

# Run stash-init
./docker/run-vcfstash.sh stash-init \
  --vcf /data/gnomad.vcf.gz \
  --output /cache \
  -y /data/params.yaml
```

#### Example 2: Annotate with Docker

```bash
# Create annotation.config
cat > data/annotation.config << 'EOF'
params {
    annotation_cmd = """
        ${params.annotation_tool_cmd} \
        --offline \
        --buffer_size ${params.vep_buffer} \
        --fork ${params.vep_forks} \
        --cache \
        --dir_cache ${params.vep_cache} \
        --fasta ${params.reference} \
        -i ${INPUT_BCF} \
        -o ${OUTPUT_BCF} \
        --format vcf \
        --canonical
    """
    must_contain_info_tag = 'CSQ'
}
EOF

# Run stash-annotate
./docker/run-vcfstash.sh stash-annotate \
  --name vep_gnomad \
  -a /data/annotation.config \
  --db /cache
```

#### Example 3: Annotate a Sample with Docker

```bash
# Copy sample VCF
cp /path/to/sample.vcf.gz data/

# Run annotate
./docker/run-vcfstash.sh annotate \
  -a /cache/stash/vep_gnomad \
  --vcf /data/sample.vcf.gz \
  --output /data/results \
  -y /data/params.yaml
```

### Docker Compose Configuration

You can customize the Docker Compose configuration by creating your own `docker-compose.yml` file:

```yaml
version: '3.8'

services:
  vcfstash:
    image: vcfstash
    volumes:
      - /path/to/reference:/reference:ro
      - /path/to/data:/data
      - /path/to/cache:/cache
    environment:
      - REFERENCE_PATH=/reference
    command: --help
```

### Parquet Output

VCFstash supports converting the final BCF file to Parquet format for efficient querying with tools like DuckDB:

```bash
vcfstash annotate -a /path/to/cache/stash/my_annotation \
    --vcf sample.vcf.gz \
    --output results \
    -y params.yaml \
    --parquet
```

## Performance Optimization

### Cache Size vs. Coverage Tradeoffs

The size of your cache affects both storage requirements and annotation speed:

| Allele Frequency | Approximate Size | Coverage | Use Case |
|------------------|------------------|----------|----------|
| 10% (0.1)        | Smallest         | Common variants only | Quick testing, limited storage |
| 5% (0.05)        | Small            | Most common variants | Good balance for most uses |
| 1% (0.01)        | Medium           | Many variants | Production use with adequate storage |
| 0.1% (0.001)     | Large            | Most variants | Comprehensive coverage, requires significant storage |

### Hardware Recommendations

For optimal performance:

1. **Storage**: Use SSD storage for the cache directory to maximize I/O performance
2. **Memory**: Allocate at least 8GB RAM for basic usage, 16-32GB for large files
3. **CPU**: Multi-core processors benefit annotation tools like VEP that support parallelization
4. **Network**: If using shared storage, ensure high-bandwidth, low-latency connections

### Parallelization Strategies

1. **Tool-level parallelization**: Configure your annotation tool to use multiple threads
   ```yaml
   vep_forks: 8  # For VEP
   ```

2. **Nextflow parallelization**: Configure Nextflow to use multiple CPUs
   ```groovy
   process {
       cpus = 8
   }
   ```

3. **Sample-level parallelization**: Process multiple samples in parallel using a workflow manager

## Troubleshooting

### Common Issues and Solutions

#### 1. Missing or Incompatible Reference Genome

**Symptoms**: Errors about reference allele mismatches or failed normalization

**Solution**: 
- Verify the reference genome path in params.yaml
- Check that the reference_md5sum matches
- Ensure the reference genome is properly indexed

#### 2. Annotation Tool Errors

**Symptoms**: Annotation process fails with tool-specific errors

**Solution**:
- Check tool installation and dependencies
- Verify tool version matches required_tool_version
- Test the annotation command directly outside VCFstash

#### 3. Cache Consistency Issues

**Symptoms**: Errors about mismatched optional checks

**Solution**:
- Ensure params.yaml values match those used during cache creation
- Check for changes in reference genome or tool versions
- Recreate the cache if necessary

#### 4. Performance Problems

**Symptoms**: Slow annotation or high resource usage

**Solution**:
- Filter input data by allele frequency
- Adjust CPU and memory allocation in Nextflow config
- Use SSD storage for cache directory
- Check for network bottlenecks if using shared storage

### Logging and Debugging

VCFstash provides detailed logging to help diagnose issues:

1. **Increase verbosity**: Use `-v` or `-vv` for more detailed logs
   ```bash
   vcfstash -vv annotate -a /path/to/cache ...
   ```

2. **Debug mode**: Use `--debug` to keep intermediate files
   ```bash
   vcfstash --debug annotate -a /path/to/cache ...
   ```

3. **Log files**: Check log files in the output directory
   - `vcfdb.log`: Main log file
   - `*.bcf.log`: Process-specific logs

4. **Nextflow reports**: Enable Nextflow reports for detailed execution information
   ```bash
   vcfstash annotate -a /path/to/cache ... --report --timeline --dag
   ```

## Annotation Tool Examples

VCFstash works with any annotation tool that can process VCF files. This section provides examples for popular annotation tools.

### VEP (Variant Effect Predictor)

#### params.yaml
```yaml
## * REQUIRED RESOURCES *
annotation_tool_cmd: "vep"
tool_version_command: "vep | grep -oP \"ensembl-vep\\s+:\\s+\\K\\d+\\.\\d+\""
bcftools_cmd: "${VCFSTASH_ROOT}/tools/bcftools"
reference: "/path/to/reference.fasta"
chr_add: "${VCFSTASH_ROOT}/resources/chr_add.txt"
temp_dir: "/tmp"

## * OPTIONAL RESOURCES *
vep_buffer: 500000
vep_forks: 4
vep_cache: "/path/to/vep_cache"

## * OPTIONAL CHECKS *
optional_checks:
  reference_md5sum: "28a3d9f0162be1d5db2011aa30458129"
```

#### annotation.config
```
params {
    annotation_cmd = """
        ${params.bcftools_cmd} view ${INPUT_BCF} | \
        ${params.annotation_tool_cmd} \
        --offline \
        --buffer_size ${params.vep_buffer} \
        --fork ${params.vep_forks} \
        --cache \
        --dir_cache ${params.vep_cache} \
        --fasta ${params.reference} \
        -i stdin \
        -o stdout \
        --format vcf \
        --canonical | \
        ${params.bcftools_cmd} view -o ${OUTPUT_BCF} -Ob --write-index
    """

    must_contain_info_tag = 'CSQ'
    required_tool_version = '113.0'
}

optional_checks {
    reference_md5sum = '28a3d9f0162be1d5db2011aa30458129'
}
```

### ANNOVAR

#### params.yaml
```yaml
## * REQUIRED RESOURCES *
annotation_tool_cmd: "/path/to/annovar/table_annovar.pl"
tool_version_command: "/path/to/annovar/table_annovar.pl | grep -oP 'Version: \\K[\\d\\.]+''"
bcftools_cmd: "${VCFSTASH_ROOT}/tools/bcftools"
reference: "/path/to/reference.fasta"
chr_add: "${VCFSTASH_ROOT}/resources/chr_add.txt"
temp_dir: "/tmp"

## * OPTIONAL RESOURCES *
annovar_db: "/path/to/annovar/humandb"
annovar_buildver: "hg38"
annovar_protocol: "refGene,clinvar_20220320,gnomad211_genome"
annovar_operation: "g,f,f"

## * OPTIONAL CHECKS *
optional_checks:
  reference_md5sum: "28a3d9f0162be1d5db2011aa30458129"
  annovar_buildver: "hg38"
```

#### annotation.config
```
params {
    annotation_cmd = """
        # Convert BCF to VCF for ANNOVAR
        ${params.bcftools_cmd} view ${INPUT_BCF} -Ov -o ${AUXILIARY_DIR}/input.vcf

        # Run ANNOVAR
        ${params.annotation_tool_cmd} \
        ${AUXILIARY_DIR}/input.vcf \
        ${params.annovar_db} \
        -buildver ${params.annovar_buildver} \
        -out ${AUXILIARY_DIR}/output \
        -remove \
        -protocol ${params.annovar_protocol} \
        -operation ${params.annovar_operation} \
        -nastring . \
        -vcfinput

        # Convert back to BCF
        ${params.bcftools_cmd} view ${AUXILIARY_DIR}/output.${params.annovar_buildver}_multianno.vcf \
        -Ob -o ${OUTPUT_BCF} --write-index
    """

    must_contain_info_tag = 'ANNOVAR_DATE'
    required_tool_version = '2020-06-08'
}

optional_checks {
    reference_md5sum = '28a3d9f0162be1d5db2011aa30458129'
    annovar_buildver = 'hg38'
}
```

### SnpEff

#### params.yaml
```yaml
## * REQUIRED RESOURCES *
annotation_tool_cmd: "snpEff"
tool_version_command: "snpEff -version | grep -oP 'SnpEff \\K[\\d\\.]+''"
bcftools_cmd: "${VCFSTASH_ROOT}/tools/bcftools"
reference: "/path/to/reference.fasta"
chr_add: "${VCFSTASH_ROOT}/resources/chr_add.txt"
temp_dir: "/tmp"

## * OPTIONAL RESOURCES *
snpeff_datadir: "/path/to/snpeff/data"
snpeff_genome: "GRCh38.105"
snpeff_config: "/path/to/snpeff/snpEff.config"

## * OPTIONAL CHECKS *
optional_checks:
  reference_md5sum: "28a3d9f0162be1d5db2011aa30458129"
  snpeff_genome: "GRCh38.105"
```

#### annotation.config
```
params {
    annotation_cmd = """
        ${params.annotation_tool_cmd} \
        -dataDir ${params.snpeff_datadir} \
        -config ${params.snpeff_config} \
        -nodownload \
        -v \
        ${params.snpeff_genome} \
        ${INPUT_BCF} | \
        ${params.bcftools_cmd} view -o ${OUTPUT_BCF} -Ob --write-index
    """

    must_contain_info_tag = 'ANN'
    required_tool_version = '5.1'
}

optional_checks {
    reference_md5sum = '28a3d9f0162be1d5db2011aa30458129'
    snpeff_genome = 'GRCh38.105'
}
```

### Custom Script Example

This example shows how to use a custom Python script for annotation:

#### params.yaml
```yaml
## * REQUIRED RESOURCES *
annotation_tool_cmd: "python3 /path/to/custom_annotator.py"
tool_version_command: "python3 /path/to/custom_annotator.py --version"
bcftools_cmd: "${VCFSTASH_ROOT}/tools/bcftools"
reference: "/path/to/reference.fasta"
chr_add: "${VCFSTASH_ROOT}/resources/chr_add.txt"
temp_dir: "/tmp"

## * OPTIONAL RESOURCES *
custom_db: "/path/to/custom/database.db"
custom_threads: 4

## * OPTIONAL CHECKS *
optional_checks:
  reference_md5sum: "28a3d9f0162be1d5db2011aa30458129"
  custom_db_version: "1.2.3"
```

#### annotation.config
```
params {
    annotation_cmd = """
        ${params.annotation_tool_cmd} \
        --input ${INPUT_BCF} \
        --output ${OUTPUT_BCF} \
        --database ${params.custom_db} \
        --threads ${params.custom_threads} \
        --write-index
    """

    must_contain_info_tag = 'CUSTOM_ANNO'
    required_tool_version = '1.0.0'
}

optional_checks {
    reference_md5sum = '28a3d9f0162be1d5db2011aa30458129'
    custom_db_version = '1.2.3'
}
```

### Minimal Example

Here's a minimal example that adds a simple annotation using bcftools:

#### params.yaml
```yaml
## * REQUIRED RESOURCES *
annotation_tool_cmd: "${VCFSTASH_ROOT}/tools/bcftools"
tool_version_command: "${VCFSTASH_ROOT}/tools/bcftools --version-only"
bcftools_cmd: "${VCFSTASH_ROOT}/tools/bcftools"
reference: "/path/to/reference.fasta"
chr_add: "${VCFSTASH_ROOT}/resources/chr_add.txt"
temp_dir: "/tmp"

## * OPTIONAL RESOURCES *
annotation_file: "/path/to/annotations.vcf.gz"

## * OPTIONAL CHECKS *
optional_checks:
  reference_md5sum: "28a3d9f0162be1d5db2011aa30458129"
```

#### annotation.config
```
params {
    annotation_cmd = """
        ${params.annotation_tool_cmd} annotate \
        -a ${params.annotation_file} \
        -c INFO \
        ${INPUT_BCF} \
        -Ob -o ${OUTPUT_BCF} --write-index
    """

    must_contain_info_tag = 'AF'
    required_tool_version = '1.20'
}

optional_checks {
    reference_md5sum = '28a3d9f0162be1d5db2011aa30458129'
}
```
