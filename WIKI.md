# VCFstash Wiki

This wiki provides detailed information about advanced usage scenarios, configuration options, and best practices for VCFstash.

## Table of Contents

1. [Quick Start Examples](#quick-start-examples)
2. [Setting Up a gnomAD-based Cache](#setting-up-a-gnomad-based-cache)
3. [Optional Checks](#optional-checks)
4. [Configuration Instructions](#configuration-instructions)
5. [Advanced Configuration Options](#advanced-configuration-options)
6. [Cache Structure](#cache-structure)
7. [Using Docker](#using-docker)
8. [Performance Optimization](#performance-optimization)
9. [Testing & Validation](#testing--validation)
10. [Troubleshooting](#troubleshooting)
11. [Annotation Tool Examples](#annotation-tool-examples)

## Quick Start Examples

### **Tier 1: Pre-built Cache (30 seconds)** ⚡

The fastest way to test VCFstash - uses existing cache images:

```bash
# 1. Pull the image (only needed once)
docker pull ghcr.io/julius-muller/vcfstash-cache:GRCh38-af0.10-vep115.2

# 2. Create test data
mkdir -p test-data
cat > test-data/sample.vcf << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	10001	.	T	C	.	PASS	.
1	10002	.	A	G	.	PASS	.
EOF

# 3. Compress and index
bgzip test-data/sample.vcf
tabix -p vcf test-data/sample.vcf.gz

# 4. Annotate instantly!
docker run --rm \
  -v $(pwd)/test-data:/data \
  -v $(pwd)/results:/results \
  ghcr.io/julius-muller/vcfstash-cache:GRCh38-af0.10-vep115.2 \
  annotate \
  -a /cache/stash/vep_gnomad \
  --vcf /data/sample.vcf.gz \
  --output /results

# 5. Check results
ls -la results/
zcat results/sample_vst.vcf.gz | head -20
```

**Expected result:** Annotated VCF with VEP CSQ tags in `results/sample_vst.vcf.gz`

### **Tier 2: Custom Production Setup (2-3 hours)** 🏭

Full production pipeline with VEP:

```bash
# 1. Setup directories and get VEP
mkdir -p data/{references,vep_cache,vcfs} results
docker pull ensemblorg/ensembl-vep:release_115.2

# 2. Download reference genome (5-10 minutes)
cd data/references
wget http://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz  
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
cd ../..

# 3. Install VEP cache (10-30 minutes)
docker run --rm -v $(pwd)/data:/data ensemblorg/ensembl-vep:release_115.2 \
  INSTALL.pl -a cf -s homo_sapiens -y GRCh38 -c /data/vep_cache/115

# 4. Get sample population data
cd data/vcfs
cat > test_pop.vcf << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=248956422>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	pop1	pop2
1	10001	.	T	C	.	PASS	AF=0.15	GT	0/1	0/0
1	10002	.	A	G	.	PASS	AF=0.25	GT	1/1	0/1
1	10003	.	C	T	.	PASS	AF=0.05	GT	0/0	0/1
EOF
bgzip test_pop.vcf && tabix -p vcf test_pop.vcf.gz
cd ../..

# 5. Setup VCFstash
git clone https://github.com/julius-muller/vcfstash.git
cd vcfstash

# 6. Build cache (5-10 minutes)
vcfstash stash-init \
  --vcf ../data/vcfs/test_pop.vcf.gz \
  --output ../data/vcfstash_cache \
  -y recipes/hg38_vep115_complete/params.yaml

vcfstash stash-annotate \
  --name vep_test \
  --db ../data/vcfstash_cache \
  -a recipes/hg38_vep115_complete/annotation.config

# 7. Test annotation
vcfstash annotate \
  -a ../data/vcfstash_cache/stash/vep_test \
  --vcf ../data/vcfs/test_pop.vcf.gz \
  --output ../results \
  -y recipes/hg38_vep115_complete/params.yaml

# 8. Check results
ls -la ../results/
zcat ../results/test_pop_vst.vcf.gz | grep "^1" | head -5
```

**Expected result:** Full VEP annotations with CSQ tags containing SIFT, PolyPhen, etc.

### **Tier 3: Development/Local Testing (15 minutes)** 👨‍💻

Local development with built-in test data:

```bash
# 1. Clone and setup
git clone https://github.com/julius-muller/vcfstash.git
cd vcfstash
python3 -m venv .venv
source .venv/bin/activate

# 2. Install VCFstash
pip install -e .

# 3. Quick test with built-in test data
vcfstash stash-init \
  --vcf tests/data/nodata/gnomad_test.bcf \
  --output /tmp/test_cache \
  -y tests/config/test_params.yaml

vcfstash stash-annotate \
  --name test_anno \
  --db /tmp/test_cache \
  -a tests/config/test_annotation.config

# 4. Test annotation
vcfstash annotate \
  -a /tmp/test_cache/stash/test_anno \
  --vcf tests/data/nodata/gnomad_test.bcf \
  --output /tmp/results \
  -y tests/config/test_params.yaml

# 5. Check results  
ls -la /tmp/results/
bcftools view /tmp/results/gnomad_test_vst.bcf | head -20

# 6. Run full test suite
python -m pytest tests/ -v
```

**Expected result:** All tests pass, annotated file contains MOCK_ANNO tags

### **Performance Comparison**

Testing same 1000-variant file:

| Tier | Setup Time | Annotation Time | Cache Hits | Use Case |
|------|------------|----------------|------------|----------|
| **Tier 1** | 30 seconds | ~5 seconds | 70-90% | Instant usage, demos |
| **Tier 2** | 2-3 hours | ~10 seconds | Custom | Production workflows |
| **Tier 3** | 15 minutes | ~2 seconds | 100% | Development, testing |

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
     -a example_annotation.config \
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
     -a example_annotation.config \
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
     -a example_annotation.config \
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
  -a example_annotation.config \
  --db ${STASH_DIR}

# Step 3: Add gnomAD genomes data
vcfstash stash-add -v \
  --db ${STASH_DIR} \
  -i /path/to/gnomad.genomes.filtered_${AF}.bcf

# Step 4: Second Annotation
vcfstash stash-annotate -v \
  --name gnomad_genex_${AF} \
  -a example_annotation.config \
  --db ${STASH_DIR}

# Step 5: Add dbSNP data
vcfstash stash-add -v \
  --db ${STASH_DIR} \
  -i /path/to/dbsnp.bcf

# Step 6: Final Annotation
vcfstash stash-annotate -v \
  --name gnomad_complete_${AF} \
  -a example_annotation.config \
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

```bash
# In example_annotation.config
required_tool_version = '115.0'
```
```yaml
# In params.yaml
tool_version_command: "vep | grep -oP \"ensembl-vep\\s+:\\s+\\K\\d+\\.\\d+\""
```

This ensures that the same version of the annotation tool is used, preventing inconsistencies in annotation results.

#### 3. Custom Checks

You can add any custom checks needed for your specific workflow:

```yaml
optional_checks:
  reference_md5sum: "28a3d9f0162be1d5db2011aa30458129"
  vep_cache_version: "115"
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

## Configuration Instructions

VCFstash works with **any annotation tool** by wrapping your existing command. Just split it into two parts:

1. **annotation.config**: Your command structure which will be applied to cache *and* input vcf files. This file is fixed at stash-annotate and unavailable at vcfstash annotate! 
2. **params.yaml**: Configurable values (paths, resources) which can be changed for each run

### Example: From VEP command to VCFstash

#### Original command:

```bash
   vep \
   --offline \
   --buffer_size 500000 \
   --fork 4 \
   --cache \
   --dir_cache /path/to/vep_cache/115/cachedir \
   --fasta /path/to/reference.fasta \
   -i sample.vcf \
   -o annotated.vcf \
   --format vcf \
   --canonical
```

### Step 1: Create annotation.config

As a first step, the original annotation command needs to be adapted to the VCFstash format and copied to the .config file.

Conventions here are:
1. The annotation tool command should be replaced by ${params.annotation_tool_cmd} and the actual command should be listed in the params.yaml file.
2. The input filename has to be replaced with the variable ${INPUT_BCF}
3. The output filename hat to be replaced with the variable ${OUTPUT_BCF}

```bash
   ${params.annotation_tool_cmd} \
   --offline \
   --buffer_size 500000 \
   --fork 4 \
   --cache \
   --dir_cache /path/to/vep_cache/115/cachedir \
   --fasta /path/to/reference.fasta \
   -i ${INPUT_BCF} \
   -o ${OUTPUT_BCF} \
   --format vcf \
   --canonical 
```

If there are any parts of the command that need to be kept configurable during vcfstash annotate (e.g., paths, parameters), they have to be replaced with the corresponding variable names starting with `params.`, and listed under 'OPTIONAL RESOURCES' within the params.yaml file:

```bash
   ${params.annotation_tool_cmd} \
   --offline \
   --buffer_size ${params.vep_buffer} \
   --fork ${params.vep_forks} \
   --cache \
   --dir_cache ${params.vep_cache}/115/cachedir \
   --fasta ${params.reference} \
   -i ${INPUT_BCF} \
   -o ${OUTPUT_BCF} \
   --format vcf \
   --canonical 
```

If the annotation tool (here vep) doesn't support bcf as input, the input file needs to be converted and ideally piped through bcftools view. Similarly, the output is expected to be an indexed bcf file. If the tool does not natively support such output, it can be piped through bcftools view with option -W for conversion. Instead of piping, the conversion could be done in multiple steps, however this is not recommended as it would require writing the intermediate files to disk. 

```bash
annotation_cmd = """
   ${params.bcftools_cmd} view ${INPUT_BCF} 
   | ${params.annotation_tool_cmd} \
   --offline \
   --buffer_size ${params.vep_buffer} \
   --fork ${params.vep_forks} \
   --cache \
   --dir_cache ${params.vep_cachedir}/115/cachedir \
   --fasta ${params.reference} \
   -i stdin \
   -o stdout \
   --format vcf \
   --canonical \
   | ${params.bcftools_cmd} view -o ${OUTPUT_BCF} -Ob --write-index
    """
```

Finally the expected tag in the output file needs to be specified using `must_contain_info_tag`:
```bash
    must_contain_info_tag = 'CSQ'
```

#### **Important requirements for annotation.config:**

- **${INPUT_BCF}**: This variable must be used as the input source. It represents an indexed BCF file that VCFstash will provide to your annotation pipeline.
- **${params.bcftools_cmd}**: Use this to convert or pipe the input as needed for your annotation tool.
- **${OUTPUT_BCF}**: Your command must output an indexed BCF file using this variable name. The final output must be in BCF format with an index.

The example above shows a typical pattern: read from ${INPUT_BCF}, pipe through your annotation tool, and write the result to ${OUTPUT_BCF} with indexing.

#### Variable Substitution in Annotation Commands
When writing your annotation command in `annotation.config`, you can use the following variables that will be automatically replaced with the correct paths: 
- Path to the input BCF file `INPUT_BCF`
- Path where the annotated BCF file should be written `OUTPUT_BCF`
- Directory where additional files can be written (auxiliary subdirectory of cwd) `AUXILIARY_DIR`

These variables can be used in either format:
- Plain format: `INPUT_BCF`, `OUTPUT_BCF`, `AUXILIARY_DIR`
- Shell variable format: `${INPUT_BCF}`, `${OUTPUT_BCF}`, `${AUXILIARY_DIR}`

### Step 2: Create params.yaml

As a second step, the parameters need to be defined in the params.yaml file. The YAML file is structured into three main sections:

#### 1. REQUIRED RESOURCES

These are essential parameters that must be defined and are used by the core functionality:

```yaml
## * REQUIRED RESOURCES *
# Do not change the key names in this section, as they are used in the code.

# Tool paths and commands
annotation_tool_cmd: "vep"
tool_version_command: "vep | grep -oP \"ensembl-vep\\s+:\\s+\\K\\d+\\.\\d+\""

# bcftools - best to leave as default since this is the tested and shipped version v1.20
bcftools_cmd: "${VCFSTASH_ROOT}/tools/bcftools"

# Reference data required for the normalization step
reference: "/path/to/reference.fasta"

# Mapping of chromosome names between the reference genome and the VCF file
chr_add: "${VCFSTASH_ROOT}/resources/chr_add.txt"

# Temporary directory for storing intermediate files
temp_dir: "/tmp"
```

#### 2. OPTIONAL RESOURCES

These parameters can be customized for each annotation run and are referenced in your annotation.config file:

```yaml
## * OPTIONAL RESOURCES *
# All keys here can be utilized in the example_annotation.config file as ${params.MYKEY}
vep_buffer: 500000
vep_forks: 4
vep_cache: "/path/to/vep_cachedir"
```

#### 3. OPTIONAL CHECKS

These parameters provide verification mechanisms to ensure consistency between cache creation and annotation:

```yaml
## * OPTIONAL CHECKS *
# Optional checks and verifications. These can be set to any key:value pair, 
# but all keys must match the values in the example_annotation.config.
optional_checks:
  reference_md5sum: "28a3d9f0162be1d5db2011aa30458129"
  # Add other optional verifications here
  # example_version: "1.2.3"
  # example_option: "value"
```

The final example of the params.yaml can be found in the repository: [test_params.yaml](https://github.com/julius-muller/vcfstash/blob/main/tests/config/test_params.yaml) and for annotation.config: [test_annotation.config](https://github.com/julius-muller/vcfstash/blob/main/tests/config/test_annotation.config).

> **Important**: It is the responsibility of the user that any optional resource listed here impacting annotation results is never changed between cache creation and annotation. The optional checks section helps enforce this consistency.

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
annotation_tool_cmd: "docker run --user $(id -u):$(id -g) -i -v /mnt/data:/mnt/data --rm ensemblorg/ensembl-vep:release_115.2 vep"
tool_version_command: "docker run --user $(id -u):$(id -g) -i -v /mnt/data:/mnt/data --rm ensemblorg/ensembl-vep:release_115.2 vep | grep -oP \"ensembl-vep\\s+:\\s+\\K\\d+\\.\\d+\""
```

## Cache Structure

The VCFstash cache is organized in a structured directory hierarchy that maintains both the normalized variants and their annotations. Understanding this structure can help you manage and troubleshoot your caches.

### Directory Structure

```
cache_directory/
├── blueprint/                  # Contains normalized variants
│   ├── vcfstash.bcf            # Normalized variants from input VCF files
│   ├── vcfstash.bcf.csi        # Index for the normalized variants
│   ├── sources.info            # Information about input VCF files
│   └── ...                     # Nextflow reports and logs
│
├── stash/                      # Contains annotated variants
│   └── [annotation_name]/      # Named annotation directory (e.g., "vep_gnomad")
│       ├── annotation.config   # Locked annotation configuration
│       ├── blueprint_snapshot.info  # Blueprint info at annotation time
│       ├── vcfstash_annotated.bcf   # Annotated variants
│       ├── vcfstash_annotated.bcf.csi  # Index for annotated variants
│       ├── annotation.yaml         # Configuration for annotation tools
│       └── ...                 # Nextflow reports and logs
│
├── workflow/                   # Contains Nextflow workflow files
│   ├── init.yaml               # Initial configuration
│   ├── main.nf                 # Nextflow workflow file
│   ├── modules/                # Nextflow modules
│   └── ...                     # Workflow logs
│
└── vcfdb.log                   # Log file for cache operations
```

### Key Files Explained

#### Blueprint Files

- **vcfstash.bcf**: The core normalized variant database containing all variants from input VCF files. This file is created during `stash-init` and updated during `stash-add`.
- **sources.info**: JSON file tracking all input VCF files that have been added to the cache, including their MD5 checksums and timestamps.

#### Stash Files

- **annotation.yaml**: Contains configuration for annotation tools, including paths, commands, and resource settings.
- **annotation.config**: The locked annotation configuration that defines exactly how variants are annotated. This file is frozen after `stash-annotate` to ensure consistency.
- **blueprint_snapshot.info**: Records the state of the blueprint at the time of annotation, ensuring traceability.
- **vcfstash_annotated.bcf**: The annotated variant database that serves as the cache for future annotations.

#### Workflow Files

- **init.yaml**: Initial configuration used by the Nextflow workflow.
- **main.nf**: The Nextflow workflow that orchestrates the normalization and annotation processes.

### Cache Lifecycle

1. **Initialization** (`stash-init`): Creates the blueprint directory and normalizes input variants
2. **Addition** (`stash-add`): Updates the blueprint with additional variants
3. **Annotation** (`stash-annotate`): Creates a named annotation in the stash directory
4. **Usage** (`annotate`): Uses the annotated cache to speed up annotation of new samples

Understanding this structure helps when:
- Troubleshooting annotation issues
- Managing multiple annotation versions
- Sharing caches between environments
- Backing up or archiving caches

The cache is designed to be portable - you can copy an entire cache directory to another location or system and use it there, as long as the paths in your `params.yaml` file are updated accordingly.

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
docker-compose -f docker/docker-compose.yml build
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
# Create example_annotation.config
cat > data/example_annotation.config << 'EOF'
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
  -a /data/example_annotation.config \
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

## Testing & Validation

VCFstash includes comprehensive tests to ensure reliability and correctness:

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run all tests
python -m pytest

# Run specific test with increased verbosity
python -m pytest -xvs tests/test_core.py
```

### Test Structure

The test suite is organized into two modules:

- `test_core.py`: Tests basic utility functions like MD5 calculation
- `test_annotate.py`: Tests the annotate command and full workflow

Each test prints information about which part of the code is being tested, making it easier to understand test coverage.

### Test Implementation

The tests are designed to run on any system after installation of the package, without requiring external annotation tools. They use:

- `tests/config/test_params.yaml`: Configuration for the test annotation tool (bcftools)
- `tests/config/test_annotation.config`: Configuration for the test annotation command

The tests create temporary directories for output and clean up after themselves. They use bcftools (included in the package) to simulate annotations, making the tests portable and reliable.

### Running Tests

To run the tests:

```bash
# Run all tests
python -m pytest

# Run a specific test file
python -m pytest tests/test_annotate.py

# Run a specific test with increased verbosity
python -m pytest -xvs tests/test_annotate.py::test_annotate_workflow
```

All tests should pass on any system where the package is installed, without requiring any external tools or configuration.

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
    required_tool_version = '115.2'
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
