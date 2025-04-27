# **VCFstash – turbo-charge your variant annotation workflow**

Stop re-annotating the same common variants over and over.  
`VCFstash` builds a **local, shareable cache** of already-annotated alleles and
lets your preferred tool (VEP, ANNOVAR, SnpEff, …) skip straight to the novel
ones. **Save more than 70% of your annotation time** with minimal changes to your existing pipeline.

---

## ✨ Key Features

- **Speed**: Typically reduces annotation run time by more than 70% by caching common variants
- **Flexibility**: Works with any annotation tool (VEP, SnpEff, ANNOVAR, custom scripts)
- **Simplicity**: Easy integration with existing pipelines
- **Efficiency**: Automatic variant normalization and deduplication
- **Reproducibility**: Consistent annotations with hashed and locked cache

---

## 🛠️ Quick Setup

```bash
# Clone and install
git clone https://github.com/julius-muller/vcfstash.git
cd vcfstash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

That's it! VCFstash includes all required tools (bcftools, Nextflow) in the repository.

---

## 🏁 3-step quick-start for setting up a new cache

1. **Create a blueprint** from a representative VCF
   ```bash
   vcfstash.py stash-init -i gnomad_4_wgs.vcf -o /my/destination/dir
   ```

2. **Run your annotation tool once** on the blueprint and freeze it
   ```bash
   vcfstash.py stash-annotate -d /my/destination/dir -n vep_gnomad -a annotation.config -y params.yaml
   ```

3. **Annotate your samples** using the cache
   ```bash
   vcfstash.py annotate -a /my/destination/dir/stash/vep_gnomad --vcf sample.vcf --output results -y params.yaml
   ```

**Optional:** Add more variants to your cache (e.g., merge gnomAD WES into WGS blueprint)
```bash
vcfstash.py stash-add -i gnomad_4_wes.vcf -d /my/destination/dir
```

Your `annotation.config` contains **exactly** the command you already use – see next section for details.

---

## 🔄 Converting your existing annotation command

VCFstash works with **any annotation tool** by wrapping your existing command. Just split it into two parts:

1. **annotation.config**: Your command structure (fixed)
2. **params.yaml**: Configurable values (paths, resources)

### Example: From VEP command to VCFstash

#### Original command:

```bash
vep --offline --buffer_size 500000 --fork 4 --cache \
    --dir_cache /path/to/vep_cache --fasta /path/to/reference.fasta \
    -i sample.vcf -o annotated.vcf --format vcf \
    --transcript_version --symbol --canonical
```

#### Step 1: Create annotation.config

```javascript
// annotation.config - Fixed command structure
params {
    annotation_cmd = """
      ${params.bcftools_cmd} view \${INPUT_BCF} 
      | ${params.annotation_tool_cmd} \
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

**Important requirements for annotation.config:**

- **${INPUT_BCF}**: This variable must be used as the input source. It represents an indexed BCF file that VCFstash will provide to your annotation pipeline.
- **${params.bcftools_cmd}**: Use this to convert or pipe the input as needed for your annotation tool.
- **${OUTPUT_BCF}**: Your command must output an indexed BCF file using this variable name. The final output must be in BCF format with an index.

The example above shows a typical pattern: read from \${INPUT_BCF}, pipe through your annotation tool, and write the result to ${OUTPUT_BCF} with indexing.

#### Step 2: Create params.yaml

```yaml
# params.yaml - Configurable values

# Tool paths
annotation_tool_cmd: "vep"
bcftools_cmd: "${VCFSTASH_ROOT}/tools/bcftools"

# Reference data (edit these!)
reference: "/path/to/reference.fasta"
reference_md5sum: "28a3d9f0162be1d5db2011aa30458129"

# Resources (edit these!)
vep_cache: "/path/to/vep_cache"
vep_buffer: 500000
vep_forks: 4
```

## 📋 Command Reference

VCFstash provides four simple commands that follow the workflow shown in the quick-start:

1. `stash-init`: Create a new cache from a representative VCF (e.g., gnomAD)
2. `stash-add`: Add more variants to an existing cache
3. `stash-annotate`: Run your annotation tool once on the cache
4. `annotate`: Annotate your samples using the cache

Full documentation with all options is available by running:
```bash
./vcfstash.py --help
```

## 🚀 Running annotation anywhere using the cache

Once your cache is set up, you can run the annotation command **anywhere** the cache is available - perfect for distributed computing environments or sharing with collaborators:

```bash
vcfstash.py annotate -a /path/to/cache/stash/my_annotation \
    --vcf sample.vcf \
    --output results \
    -y params.yaml
```

Just make sure your `params.yaml` file is properly configured for the new environment. This portability means you can:

- Run annotations on different compute clusters
- Share pre-built annotation caches with team members
- Integrate with any workflow management system

## 🖥️ Resource Management

VCFstash supports optional resource management through Nextflow configuration. Create a `nextflow.config` file to control CPU, memory, and other resources:

```groovy
// Process configuration
process {
    executor = 'local'
    cpus = 4               // Default CPUs for all processes
    memory = '8 GB'        // Default memory for all processes

    // Process-specific settings
    withName: 'RenameAndNormalizeVCF' {
        memory = '20 GB'   // More memory for this specific process
    }
}
```

This configuration is optional - VCFstash will work with default settings, but customizing resources can improve performance for your specific environment.

## 💡 Performance Tips

- **Start with quality data**: Use a large, comprehensive variant database (like gnomAD) for initialization
- **Keep cache size manageable**: Using only common alleles (e.g., 10% allele frequency or higher) provides excellent performance while avoiding huge caches
- **Benefit from scale**: For large cohorts, the speedup increases with each additional sample - the more you use it, the more time you save
- **Optimize I/O**: Consider using SSD storage for the cache directory to maximize performance
- **Parallelize wisely**: Adjust CPU and memory settings in the Nextflow configuration based on your available resources

## 🧪 Testing & Validation

VCFstash includes comprehensive tests to ensure reliability and correctness:

```bash
# Run all tests
python -m pytest

# Run specific test with increased verbosity
python -m pytest -xvs tests/test_validation.py
```

To validate your own setup or create new test references:

```bash
python tests/update_reference.py --golden
```

This creates reference outputs that serve as the gold standard for tests, helping you verify your installation and configuration are working correctly.

## License

VCFstash is available under the MIT License - see the LICENSE file for details.


