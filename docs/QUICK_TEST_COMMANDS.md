# Quick Test Commands for Each Tier 🚀

## **Tier 1: Pre-built Cache Images** ⚡
*Fastest way to test VCFcache (30 seconds)*

### Test with minimal sample data:
```bash
# 1. Pull the image (one-time, ~500MB)
docker pull ghcr.io/julius-muller/vcfcache-cache:GRCh38-af0.10-vep115.2

# 2. Create test VCF
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
  ghcr.io/julius-muller/vcfcache-cache:GRCh38-af0.10-vep115.2 \
  annotate \
  -a /cache/cache/vep_gnomad \
  --vcf /data/sample.vcf.gz \
  --output /results

# 5. Check results
ls -la results/
zcat results/sample_vst.vcf.gz | head -20
```

**Expected result:** Annotated VCF in `results/sample_vst.vcf.gz` with VEP annotations!

---

## **Tier 2: Custom Recipe Production** 🏭
*Full production setup (2-3 hours for setup, then fast)*

### Quick setup with docker for VEP:
```bash
# 1. Setup directories
mkdir -p data/{references,vep_cache,vcfs} results
cd data

# 2. Get VEP image
docker pull ensemblorg/ensembl-vep:release_115.2

# 3. Download reference (5 min download)
cd references
wget -q http://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz  
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
cd ..

# 4. Install VEP cache (10-30 min)
docker run --rm -v $(pwd):/data ensemblorg/ensembl-vep:release_115.2 \
  INSTALL.pl -a cf -s homo_sapiens -y GRCh38 -c /data/vep_cache/115

# 5. Create test population VCF
cd vcfs
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

# 6. Clone VCFcache
git clone https://github.com/julius-muller/vcfcache.git
cd vcfcache

# 7. Build cache (5-10 min)
vcfcache blueprint-init \
  --vcf ../data/vcfs/test_pop.vcf.gz \
  --output ../data/vcfcache_cache \
  -y recipes/hg38_vep115_complete/params.yaml

vcfcache cache-build \
  --name vep_test \
  --db ../data/vcfcache_cache \
  -a recipes/hg38_vep115_complete/annotation.config

# 8. Test annotation
vcfcache annotate \
  -a ../data/vcfcache_cache/cache/vep_test \
  --vcf ../data/vcfs/test_pop.vcf.gz \
  --output ../results \
  -y recipes/hg38_vep115_complete/params.yaml

# 9. Check results
ls -la ../results/
zcat ../results/test_pop_vst.vcf.gz | grep "^1" | head -5
```

**Expected result:** Custom VEP-annotated cache with full CSQ annotations!

---

## **Tier 3: Development/Local** 👨‍💻
*Development environment (15 minutes)*

### Local installation test:
```bash
# 1. Clone and setup
git clone https://github.com/julius-muller/vcfcache.git
cd vcfcache

# 2. Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# 3. Install VCFcache
pip install -e .

# 4. Quick test with built-in test data
vcfcache blueprint-init \
  --vcf tests/data/nodata/gnomad_test.bcf \
  --output /tmp/test_cache \
  -y tests/config/test_params.yaml

vcfcache cache-build \
  --name test_anno \
  --db /tmp/test_cache \
  -a tests/config/test_annotation.config

# 5. Test annotation
vcfcache annotate \
  -a /tmp/test_cache/cache/test_anno \
  --vcf tests/data/nodata/gnomad_test.bcf \
  --output /tmp/results \
  -y tests/config/test_params.yaml

# 6. Check results
ls -la /tmp/results/
bcftools view /tmp/results/gnomad_test_vst.bcf | head -20

# 7. Run tests
python -m pytest tests/ -v
```

**Expected result:** All tests pass, annotated test file with MOCK_ANNO tags!

---

## **Performance Comparison** 📊

Test the same 1000-variant VCF across tiers:

### **Tier 1** (Pre-built): 
- **Setup time**: 30 seconds (just pull image)
- **Annotation time**: ~5 seconds
- **Cache coverage**: Contains gnomAD chr1 variants at specified AF threshold
- **Cache hits**: Varies by sample origin (30-80% typical for population studies)

### **Tier 2** (Custom):  
- **Setup time**: 2-3 hours (one-time)
- **Annotation time**: ~10 seconds (first run), ~5 seconds (subsequent)
- **Cache hits**: Depends on your input data

### **Tier 3** (Development):
- **Setup time**: 15 minutes
- **Annotation time**: ~2 seconds (test data)
- **Cache hits**: ~100% (test data designed for this)

---

## **Troubleshooting Quick Fixes** 🔧

### **Docker Issues:**
```bash
# Permission denied
sudo chmod 666 /var/run/docker.sock

# Out of space
docker system prune -a

# Can't pull image
docker login ghcr.io
```

### **VCF Format Issues:**
```bash
# Fix VCF format
bcftools view input.vcf -Oz -o input.vcf.gz
tabix -p vcf input.vcf.gz
```

### **Missing Tools:**
```bash
# Ubuntu/Debian
sudo apt install bcftools samtools tabix

# macOS
brew install bcftools samtools htslib
```

### **Path Issues:**
```bash
# Make sure files exist
ls -la /path/to/your/file

# Use absolute paths
realpath your_file.vcf.gz
```