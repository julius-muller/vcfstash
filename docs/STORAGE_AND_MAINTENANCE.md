# Storage Locations & Maintenance Guide

## 📍 Storage Locations by Tier

### **Tier 1: Pre-built Cache Images** 
**Target Users**: Immediate usage, no setup

**Storage Locations:**
- **Docker Registry**: `ghcr.io/julius-muller/vcfcache-cache:*`
- **Local Docker**: `~/.docker/` (Docker images)
- **Runtime Cache**: `/cache/` (inside container)
- **User Data**: Mounted volumes (your choice)

**What's Stored:**
```
Container /cache/
├── cache/vep_gnomad/          # Pre-built annotated cache
├── params.yaml                # Configuration used
└── reference/                 # Reference genome files
```

### **Tier 2: Custom Recipe Production**
**Target Users**: Production workloads, customized pipelines

**Storage Locations:**
- **VEP Cache**: `/data/vep_cache/115/` (3-5 GB)
- **Reference**: `/data/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa*` (3 GB)
- **VCFcache Cache**: `/data/vcfcache_cache/` (varies by input data)
- **Input Data**: `/data/vcfs/` (gnomAD, etc.)

**Directory Structure:**
```
/data/
├── vep_cache/115/             # VEP annotation cache (3-5 GB)
├── references/                # Reference genomes (3 GB)
│   └── Homo_sapiens.GRCh38.dna.primary_assembly.fa*
├── vcfs/                      # Input VCF files (varies)
│   └── gnomad.*.vcf.gz
└── vcfcache_cache/            # VCFcache cache output
    ├── blueprint/             # Normalized variants
    └── cache/                 # Annotated variants
        └── vep_gnomad/
```

### **Tier 3: Development/Local**
**Target Users**: Developers, testing, custom environments

**Storage Locations:**
- **Code**: `~/vcfcache/` (Git repo)
- **Virtual Env**: `~/vcfcache/.venv/`
- **Test Data**: `~/vcfcache/tests/data/`
- **Local Cache**: User-defined locations

---

## 🔧 Maintenance Checklist

### **Weekly Tasks** (Automated via GitHub Actions)

- [x] **Build fresh cache images** (Sunday 2 AM UTC)
  ```bash
  # Triggered automatically or manual dispatch:
  # Go to: GitHub → Actions → "Build and Publish cache image" → Run workflow
  ```

### **Monthly Tasks**

- [ ] **Update VEP version** (when new releases available)
  - Check: https://github.com/Ensembl/ensembl-vep/releases
  - Update: `recipes/hg38_vep115_complete/params.yaml`
  - Update: `.github/workflows/build_cache.yml`

- [ ] **Update gnomAD data** (when new versions released)
  - Check: https://gnomad.broadinstitute.org/downloads
  - Update GitHub Actions to use new URLs

- [ ] **Check Docker base images**
  - Ubuntu 25.04 security updates
  - VEP Docker image updates

### **Quarterly Tasks**

- [ ] **Reference genome updates** (rare)
  - GRCh38 patches: Check NCBI/Ensembl
  - Update MD5 checksums in configs

- [ ] **Storage cleanup**
  ```bash
  # Clean old Docker images (local development)
  docker system prune -a
  
  # Clean old GitHub packages (keep last 5 versions)
  # Done via GitHub Container Registry settings
  ```

### **As-Needed Tasks**

- [ ] **Recipe updates** (new annotation tools, flags)
- [ ] **Bug fixes** (if users report issues)
- [ ] **Performance optimizations**

---

## 📊 Storage Requirements

| Component | Size | Frequency | Notes |
|-----------|------|-----------|--------|
| **Tier 1 Images** | 500MB-3GB | Weekly builds | Auto-deleted old versions |
| **VEP Cache** | 3-5 GB | Monthly updates | Download once, reuse |
| **Reference Genome** | 3 GB | Rare updates | GRCh38 stable |
| **gnomAD Data** | 100GB-1TB | Quarterly | Optional, user choice |
| **VCFcache Cache** | 1-50 GB | Per project | Depends on input data |

---

## 🚨 Monitoring & Alerts

### **GitHub Actions Health**
- **Monitor**: Weekly cache builds succeed
- **Alert**: Build failures → Check logs
- **Fix**: Usually VEP/network issues

### **Docker Registry Limits**
- **Monitor**: GitHub Container Registry storage
- **Alert**: Approaching limits → Delete old versions
- **Fix**: Keep only latest 3 versions per AF threshold

### **User Feedback Channels**
- **GitHub Issues**: Bug reports, feature requests
- **Discussions**: Usage questions
- **Email**: Critical issues

---

## 🔄 Disaster Recovery

### **If Cache Builds Fail**
1. Check GitHub Actions logs
2. Test VEP Docker image manually
3. Verify gnomAD URLs still valid
4. Rollback to previous working config

### **If Docker Registry Issues**
1. GitHub Container Registry outages → Use docker.io
2. Storage limits → Clean old versions
3. Access issues → Check GitHub tokens

### **If Reference Data Corrupted**
1. Re-download from Ensembl/NCBI
2. Verify MD5 checksums
3. Update all configs with new checksums