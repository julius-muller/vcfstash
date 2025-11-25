# BCF Storage Strategy

## Tiered Storage Approach

VCFstash uses a **hybrid storage strategy** to optimize for cost, speed, and convenience.

### Storage Tiers

```
┌─────────────────────────────────────────┐
│ Tier 1: GitHub Releases                 │
│ - Files: chr1, small BCFs (<100MB)      │
│ - Cost: FREE                             │
│ - Speed: Fast                            │
│ - Access: Built-in to workflows         │
└─────────────────────────────────────────┘

┌─────────────────────────────────────────┐
│ Tier 2: Cloud Storage (Public Bucket)   │
│ - Files: Genome-wide (>100MB)           │
│ - Cost: ~$0.02/GB/month                 │
│ - Speed: Very fast (CDN)                │
│ - Access: HTTPS, no auth needed         │
└─────────────────────────────────────────┘

┌─────────────────────────────────────────┐
│ Tier 3: On-Demand Generation            │
│ - Files: Custom configurations          │
│ - Cost: CI minutes only                 │
│ - Speed: 10-30 minutes                  │
│ - Access: Generated as needed           │
└─────────────────────────────────────────┘
```

## File Size Estimates

| Configuration | Variants | File Size | Recommended Storage |
|--------------|----------|-----------|---------------------|
| chr1, AF≥0.10 | ~50K | ~10 MB | GitHub Releases |
| chr1, AF≥0.05 | ~200K | ~40 MB | GitHub Releases |
| chr1, AF≥0.01 | ~1M | ~200 MB | GitHub Releases or Cloud |
| Genome, AF≥0.10 | ~1.2M | ~240 MB | Cloud Storage |
| Genome, AF≥0.05 | ~5M | ~1 GB | Cloud Storage |
| Genome, AF≥0.01 | ~20M | ~4 GB | Cloud Storage |

## GitHub Releases Setup (Current)

### Configuration

Files are automatically uploaded by the workflow:

```yaml
# .github/workflows/precompute_gnomad_bcf.yml
- name: Create GitHub Release
  uses: softprops/action-gh-release@v2
  with:
    tag_name: gnomad-bcf-v4.1-af0.10-chr1
    files: |
      gnomad_v4.1_GRCh38_chr1_af0.10.bcf
      gnomad_v4.1_GRCh38_chr1_af0.10.bcf.csi
```

### Access

```bash
# Download from releases
wget https://github.com/YOUR_USER/vcfstash/releases/download/TAG/file.bcf
```

### Limits

- **File size**: 2GB per file
- **Total storage**: Generous for public repos (no hard limit)
- **Bandwidth**: Unlimited for public repos

## Cloud Storage Setup (Optional)

### Google Cloud Storage (Recommended)

#### 1. Create Public Bucket

```bash
# Install gcloud CLI
curl https://sdk.cloud.google.com | bash

# Authenticate
gcloud auth login

# Create bucket
gsutil mb -p your-project-id -l us-central1 gs://vcfstash-bcf-public

# Make bucket public
gsutil iam ch allUsers:objectViewer gs://vcfstash-bcf-public

# Enable uniform access
gsutil uniformbucketlevelaccess set on gs://vcfstash-bcf-public
```

#### 2. Upload BCF Files

```bash
# Upload with public-read ACL
gsutil -h "Cache-Control:public, max-age=86400" \
  cp gnomad_v4.1_GRCh38_af0.10.bcf \
  gs://vcfstash-bcf-public/v4.1/
```

#### 3. Access Files

```bash
# Direct HTTPS access (no auth needed!)
wget https://storage.googleapis.com/vcfstash-bcf-public/v4.1/gnomad_v4.1_GRCh38_af0.10.bcf
```

### AWS S3 (Alternative)

```bash
# Create public bucket
aws s3 mb s3://vcfstash-bcf-public --region us-east-1

# Upload with public-read ACL
aws s3 cp gnomad_v4.1_GRCh38_af0.10.bcf \
  s3://vcfstash-bcf-public/v4.1/ \
  --acl public-read

# Access at:
# https://vcfstash-bcf-public.s3.amazonaws.com/v4.1/gnomad_v4.1_GRCh38_af0.10.bcf
```

## Workflow Integration

### Update workflows to support cloud storage:

```yaml
# .github/workflows/build_cache_hail.yml
- name: Check for pre-computed BCF
  run: |
    # Try GitHub Releases first
    if wget --spider "https://github.com/${{ github.repository }}/releases/download/${TAG}/${BCF}"; then
      wget "https://github.com/${{ github.repository }}/releases/download/${TAG}/${BCF}"
    # Fall back to cloud storage
    elif wget --spider "https://storage.googleapis.com/vcfstash-bcf-public/v4.1/${BCF}"; then
      wget "https://storage.googleapis.com/vcfstash-bcf-public/v4.1/${BCF}"
    # Generate on-demand
    else
      python scripts/export_gnomad_hail.py --output "${BCF}" ...
    fi
```

## Cost Analysis

### GitHub Releases (Public Repo)

| Item | Cost |
|------|------|
| Storage | FREE |
| Bandwidth | FREE |
| Total | $0/month |

### Google Cloud Storage

| Item | Monthly Cost (100GB) |
|------|---------------------|
| Storage (Standard) | $2.00 |
| Bandwidth (100GB/month) | $1.20 |
| Operations (10K reads) | $0.40 |
| Total | ~$3.60/month |

### AWS S3

| Item | Monthly Cost (100GB) |
|------|---------------------|
| Storage (Standard) | $2.30 |
| Bandwidth (100GB/month) | $9.00 |
| Operations (10K GET) | $0.40 |
| Total | ~$11.70/month |

## Recommended Strategy

### For Public Projects (Recommended)

```
Small BCFs (chr1)           → GitHub Releases (FREE)
Medium BCFs (<1GB)          → GitHub Releases (FREE)
Large BCFs (genome-wide)    → GCS Public Bucket (~$3-5/month)
```

### For Private Projects

```
All BCFs → GCS/S3 Private Bucket (pay per use)
```

### For Organizations with Existing Infrastructure

```
All BCFs → Your existing CDN/object storage
```

## File Organization

### Recommended Directory Structure

```
GitHub Releases:
├── gnomad-bcf-v4.1-af0.10-chr1/
│   ├── gnomad_v4.1_GRCh38_chr1_af0.10.bcf
│   └── gnomad_v4.1_GRCh38_chr1_af0.10.bcf.csi

Cloud Storage (gs://vcfstash-bcf-public/):
├── v4.1/
│   ├── gnomad_v4.1_GRCh38_af0.10.bcf
│   ├── gnomad_v4.1_GRCh38_af0.05.bcf
│   └── gnomad_v4.1_GRCh38_af0.01.bcf
├── v4.0/
│   └── gnomad_v4.0_GRCh38_af0.10.bcf
└── metadata/
    └── file_index.json
```

## Maintenance

### Automated Cleanup

```bash
# Delete old releases (keep last 3 versions)
gh release list | tail -n +4 | awk '{print $1}' | xargs -I {} gh release delete {}

# Delete old cloud storage files
gsutil -m rm gs://vcfstash-bcf-public/v3.*/\*
```

### Monitoring

```bash
# Check storage usage
gsutil du -sh gs://vcfstash-bcf-public

# Check bandwidth
gcloud logging read "resource.type=gcs_bucket" --limit 100
```

## Migration Path

If you need to move from GitHub Releases to Cloud Storage:

```bash
# 1. Download from releases
for tag in $(gh release list | awk '{print $1}'); do
  gh release download "$tag" -D /tmp/bcf-files/
done

# 2. Upload to cloud storage
gsutil -m cp -r /tmp/bcf-files/* gs://vcfstash-bcf-public/

# 3. Update workflows to use new URLs
# (see Workflow Integration section above)
```

## Best Practices

1. **Start with GitHub Releases** - Free and easy for small files
2. **Monitor storage usage** - GitHub has generous limits but check periodically
3. **Use cloud storage for large files** - More cost-effective at scale
4. **Enable CDN caching** - Faster downloads worldwide
5. **Version your BCF files** - Keep multiple gnomAD versions available
6. **Document access URLs** - Make it easy for users to find files

## Security Considerations

### Public Access

- ✅ Safe for gnomAD data (already public)
- ✅ No PHI or sensitive information
- ✅ Read-only access

### Private Access (if needed)

```bash
# GCS with signed URLs (expires after 1 hour)
gsutil signurl -d 1h service-account-key.json \
  gs://vcfstash-bcf-private/file.bcf
```
