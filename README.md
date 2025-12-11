# VCFcache – Cache once, annotate fast

Cache common variants once, reuse them for every sample. VCFcache builds a normalized blueprint, annotates it once, and reuses those results so only rare/novel variants are annotated at runtime.

---

## What’s new
- Single lean Docker image: `ghcr.io/julius-muller/vcfcache:<version>`
- Caches/blueprints live on Zenodo and are resolved by alias via `public_caches.yaml`
- Alias schema:
  - Blueprint: `bp-<genome>-<source>-<release>-<filt>`
  - Cache: `cache-<genome>-<source>-<release>-<filt>-<tool>-<tool_version>-<preset>`
- Manifest entries carry `version` (vcfcache version used), no image tags
- Bundled bcftools 1.22; no Nextflow/JVM required

---

## Quick start (pip install)

Requirements: Python >= 3.11, bcftools >= 1.20

```bash
# Install vcfcache
pip install vcfcache

# Run comprehensive demo (tests all 4 commands)
vcfcache demo

# Use it
vcfcache --help
```

**Note:** bcftools >= 1.20 is a runtime dependency and must be installed separately:
- Ubuntu/Debian: `sudo apt-get install bcftools`
- macOS: `brew install bcftools`
- Conda: `conda install -c bioconda bcftools`
- Or set `VCFCACHE_BCFTOOLS=/path/to/bcftools` to use a specific version

The demo command runs a complete workflow testing:
1. `blueprint-init` - Create cache from variants
2. `blueprint-extend` - Add more variants
3. `cache-build` - Annotate the cache
4. `annotate` - Use cache to annotate sample VCF

---

## Quick start (from source)

Requirements: Python 3.13+, uv, bcftools >= 1.20

```bash
git clone https://github.com/julius-muller/vcfcache.git
cd vcfcache
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
vcfcache --help
python -m pytest tests -q   # optional
```

---

## Quick start (Docker + public cache)

The CLI can fetch caches from Zenodo by alias using the bundled `public_caches.yaml`.

```bash
docker pull ghcr.io/julius-muller/vcfcache:v0.3.0

docker run --rm -v $(pwd):/work ghcr.io/julius-muller/vcfcache:v0.3.0 \
  annotate \
    -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic \
    --vcf /work/tests/data/nodata/sample4.bcf \
    --output /work/out \
    --force
```

List public caches:
```bash
vcfcache list --public-caches
```

Show the recorded annotation command without running it:
```bash
vcfcache annotate -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic --show-command
```

---

## Build your own cache

1) Normalize/deduplicate into a blueprint:
```bash
vcfcache blueprint-init \
  --vcf gnomad_subset.bcf \
  --output ./cache \
  -y params.yaml \
  --force
```

2) Annotate the blueprint (cache-build) with your recipe:
```bash
vcfcache cache-build \
  --name vep_custom \
  --db ./cache \
  -a annotation.yaml \
  -y params.yaml \
  --force
```

3) Use the cache on samples:
```bash
vcfcache annotate \
  -a ./cache/cache/vep_custom \
  --vcf sample.bcf \
  --output ./results \
  --force
```

---

## Manifest & aliases

- Manifest: `public_caches.yaml`
- Fields (minimum): `alias`, `type` (blueprint|cache), `version`, `genome`, `source`, `release`, `filt`, `doi`, `updated_at`, `md5`, plus `tool/tool_version/preset` for caches.
- Aliases:
  - Blueprint: `bp-hg38-gnomad-4.1joint-AF0100`
  - Cache: `cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic`
- `vcfcache annotate -a <alias>` auto-downloads from Zenodo and extracts to `~/.cache/vcfcache/caches/<alias>`.
- `vcfcache pull --doi <doi> --dest <dir>` for manual download/extract.
- `vcfcache push --cache-dir <dir>` uploads to Zenodo (requires `ZENODO_TOKEN`; uses real Zenodo by default).

---

## Configuration (YAML)

Two files, stored with the cache for reproducibility:

1) `params.yaml` — resource paths/tools:
```yaml
bcftools_cmd: "bcftools"         # bundled bcftools by default
annotation_tool_cmd: "vep"
tool_version_command: "vep --version"
temp_dir: "/tmp"
optional_checks: {}
```

2) `annotation.yaml` — immutable annotation recipe:
```yaml
annotation_cmd: |
  ${params.bcftools_cmd} view ${INPUT_BCF} | \
  ${params.annotation_tool_cmd} --offline --cache \
    --dir_cache ${params.vep_cache} --format vcf --vcf \
    -i STDIN -o STDOUT | \
  ${params.bcftools_cmd} view -o ${OUTPUT_BCF} -Ob -W
must_contain_info_tag: CSQ
required_tool_version: "115.2"
optional_checks:
  vep_cache_version: "115"
  genome_build: "GRCh38"
```

Special variables: `${INPUT_BCF}`, `${OUTPUT_BCF}`, `${AUXILIARY_DIR}`, and `${params.*}`.

---

## bcftools Configuration

vcfcache requires bcftools >= 1.20. The system will automatically search for bcftools in your PATH, but you can override this in two ways:

**1. Environment variable (recommended for non-standard installations):**
```bash
# Override system bcftools
export VCFCACHE_BCFTOOLS=/path/to/bcftools-1.22/bin/bcftools
vcfcache annotate ...

# Or inline for a single command
VCFCACHE_BCFTOOLS=/opt/bcftools/bin/bcftools vcfcache demo --smoke-test
```

**2. Via params.yaml (for annotation workflows):**
```yaml
bcftools_cmd: "/usr/local/bin/bcftools-1.22"
```

The environment variable takes precedence over PATH search, but `bcftools_cmd` in params.yaml is used within annotation commands. Use the environment variable when you need to override the system bcftools globally, and use `bcftools_cmd` in params.yaml for annotation-specific bcftools configuration.

---

## Testing

- From source: `python -m pytest tests -q`
- In Docker (dev stage):  
  `docker build --target test -f docker/Dockerfile.vcfcache -t vcfcache:test .`  
  `docker run --rm vcfcache:test /opt/venv/bin/python -m pytest tests -q`

---

## Docker

- Lean runtime: `ghcr.io/julius-muller/vcfcache:v0.3.0` (bundled bcftools 1.22)
- Multi-stage build: `docker/Dockerfile.vcfcache`
  - `--target test` installs dev deps and can run pytest
  - `--target final` is the shipped image (small)

---

## Notes

- No Nextflow/JVM needed.
- Caches are portable tarballs; manifest `version` is informational to track which vcfcache built them.
- bcftools is bundled and used by default; override via `bcftools_cmd` if you need a different binary.
