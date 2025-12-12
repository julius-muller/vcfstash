# VCFcache Wiki

Practical reference for building, distributing, and consuming VCFcache blueprints and caches.

---

## Table of Contents
1. Quick Start
2. Aliases & Manifest
3. Cache Layout
4. Building Caches
5. Using Public Caches (Zenodo)
6. Configuration (YAML)
7. Docker Build/Testing
8. Troubleshooting

---

## 1) Quick Start

From source (dev):
```bash
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
vcfcache --help
python -m pytest tests -q
```

Docker + public cache:
```bash
docker pull ghcr.io/julius-muller/vcfcache:v0.3.0
docker run --rm -v $(pwd):/work ghcr.io/julius-muller/vcfcache:v0.3.0 \
  annotate \
    -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic \
    --vcf /work/tests/data/nodata/sample4.bcf \
    --output /work/out \
    --force
```

List available public items from Zenodo:
```bash
vcfcache list caches
vcfcache list blueprints
```

---

## 2) Aliases & Discovery

Alias schema:
- Blueprint: `bp-<genome>-<source>-<release>-<filt>`
- Cache: `cache-<genome>-<source>-<release>-<filt>-<tool>-<tool_version>-<preset>`

Helper commands:
- Discover: `vcfcache list [blueprints|caches]`
- Resolve/download: `vcfcache annotate -a <alias>` (auto downloads from Zenodo to `~/.cache/vcfcache/caches/`)
- Upload: `vcfcache push --cache-dir <dir>` (requires `ZENODO_TOKEN`; use `--test` + `ZENODO_SANDBOX_TOKEN` for sandbox)

---

## 3) Cache Layout (on disk)

```
<cache_root>/
├── blueprint/
│   ├── vcfcache.bcf(.csi)
│   └── sources.info
├── cache/
│   └── <alias or name>/
│       ├── vcfcache_annotated.bcf(.csi)
│       ├── annotation.yaml
│       └── params.snapshot.yaml
└── work/ (temp, safe to delete)
```

Pure Python implementation. Docker image includes bcftools; for pip install, bcftools >= 1.20 must be installed separately.

---

## 4) Building Caches

Blueprint:
```bash
vcfcache blueprint-init \
  --vcf input.bcf \
  --output ./cache_root \
  -y params.yaml \
  --force
```

Extend blueprint:
```bash
vcfcache blueprint-extend \
  --db ./cache_root \
  --vcf extra.bcf \
  -y params.yaml
```

Cache-build (annotate blueprint):
```bash
vcfcache cache-build \
  --name vep_custom \
  --db ./cache_root \
  -a annotation.yaml \
  -y params.yaml \
  --force
```

---

## 5) Using Public Caches (Zenodo)

Annotate with an alias (auto download/extract):
```bash
vcfcache annotate \
  -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic \
  --vcf sample.bcf \
  --output outdir \
  --force
```

Show recorded annotation command:
```bash
vcfcache annotate -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic --show-command
```

---

## 6) Configuration (YAML)

`params.yaml` (tools/resources):
```yaml
bcftools_cmd: "bcftools"
annotation_tool_cmd: "vep"
tool_version_command: "vep --version"
temp_dir: "/tmp"
optional_checks: {}
```

`annotation.yaml` (immutable recipe):
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

Special vars: `${INPUT_BCF}`, `${OUTPUT_BCF}`, `${AUXILIARY_DIR}`, `${params.*}`.

**bcftools Override:**

vcfcache requires bcftools >= 1.20. Override system bcftools with:
```bash
export VCFCACHE_BCFTOOLS=/path/to/bcftools-1.22
vcfcache annotate ...
```

Or in params.yaml:
```yaml
bcftools_cmd: "/opt/bcftools/bin/bcftools"
```

The environment variable affects vcfcache's internal bcftools checks, while `bcftools_cmd` in params.yaml is used within annotation commands.

---

## 7) Docker Build/Testing

Multi-stage Dockerfile: `docker/Dockerfile.vcfcache`
- `--target test` installs dev deps and can run pytest:
  ```bash
  docker build --target test -f docker/Dockerfile.vcfcache -t vcfcache:test .
  docker run --rm vcfcache:test /opt/venv/bin/python -m pytest tests -q
  ```
- `--target final` is the lean runtime image:
  ```bash
  docker build --target final -f docker/Dockerfile.vcfcache -t ghcr.io/julius-muller/vcfcache:v0.3.0 .
  ```

---

## 8) Troubleshooting

- **bcftools not found or too old** (pip install): Install bcftools >= 1.20, or set `VCFCACHE_BCFTOOLS=/path/to/newer/bcftools` to override system version. Docker image already includes bcftools.
- **Alias not found on Zenodo**: Ensure the Zenodo record has keywords `vcfcache`, `cache`/`blueprint`, and the alias; verify with `vcfcache list`.
- **Zenodo upload**: Export `ZENODO_TOKEN`; set `ZENODO_SANDBOX=1` to target the sandbox API.
- **Large Docker images**: Use `--target final` for runtime; `--target test` is only for CI/dev. 
