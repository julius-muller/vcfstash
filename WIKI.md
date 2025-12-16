# VCFcache Wiki

User-facing documentation for building and using VCFcache. This page is intentionally “walkthrough style”: each section includes enough context to be understood on its own.

---

## Table of Contents
1. What VCFcache does
2. Quick Start
3. Key concepts (Blueprint vs Cache)
4. Finding caches (Zenodo + local)
5. Inspecting a cache (what params does it require?)
6. Building your own cache (end-to-end)
7. Using a cache to annotate samples
8. Configuration reference (`params.yaml` + `annotation.yaml`)
9. Docker usage
10. Troubleshooting
11. CLI reference (all commands & flags)

---

## 1) What VCFcache does

VCFcache accelerates variant annotation by caching annotation results for variants that are seen often (typically from a population dataset such as gnomAD). Instead of re-running the annotation tool for those variants on every sample, VCFcache reuses the cached annotations and only runs the annotation tool for variants that are not found in the cache.

VCFcache is genome-agnostic: you can use it for human, mouse, plants, or any model organism, as long as the blueprint/caches and your samples use compatible reference builds and contig naming (and your annotation pipeline is configured accordingly).

The workflow is:
1. Create a **blueprint** (a normalized, deduplicated list of variants).
2. Annotate the blueprint once to create an annotated **cache**.
3. Annotate samples by doing “cache lookup + annotate leftovers”.

### Design choices (benefits and trade-offs)

- **Tool-agnostic annotation**: the annotation step is defined by an `annotation.yaml` “recipe” that can call VEP or any other annotator you can run as a shell pipeline.
- **Reproducibility by construction**: caches store the immutable recipe (`annotation.yaml`) and a snapshot of the machine-specific parameters used (`params.snapshot.yaml`).
- **Fast runtime, predictable costs**: heavy work happens once during cache creation; per-sample runtime scales with novelty.

Trade-off:
- A cache is only valid for the exact annotation recipe/tooling it was built with (same semantics and compatible tool versions).

---

## 2) Quick Start

This section is for “I want to see it work” with minimal setup.

### Install (pip)

Requires Python 3.11+ and `bcftools >= 1.20`.
```bash
pip install vcfcache
vcfcache --help
```

### Use Docker + a public cache (fastest)

The Docker image includes `bcftools`, so no host install is required.
```bash
docker pull ghcr.io/julius-muller/vcfcache:latest
docker run --rm -v $(pwd):/work ghcr.io/julius-muller/vcfcache:latest \
  annotate \
    -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic \
    --vcf /work/sample.bcf \
    --output /work/out \
    --force
```

If you don’t know what caches exist yet, use:
```bash
vcfcache list caches
```

---

## 3) Key concepts (Blueprint vs Cache)

This section explains what you will see on disk and what each thing is used for.

### Blueprint (the variant set)

A blueprint is the canonical list of variants you want to cache. It removes sample genotypes and standardizes variant representation so the cache lookup is stable and repeatable. You build it once from a population callset, and optionally extend it later.

Location:
`<cache_root>/blueprint/vcfcache.bcf(.csi)`

### Cache (the annotated blueprint)

A cache is an annotated blueprint. A single blueprint can have multiple caches (e.g., VEP “basic”, VEP “everything”, SnpEff, custom pipelines). Each cache is stored under a distinct `annotation_name`.

Location:
`<cache_root>/cache/<annotation_name>/vcfcache_annotated.bcf(.csi)`

### Cache root layout

```
<cache_root>/
├── blueprint/
│   ├── vcfcache.bcf(.csi)
│   └── sources.info
├── cache/
│   └── <annotation_name>/
│       ├── vcfcache_annotated.bcf(.csi)
│       ├── annotation.yaml
│       └── params.snapshot.yaml
└── work/ (temp, safe to delete)
```

Why the extra files in `cache/<annotation_name>/`?
- `annotation.yaml` is the immutable recipe that defines how annotation was produced.
- `params.snapshot.yaml` is the resolved parameter set used when building the cache (paths/tool commands/etc).

---

## 4) Finding caches (Zenodo + local)

This section covers how to discover available caches and where downloads are stored.

### List public items (Zenodo)

```bash
vcfcache list caches
vcfcache list blueprints
```

You can filter the public list:
```bash
vcfcache list caches --genome GRCh37 --source gnomad
```

Zenodo sandbox (useful for testing uploads) is selected with `--debug`:
```bash
vcfcache list caches --debug
```

Shortcuts (same as `vcfcache list ...`):
```bash
vcfcache caches
vcfcache blueprints
```

### Sharing and discovery notes (Zenodo)

- If you maintain a cache for your team or project, **publishing it to Zenodo** is a convenient way to share it: it gives you a DOI, versioning, and a stable download endpoint.
- Public discovery via `vcfcache list caches` works by searching Zenodo records with keywords (`vcfcache` + `cache`/`blueprint` + other keywords like genome/source/tool). When you upload via `vcfcache push`, VCFcache automatically adds the keywords needed for later discovery.

### Where downloads are stored (`VCFCACHE_DIR`)

When VCFcache downloads a public cache/blueprint, it extracts it under a base directory:
- default: `~/.cache/vcfcache`
- override: `VCFCACHE_DIR=/path/to/large/disk`

```bash
export VCFCACHE_DIR=/data/vcfcache
vcfcache list caches
```

### List local items (what you already downloaded/built)

List local caches/blueprints in the default base directory:
```bash
vcfcache caches --local
vcfcache blueprints --local
```

Or point to another base directory:
```bash
vcfcache caches --local --path /data/vcfcache
```

---

## 5) Inspecting a cache (what params does it require?)

Caches depend on their `annotation.yaml` recipe, which can reference `${params.*}` keys that must exist in your `params.yaml` at runtime (or be compatible with the stored snapshot).

Use `--inspect` to see:
- constraints from `annotation.yaml` (e.g., required tool version, required INFO tag)
- which `${params.*}` keys are referenced by the recipe
- which values are present in `params.snapshot.yaml`
- a minimal `params.yaml` template you can copy

Inspect either the cache root or a specific annotation directory:
```bash
vcfcache list --inspect /path/to/cache_root
vcfcache list --inspect /path/to/cache_root/cache/<annotation_name>
```

This is the quickest way to answer “what do I need in my `params.yaml` to use this cache correctly?”

---

## 6) Building your own cache (end-to-end)

This section walks through creating a cache from scratch.

### Prerequisites

- Input population file in BCF/VCF form (ideally BCF) and indexed (e.g., `.csi`).
- `bcftools >= 1.20`.
- An annotation tool/pipeline you can run non-interactively (VEP is a common choice).

### Step 1: Create a blueprint

This creates `<output>/blueprint/vcfcache.bcf` from your population file:
```bash
vcfcache blueprint-init \
  --vcf input.bcf \
  --output ./cache_root \
  -y params.yaml \
  --force
```

If you want multiallelic splitting during blueprint creation:
```bash
vcfcache blueprint-init --vcf input.bcf --output ./cache_root -y params.yaml --normalize --force
```

### Step 2 (optional): Extend the blueprint

Use this when you have additional population VCF/BCF files and want to merge them into the same blueprint:
```bash
vcfcache blueprint-extend --db ./cache_root --vcf extra.bcf -y params.yaml
```

### Step 3: Build an annotated cache from the blueprint

`cache-build` runs the annotation recipe on the blueprint once and writes a new cache under `cache/<name>/`:
```bash
vcfcache cache-build \
  --db ./cache_root \
  --name vep_basic \
  -a annotation.yaml \
  -y params.yaml \
  --force
```

What gets stored (per cache):
- `vcfcache_annotated.bcf(.csi)` — annotated blueprint
- `annotation.yaml` — recipe (copied, treated as immutable)
- `params.snapshot.yaml` — params snapshot used to build the cache

### Step 4 (optional): Share your cache via Zenodo (upload)

If you want to share your blueprint or cache with collaborators (or publish it publicly), you can upload the directory to Zenodo.

Prerequisites:
- Set `ZENODO_TOKEN` (production) or `ZENODO_SANDBOX_TOKEN` (sandbox).
- Start with the sandbox for testing: published Zenodo records cannot be deleted (only versioned).

Upload a cache root directory (auto-detects blueprint vs cache):
```bash
export ZENODO_SANDBOX_TOKEN=...
vcfcache push --cache-dir ./cache_root --test --publish
```

For production:
```bash
export ZENODO_TOKEN=...
vcfcache push --cache-dir ./cache_root --publish
```

After publishing, VCFcache prints a DOI. Others can download it by DOI:
```bash
vcfcache cache-build --doi <DOI>
```

Tip:
- Use `--metadata <file>` to set a nicer Zenodo title/description/creators/keywords.
- Use `--test` for sandbox uploads; use `--debug` for sandbox listing/downloading.

### Alternative: download a pre-built cache by DOI

If you already have a DOI for a pre-built cache, you can download/extract it without setting any environment variables:
```bash
vcfcache cache-build --doi <DOI> -o /data/vcfcache
```

---

## 7) Using a cache to annotate samples

This section covers how to annotate a sample file once you have a cache.

### Annotate via a public alias (auto-download on first use)

```bash
vcfcache annotate \
  -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic \
  --vcf sample.bcf \
  --output outdir \
  --force
```

### Annotate via a local cache directory

`-a` must point to the **annotation directory** (the one containing `annotation.yaml` and `vcfcache_annotated.bcf`):
```bash
vcfcache annotate \
  -a /path/to/cache_root/cache/<annotation_name> \
  --vcf sample.bcf \
  --output outdir \
  --force
```

### Useful annotate options

- Print the stored annotation command (does not run annotation):
  ```bash
  vcfcache annotate -a /path/to/cache_root/cache/<annotation_name> --show-command
  ```
- Run an uncached annotation pass (useful for benchmarking/debugging):
  ```bash
  vcfcache annotate -a /path/to/cache_root/cache/<annotation_name> --vcf sample.bcf --output outdir --uncached --force
  ```
- Convert output to Parquet (for downstream querying):
  ```bash
  vcfcache annotate -a /path/to/cache_root/cache/<annotation_name> --vcf sample.bcf --output outdir --parquet --force
  ```

---

## 8) Configuration reference (`params.yaml` + `annotation.yaml`)

VCFcache deliberately splits configuration into:
- `params.yaml`: environment/tool paths (site-specific, editable per machine)
- `annotation.yaml`: immutable recipe describing the annotation pipeline (cache-specific)

This separation is a key design choice: it lets you distribute caches safely while still allowing users to adapt paths/tool wrappers locally.

### `params.yaml` (tools/resources)

Example:
```yaml
params:
  bcftools_cmd: "bcftools"
  annotation_tool_cmd: "vep"
  tool_version_command: "vep --version"
  temp_dir: "/tmp"
  optional_checks: {}
```

### `annotation.yaml` (immutable recipe)

Example skeleton:
```yaml
annotation_cmd: |
  ${params.bcftools_cmd} view ${INPUT_BCF} | \
  ${params.annotation_tool_cmd} ... | \
  ${params.bcftools_cmd} view -o ${OUTPUT_BCF} -Ob -W
must_contain_info_tag: CSQ
required_tool_version: "115.2"
```

Variables you can use:
- `${INPUT_BCF}`, `${OUTPUT_BCF}`, `${AUXILIARY_DIR}`
- `${params.*}` (values from your `params.yaml`)

### Capturing non-BCF outputs (`${AUXILIARY_DIR}` is required)

Many annotation tools can produce extra outputs besides the annotated BCF (HTML reports, JSON/TSV sidecars, warnings/stderr logs, summary stats, plugin outputs, etc.). VCFcache will **only keep those extra files** if your `annotation_cmd` writes them into `${AUXILIARY_DIR}`.

Important implications:
- Where it ends up: during `cache-build` it is stored inside the cache annotation directory (next to `vcfcache_annotated.bcf`); during `annotate` it is stored in your chosen `--output` directory as `./auxiliary/`.
- If you run the tool in Docker/Apptainer, ensure `${AUXILIARY_DIR}` is a path that exists **inside** the container (VCFcache creates it on the host and substitutes the path into the command). Your wrapper/command must actually write/copy outputs there.
- If the tool writes to the current working directory or to some internal temp directory, those files may be lost after the workflow finishes.
- If nothing is written to `${AUXILIARY_DIR}`, VCFcache removes the empty directory automatically.

Example (VEP-style outputs captured into `${AUXILIARY_DIR}`):
```yaml
annotation_cmd: |
  ${params.bcftools_cmd} view ${INPUT_BCF} -Ou | \
  ${params.annotation_tool_cmd} \
    --format vcf --vcf --offline \
    --stats_file ${AUXILIARY_DIR}/vep_stats.html \
    --warning_file ${AUXILIARY_DIR}/vep_warnings.txt \
    2> ${AUXILIARY_DIR}/vep_stderr.txt | \
  ${params.bcftools_cmd} view -o ${OUTPUT_BCF} -Ob -W
```

Practical tip:
- Use `vcfcache list --inspect <cache>` to see the exact `${params.*}` keys your cache requires and generate a minimal template.

### bcftools requirement / override

VCFcache requires `bcftools >= 1.20`. If your system `bcftools` is too old, point VCFcache to a newer binary:
```bash
export VCFCACHE_BCFTOOLS=/path/to/bcftools-1.22
```

---

## 9) Docker usage

Docker is a convenient way to avoid installing `bcftools` on the host and to run VCFcache consistently across machines.

### Run the published image

```bash
docker pull ghcr.io/julius-muller/vcfcache:latest
docker run --rm ghcr.io/julius-muller/vcfcache:latest --help
```

### Build locally

```bash
docker build --target final -f docker/Dockerfile.vcfcache -t vcfcache:local .
```

---

## 10) Troubleshooting

- **bcftools not found or too old**: install `bcftools >= 1.20`, or set `VCFCACHE_BCFTOOLS=...`.
- **Nothing shows up in `--local` listing**: confirm `VCFCACHE_DIR` and directory layout (`<base>/caches/*` and `<base>/blueprints/*`).
- **Alias not found on Zenodo**: verify with `vcfcache list caches` and ensure you used the right Zenodo environment (`--debug` uses sandbox).
- **Downloads too large for $HOME**: set `VCFCACHE_DIR` to a large disk.

---

## 11) CLI reference (all commands & flags)

This section is a compact reference for the CLI. It complements the walkthrough above.

### `vcfcache blueprint-init`

Create a blueprint from a local VCF/BCF, or download a blueprint from Zenodo.

- `-i/--vcf`: input VCF/BCF (local blueprint creation).
- `--doi`: download blueprint by DOI (mutually exclusive with `--vcf`).
- `-o/--output`: output directory for the blueprint (default: `./cache`).
- `-y/--yaml`: params YAML used during local blueprint creation (tool paths, temp dir).
- `-t/--threads`: threads for `bcftools` operations (default: 1).
- `-n/--normalize`: split multiallelic variants during blueprint creation.
- `-f/--force`: overwrite output directory if it exists.
- `--debug`: keep intermediate work dirs and use Zenodo sandbox for list/download operations.

### `vcfcache blueprint-extend`

Extend an existing blueprint with additional variants.

- `--db`: blueprint directory (cache root containing `blueprint/`).
- `--vcf`: additional VCF/BCF to merge into the blueprint.
- `-y/--yaml`: params YAML used for local blueprint operations.
- `-t/--threads`: threads for `bcftools` operations (default: 1).
- `-n/--normalize`: split multiallelic variants during extension.
- `-f/--force`: overwrite temporary outputs if present.
- `--debug`: keep intermediate work dirs.

### `vcfcache cache-build`

Two modes:
- Build a cache from a blueprint (`--db` or `--doi` pointing to a blueprint): requires `-a/--anno-config` + `-n/--name`.
- Download a pre-built cache (`--doi` pointing to a cache): forbids `-a` and ignores `-n`.

Flags:
- `-d/--db`: local blueprint directory.
- `--doi`: DOI for blueprint or cache.
- `-o/--output`: base directory for DOI downloads (overrides `VCFCACHE_DIR` for this command).
- `-n/--name`: cache name when building from a blueprint.
- `-a/--anno-config`: annotation recipe YAML (required when building from a blueprint).
- `-y/--params`: params YAML (copied to `params.snapshot.yaml`).
- `-t/--threads`: threads for `bcftools` operations (default: 1).
- `-f/--force`: overwrite existing cache outputs.
- `--debug`: keep intermediate work dirs and use Zenodo sandbox for DOI operations.

### `vcfcache annotate`

Annotate a sample VCF/BCF using a cache.

- `-a/--annotation_db`: cache annotation directory (local path) or cache alias (Zenodo).
- `-i/--vcf`: input sample VCF/BCF.
- `-o/--output`: output directory.
- `-y/--yaml`: params YAML for runtime (copied to workflow `params.snapshot.yaml`).
- `-t/--threads`: threads for `bcftools` (default: 1).
- `-f/--force`: overwrite outputs.
- `--uncached`: run annotation without using the cache (debug/benchmark).
- `--parquet`: also convert output to Parquet for downstream querying.
- `--show-command`: print stored annotation command and exit.
- `--list`: list annotation caches under `-a` (when `-a` points to a directory of caches).
- `--debug`: keep intermediate files and use sandbox for alias resolution/downloads.
- `-v/--verbose`: increase logging (`-v` INFO, `-vv` DEBUG).

### `vcfcache list` / `vcfcache caches` / `vcfcache blueprints`

Discover public items on Zenodo or list local items.

- `caches|blueprints`: choose which type to list (default is `blueprints` for `vcfcache list`).
- `--genome`: filter public Zenodo results by genome keyword (e.g. `GRCh37`, `GRCh38`).
- `--source`: filter public Zenodo results by source keyword (e.g. `gnomad`).
- `--debug`: query Zenodo sandbox instead of production.
- `--local`: list local items instead of querying Zenodo.
- `--path`: base directory to search for local listing (defaults to `VCFCACHE_DIR` or `~/.cache/vcfcache`).
- `--inspect`: inspect a local cache/blueprint directory and print required `${params.*}` keys.

### `vcfcache push`

Upload a local blueprint or cache directory to Zenodo (great for sharing caches).

- `--cache-dir`: path to a cache root directory (auto-detects blueprint vs cache).
- `--dest`: upload destination (currently `zenodo`).
- `--test`: upload to Zenodo sandbox (uses `ZENODO_SANDBOX_TOKEN`).
- `--metadata`: optional YAML/JSON file to set Zenodo metadata (title/description/creators/keywords).
- `--publish`: publish the record immediately after upload.

Environment variables:
- `ZENODO_TOKEN`: production token.
- `ZENODO_SANDBOX_TOKEN`: sandbox token.
- `ZENODO_SANDBOX=1`: affects some operations; for VCFcache CLI, `--debug` and `--test` are the main toggles.

### `vcfcache demo`

Convenience command for smoke testing and benchmarking.

- `--smoke-test`: runs a comprehensive workflow demo.
- `-q/--quiet`: suppress detailed output.
- Benchmark mode: provide `-a <cache>` and `--vcf <file>` (and `-y <params>`) to compare cached vs uncached behavior.

---

## 12) Variant Preservation Feature

### Overview

**VCFcache preserves ALL user variants** in cached mode, even if the annotation tool drops some during annotation. This is a deliberate design choice to ensure users don't lose their data.

### Background: The Annotation Tool Dropping Problem

Many annotation tools (like VEP) silently drop variants from non-standard contigs or variants they cannot process. For example:

```bash
# Input VCF: 1000 variants (including variants on chrUn, chrEBV, etc.)
vep --input input.vcf --output output.vcf
# Output VCF: 950 variants (50 variants on non-standard contigs were dropped!)
```

**This is frustrating because:**
- Users lose variants without clear warning
- Dropped variants may still be scientifically interesting
- No way to recover the dropped variants later

### VCFcache's Solution: Variant Preservation

When using vcfcache in **cached mode**, all input variants are preserved:

```bash
# Uncached annotation (standard VEP behavior)
vcfcache annotate --uncached -a cache/ --vcf input.vcf --output uncached/
# Output: 950 annotated variants (50 dropped by VEP)

# Cached annotation (vcfcache preserves all variants)
vcfcache annotate -a cache/ --vcf input.vcf --output cached/
# Output: 1000 variants (950 annotated + 50 preserved without annotation)
```

### How It Works

**4-Step Cached Workflow:**

1. **Normalize input**: Split multiallelic variants
2. **Filter missing**: Extract variants not in cache
3. **Annotate missing**: Run annotation tool (may drop some variants)
4. **Merge annotations**: Merge annotated variants back into normalized input

**Key behavior:**
- If annotation tool drops variants in step 3, those variants remain in the final output (without annotation)
- A warning is logged showing how many variants were dropped:
  ```
  WARNING: Annotation tool dropped 50 variants from input.
  These variants are preserved in cached output (without CSQ annotation)
  but would be missing in uncached output.
  This is a FEATURE - vcfcache preserves all your variants!
  ```

### Comparing Cached vs Uncached Outputs

Because cached preserves all variants while uncached doesn't, **total variant counts may differ**:

```bash
# Cached output
bcftools view cached/sample_vc.bcf | wc -l
# 1000 variants

# Uncached output
bcftools view uncached/sample_vc.bcf | wc -l
# 950 variants
```

**To compare only annotated variants** (the fair comparison):

```bash
# Filter both sides to annotated variants only
bcftools view -i 'INFO/CSQ!=""' cached/sample_vc.bcf > cached_annotated.vcf
bcftools view -i 'INFO/CSQ!=""' uncached/sample_vc.bcf > uncached_annotated.vcf

# Now compare
bcftools view -H cached_annotated.vcf | md5sum
bcftools view -H uncached_annotated.vcf | md5sum
# MD5 hashes should be identical (or very similar - minor order differences acceptable)
```

### Identifying Preserved Variants

To find variants that were preserved without annotation:

```bash
# Variants without CSQ annotation
bcftools view -i 'INFO/CSQ=""' cached/sample_vc.bcf

# Count unannotated variants
bcftools view -i 'INFO/CSQ=""' cached/sample_vc.bcf | grep -v "^#" | wc -l
```

### Use Cases

**When variant preservation is valuable:**
- Research on non-standard contigs (chrUn, alternative haplotypes, etc.)
- Quality control - knowing which variants couldn't be annotated
- Post-processing - manually annotate dropped variants later
- Reproducibility - keeping complete variant sets

**When you want only annotated variants:**
Simply filter the output to annotated variants only:

```bash
bcftools view -i 'INFO/CSQ!=""' cached/sample_vc.bcf -o annotated_only.bcf -Ob
```

### Testing and Validation

VCFcache includes comprehensive tests to ensure:
1. Cached annotation preserves all input variants
2. Annotated variant sets are identical between cached and uncached (when filtered to annotated)
3. Warnings are issued when annotation tools drop variants
4. MD5 comparison (when filtering to annotated variants) shows consistency

See `tests/test_contig_mismatches.py` for detailed test implementation.

### Configuration

This behavior is **always enabled** in cached mode and cannot be disabled. It is a core feature of vcfcache's design philosophy: **preserve user data whenever possible**.

If you need behavior identical to uncached mode, simply filter the output:
```bash
vcfcache annotate -a cache/ --vcf input.vcf --output output/
bcftools view -i 'INFO/CSQ!=""' output/sample_vc.bcf -o filtered.bcf -Ob
```

### Performance

Variant preservation has **minimal performance overhead** because:
- No additional filtering or splitting operations required
- Standard merge operation handles all variants
- Warning calculation is a simple count comparison

The overhead is negligible compared to the annotation step itself.
