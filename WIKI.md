# VCFcache Wiki

User-facing documentation for building and using VCFcache. This page is intentionally “walkthrough style”: each section includes enough context to be understood on its own.

---

## Table of Contents
1. What VCFcache does
2. Performance model (runtime efficiency)
3. Quick Start
4. Key concepts (Blueprint vs Cache)
5. Finding caches (Zenodo + local)
6. Inspecting a cache (what params does it require?)
7. Building your own cache (end-to-end)
8. Using a cache to annotate samples
9. Configuration reference (`params.yaml` + `annotation.yaml`)
10. Docker usage
11. Troubleshooting
12. CLI reference (all commands & flags)

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

## 2) Performance model (runtime efficiency)

VCFcache dramatically reduces per-sample annotation time by caching results for common variants and reusing them instead of re-computing each time. This section explains the runtime model and expected speed gains.

### How cache lookups work

Cache lookups are extremely fast – VCFcache uses bcftools/htslib to query a pre-indexed binary BCF cache file, so retrieving an annotation for a given variant is effectively a **constant-time operation**. Thanks to the BCF index, bcftools provides near-instant random access to any variant's data. Using binary BCF format instead of text VCF further reduces overhead (2-3× faster for queries).

In short: **looking up a cached annotation incurs minimal, constant cost regardless of cache size**, whereas annotating a variant from scratch is orders of magnitude slower.

### Why conventional annotation is slow

Conventional annotation pipelines (e.g., Ensembl VEP with plugins like SpliceAI, dbNSFP, or CADD) are **I/O bound, not compute bound**. They don't perform heavy calculations - they're essentially lookup engines that retrieve pre-computed scores and predictions from large external data files.

The bottleneck is **data access**: For each variant, the pipeline typically queries multiple large external files (tabix-indexed VCFs, BED files, BigWig tracks, database files, etc.). Each query involves:
- Seeking to the genomic position in a large file
- Decompressing data blocks (for bgzipped files)
- Parsing the relevant records
- Often querying wide genomic intervals or scanning multiple resources

Even though the annotation tool itself is fast, **pulling data from many large files is slow**. The algorithms behind tools like tabix and BigWig are optimized for one-off lookups (as in genome browsers) but not for processing thousands of sequential queries in batch. Each per-variant query is handled as an independent search, repeating the same I/O operations over and over.

**Example**: If a variant requires data from 5 different external files (e.g., VEP cache, dbNSFP, CADD, SpliceAI, conservation scores), and each file access takes 2 seconds due to file size and I/O, that's 10 seconds per variant - but the actual computation is negligible.

### Runtime model variables

Let us define a few variables to formalize the runtime for a given sample:

- **N** – number of variants in the sample
- **f** – fraction of variants that hit the cache (cache hit rate, typically 0.70-0.90)
- **t_ann** – average time to fully annotate one variant without caching (per-variant annotation time)
- **t_lookup** – time to fetch one variant's annotation from cache (very small: ~0.001-0.01s)
- **t_overhead** – one-time fixed overhead for preparing the job and intersecting with cache

### Runtime comparison

**Baseline (no cache)**: All variants annotated from scratch
```
T_baseline = N × t_ann
```

**With VCFcache**: Cached variants handled much faster
```
T_cached = t_overhead + (f × N × t_lookup) + ((1-f) × N × t_ann)
```

Because cache lookup is extremely fast (t_lookup ≪ t_ann), this simplifies to:
```
T_cached ≈ t_overhead + (1-f) × N × t_ann
```

For negligible overhead:
```
T_cached ≈ (1-f) × N × t_ann
```

### Expected speed-up

The speed-up factor is:
```
Speed-up ≈ T_baseline / T_cached ≈ 1 / (1-f)
```

**Practical examples**:
- **f = 0.80** (80% cache hits) → **~5× faster**
- **f = 0.90** (90% cache hits) → **~10× faster**
- **f = 0.95** (95% cache hits) → **~20× faster**

In practice, using a large population cache (like one built from gnomAD) often yields **70-90% cache hit rates** for typical sample VCFs, resulting in **5-10× speed-ups**.

### Impact of pipeline complexity

The benefit of caching **grows with pipeline complexity**. As the annotation pipeline incorporates more resource-intensive analyses (e.g., SpliceAI splice predictions or deep-learning scores), the average annotation time t_ann per variant increases, making the baseline runtime grow linearly.

However, **cache lookup time (t_lookup) and overhead (t_overhead) remain flat** regardless of pipeline complexity. No matter how many annotations are in the pipeline, fetching a cached result is still just a quick indexed lookup.

Therefore: **the more complex (and slow) the annotation pipeline, the more advantageous caching becomes**.

### Scalability: cache size independence

Cache lookup performance is **independent of cache size**. A cache with ten million variants can be queried just as quickly as a cache with ten thousand variants. Each lookup is a direct index access (like a dictionary lookup by position) with negligible cost.

This ensures vcfcache remains scalable – you can build very comprehensive caches covering broad population datasets without worrying about slowing down per-variant lookup performance. **Each sample's runtime scales primarily with the novelty of its variants, not with the size of the cache**.

### Contig naming considerations

VCFcache requires **matching contig naming** between your cache and input samples. For example, if your cache uses `chr1`, `chr2`, etc., your samples must also use `chr1`, `chr2`, etc. Mixing naming conventions (e.g., cache with `chr1` and sample with `1`) will result in an error. At annotation start, vcfcache reports the contig overlap and fails fast if there is no overlap.

### Genome build compatibility

The genome build must be explicitly set in both `params.yaml` and `annotation.yaml` (e.g., `GRCh38`, `GRCh37`). VCFcache validates that these values match and logs the genome build during workflow startup.

### When caching provides less benefit

For **very small datasets** or **very simple annotation pipelines**, the fixed overhead (t_overhead) may dominate the runtime, reducing or eliminating speed gains. Specifically:

- **Small N** (few variants): If your sample has only dozens or hundreds of variants, the overhead of setting up the caching workflow may be comparable to just annotating them directly
- **Simple pipelines with minimal data access** (small t_ann): If your annotation pipeline only adds a few basic fields without accessing large external data sources (e.g., VEP without any plugins, or simple dbSNP ID lookups), there's less data access cost to save. The benefit of caching grows with the number and size of external data sources your pipeline accesses (BigWig files, large tabix-indexed databases, VEP plugins like SpliceAI/dbNSFP/CADD, etc.)
- **Testing with tiny datasets**: If you're testing vcfcache with a small test VCF, you may not see significant speed-ups and might even observe slower performance due to overhead

**Key insight**: Annotation pipelines are essentially lookup engines - they don't do heavy computation, they access data from large external files. The per-variant time (t_ann) is dominated by I/O cost: seeking through BigWig files, querying tabix-indexed VCFs, loading plugin data, etc. Even if it takes 10 seconds to pull all the required data for a single variant from various large files, fetching that same pre-computed result from the cache takes milliseconds - **a 100× speed-up per cached variant**. The more external data sources your pipeline queries, the larger t_ann becomes, and the more dramatic the caching benefit.

**Recommendation**: VCFcache is designed for production annotation workflows where you have many samples (hundreds to thousands) with typical variant counts (thousands to millions) and annotation pipelines that access multiple external data sources. The benefits compound over many samples - the one-time cache build cost is amortized across all future samples.

### Summary

VCFcache shifts the runtime model to **"annotate once, reuse often"**, delivering:
- **Predictable speed-ups**: 2-10× for typical samples (60-90% cache hits)
- **Greater gains for complex pipelines**: More expensive annotations = higher relative benefit
- **Constant-time lookups**: Cache size doesn't affect query performance
- **Scalable design**: Runtime determined by variant novelty, not cache size
- **Best for production workflows**: Benefits compound when annotating many samples

---

## 3) Quick Start

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
    --output /work/sample_vc.bcf \
    --stats-dir /work/out \
    --force
```

If you don't know what caches exist yet, use:
```bash
vcfcache list caches
```

Note: Run `vcfcache demo` without arguments to see available demo modes.

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

## 5) Finding caches (Zenodo + local)

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
vcfcache list caches --local
vcfcache list blueprints --local
```

Or point to a custom directory:
```bash
vcfcache list caches --local /data/vcfcache
```

---

## 6) Inspecting a cache (what params does it require?)

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

## 7) Building your own cache (end-to-end)

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

## 8) Using a cache to annotate samples

This section covers how to annotate a sample file once you have a cache.

### Annotate via a public alias (auto-download on first use)

```bash
vcfcache annotate \
  -a cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic \
  --vcf sample.bcf \
  --output sample_vc.bcf \
  --stats-dir outdir \
  --force
```

### Annotate via a local cache directory

`-a` must point to the **annotation directory** (the one containing `annotation.yaml` and `vcfcache_annotated.bcf`).
Start with `--requirements`: it explains what tool and parameters the cache expects, and validates tool versions.
```bash
vcfcache annotate -a /path/to/cache_root/cache/<annotation_name> --requirements
vcfcache annotate \
  -a /path/to/cache_root/cache/<annotation_name> \
  --vcf sample.bcf \
  --output sample_vc.bcf \
  --stats-dir outdir \
  --force
```

### Useful annotate options

- Print cache requirements and validate your environment (does not run annotation).  
  If `-y/--yaml` is provided, it checks that file; otherwise it uses the cache’s `params.snapshot.yaml`:
  ```bash
  vcfcache annotate -a /path/to/cache_root/cache/<annotation_name> --requirements
  vcfcache annotate -a /path/to/cache_root/cache/<annotation_name> --requirements -y params.yaml
  ```
- Run an uncached annotation pass (useful for benchmarking/debugging):
  ```bash
  vcfcache annotate -a /path/to/cache_root/cache/<annotation_name> --vcf sample.bcf --output sample_vc.bcf --stats-dir outdir --uncached --force
  ```
- Convert output to Parquet (for downstream querying):
  ```bash
  vcfcache annotate -a /path/to/cache_root/cache/<annotation_name> --vcf sample.bcf --output sample_vc.bcf --stats-dir outdir --parquet --force
  ```
- Preserve variants without annotation in output (by default, vcfcache mirrors annotation tool behavior):
  ```bash
  vcfcache annotate -a /path/to/cache_root/cache/<annotation_name> --vcf sample.bcf --output sample_vc.bcf --stats-dir outdir --preserve-unannotated --force
  ```
- Skip splitting multiallelic variants for small speedup (~6% of runtime):
  ```bash
  vcfcache annotate -a /path/to/cache_root/cache/<annotation_name> --vcf sample.bcf --output sample_vc.bcf --stats-dir outdir --skip-split-multiallelic --force
  ```
  **WARNING**: Use `--skip-split-multiallelic` ONLY if you are certain your input VCF has no multiallelic variants (no commas in ALT field). If multiallelic variants are present, this will cause format inconsistencies between cached and uncached outputs. By default, vcfcache always splits multiallelic variants using `bcftools norm -m-` to ensure consistent output format.

### Input preprocessing

VCFcache automatically preprocesses input VCF/BCF files before annotation:

1. **Multiallelic splitting** (unless `--skip-split-multiallelic` is used): Splits multiallelic variants into biallelic records using `bcftools norm -m-`
2. **Spanning deletion removal**: Removes variants with `ALT=*` (spanning deletion alleles), which are VCF placeholder notation that annotation tools like VEP cannot process

These preprocessing steps ensure:
- Consistent input format for annotation tools
- Identical outputs between cached and uncached modes
- No unannotated variants in the output (annotation tools skip `ALT=*` variants)

---

## 9) Configuration reference (`params.yaml` + `annotation.yaml`)

VCFcache deliberately splits configuration into:
- `params.yaml`: environment/tool paths (site-specific, editable per machine)
- `annotation.yaml`: immutable recipe describing the annotation pipeline (cache-specific)

This separation is a key design choice: it lets you distribute caches safely while still allowing users to adapt paths/tool wrappers locally.

### `params.yaml` (tools/resources)

Example:
```yaml
bcftools_cmd: "bcftools"
annotation_tool_cmd: "vep"
tool_version_command: "vep --version"
temp_dir: "/tmp"
threads: 1
genome_build: "GRCh38"
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
genome_build: "GRCh38"
```

Variables you can use:
- `${INPUT_BCF}`, `${OUTPUT_BCF}`, `${AUXILIARY_DIR}`
- `${params.*}` (values from your `params.yaml`)

### Capturing non-BCF outputs (`${AUXILIARY_DIR}` is required)

Many annotation tools can produce extra outputs besides the annotated BCF (HTML reports, JSON/TSV sidecars, warnings/stderr logs, summary stats, plugin outputs, etc.). VCFcache will **only keep those extra files** if your `annotation_cmd` writes them into `${AUXILIARY_DIR}`.

Important implications:
- Where it ends up: during `cache-build` it is stored inside the cache annotation directory (next to `vcfcache_annotated.bcf`); during `annotate` it is stored under your stats directory at `./<output_file>_vcstats/auxiliary/`.
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

## 10) Docker usage

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

## 11) Troubleshooting

- **bcftools not found or too old**: install `bcftools >= 1.20`, or set `VCFCACHE_BCFTOOLS=...`.
- **Nothing shows up in `--local` listing**: confirm `VCFCACHE_DIR` and directory layout (`<base>/caches/*` and `<base>/blueprints/*`).
- **Alias not found on Zenodo**: verify with `vcfcache list caches` and ensure you used the right Zenodo environment (`--debug` uses sandbox).
- **Downloads too large for $HOME**: set `VCFCACHE_DIR` to a large disk.
- **No speed-up with small test datasets**: See [Performance model - When caching provides less benefit](#when-caching-provides-less-benefit). VCFcache is designed for production workflows with many samples.

### Known issues with annotation tools

**VEP non-deterministic output (affects VEP ≥ 113, not a VCFcache issue)**

In some VEP versions, annotation output may vary non-deterministically depending on the specific set of variants in the input, even when annotating the same individual variants. This affects both cached and uncached VCFcache runs identically, and also affects standard VEP usage with different options.

This means that while VCFcache guarantees functional equivalence (same variants receive same annotations), exact MD5 checksum identity between cached and uncached outputs cannot be guaranteed when using affected VEP versions.

**Impact**: This does not affect the correctness or usability of VCFcache - annotations are still accurate and consistent for each variant. It only affects byte-for-byte reproducibility testing via MD5 checksums.

**Reference**: This is a known VEP issue tracked at [ensembl-vep#1959](https://github.com/Ensembl/ensembl-vep/issues/1959).

**Recommendation**: For validation purposes, compare annotation content semantically (e.g., verify that the same variants have the same CSQ tags) rather than relying solely on MD5 checksums when using affected VEP versions.

---

## 12) CLI reference (all commands & flags)

This section is a compact reference for the CLI. It complements the walkthrough above.

**Note**: Commands with `--debug` flag now display detailed timing information when operations complete.

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
- `-o/--output`: output BCF file (use `-`/`stdout` to stream).
- `--stats-dir`: optional directory for logs/workflow/auxiliary. If provided, stats are written under `<stats_dir>/<input_basename>_vcstats`; if omitted, they default to `<cwd>/<input_basename>_vcstats`.
- `--no-stats`: disable stats/logs output (also disables `vcfcache compare`).
- `--md5-all`: compute full MD5 of all variants (no header) and store in stats (slow; may differ between runs).
- `-y/--yaml`: params YAML for runtime (copied to workflow `params.snapshot.yaml`).
- `-t/--threads`: threads for `bcftools` (default: 1).
- `-f/--force`: overwrite outputs.
- `--uncached`: run annotation without using the cache (debug/benchmark).
- `--parquet`: also convert output to Parquet for downstream querying.
- `--requirements`: print cache requirements and the stored annotation command, then exit.
- `--list`: list annotation caches under `-a` (when `-a` points to a directory of caches).
- `--debug`: keep intermediate files and use sandbox for alias resolution/downloads.
- `-v/--verbose`: increase logging (`-v` INFO, `-vv` DEBUG).

### `vcfcache list`

Discover public items on Zenodo or list local items.

**Required positional argument**:
- `{blueprints,caches}`: choose which type to list.

**Options**:
- `--genome GENOME`: filter public Zenodo results by genome keyword (e.g., `GRCh37`, `GRCh38`).
- `--source SOURCE`: filter public Zenodo results by source keyword (e.g., `gnomad`).
- `--local [PATH]`: list local items instead of querying Zenodo. Optionally specify a custom path (default: `VCFCACHE_DIR` or `~/.cache/vcfcache`).
- `--inspect PATH`: inspect a local cache/blueprint directory and print required `${params.*}` keys.
- `--debug`: query Zenodo sandbox instead of production.

**Examples**:
```bash
vcfcache list caches                    # List public caches from Zenodo
vcfcache list blueprints --local        # List local blueprints
vcfcache list caches --local /data      # List caches at custom path
vcfcache list --inspect /path/to/cache  # Inspect specific cache
```

**Note**: The old convenience aliases `vcfcache caches` and `vcfcache blueprints` have been removed. Use `vcfcache list caches` and `vcfcache list blueprints` instead.

### `vcfcache push`

Upload a local blueprint or cache directory to Zenodo (great for sharing caches).

- `--cache-dir`: path to a vcfcache base directory (blueprint upload) **or** a specific cache directory under `<base>/cache/<cache_name>` (cache+blueprint upload).
- `--dest`: upload destination (currently `zenodo`).
- `--test`: upload to Zenodo sandbox (uses `ZENODO_SANDBOX_TOKEN`).
- `--metadata`: YAML/JSON file with Zenodo metadata (title/description/creators/keywords, etc.). Required for real uploads.
- `--publish`: publish the record immediately after upload.
- `--yes`: skip confirmation prompt.

**What gets uploaded**
- **Blueprint upload** (base dir): tarball contains the full `blueprint/` plus an empty `cache/` directory.
- **Cache upload** (`<base>/cache/<name>`): tarball contains `blueprint/` and only the selected cache under `cache/<name>`.

**Zenodo metadata tips**
- Include keywords: `vcfcache`, `blueprint` or `cache`, and the alias (`bp-...` / `cache-...`) so `vcfcache list` can find it.
- The metadata templates in `resources/zenodo/` already include these keywords.

Environment variables:
- `ZENODO_TOKEN`: production token.
- `ZENODO_SANDBOX_TOKEN`: sandbox token.
- `ZENODO_SANDBOX=1`: affects some operations; for VCFcache CLI, `--debug` and `--test` are the main toggles.

### `vcfcache demo`

Convenience command for smoke testing and benchmarking.

**Two modes**:

1. **Smoke test mode**: Runs comprehensive test of all 4 commands using bundled demo data
   ```bash
   vcfcache demo --smoke-test
   ```

2. **Benchmark mode**: Compares cached vs uncached annotation performance
   ```bash
   vcfcache demo -a <cache> --vcf <file> -y <params> [--output <file>] [--stats-dir <dir>]
   ```

**Options**:
- `--smoke-test`: run smoke test mode.
- `-a/--annotation_db CACHE`: benchmark mode - path to annotation cache (requires `--vcf` and `-y`).
- `--vcf VCF`: benchmark mode - path to VCF/BCF file to annotate (requires `-a` and `-y`).
- `-y/--yaml PARAMS`: params YAML file (required for benchmark mode).
- `--output FILE`: benchmark mode - output BCF path (default: temporary file in /tmp).
- `--stats-dir DIR`: benchmark mode - stats directory (optional). If omitted, annotate defaults to `<cwd>/<input_basename>_vcstats`.
- `--debug`: keep intermediate files for inspection; also enables detailed timing output.
- `-q/--quiet`: suppress detailed output (show only essential information).

**Compressed help**: Running `vcfcache demo` without any options shows a brief usage guide. Use `vcfcache demo --help` for full help.
