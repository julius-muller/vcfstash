# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added
- Enforced genome build metadata in configuration:
  - `genome_build` is now required in both `params.yaml` and `annotation.yaml`
  - Workflow startup logs the genome build from both configs
  - Mismatched genome builds now fail fast

### Changed
- **BREAKING:** Removed all contig renaming/variant-generation support.
  - The `--generate-contig-variants` option is removed.
  - Annotation now logs contig overlap at startup and fails if there is no overlap.
- **BREAKING:** `vcfcache annotate` now writes directly to an output BCF file instead of creating an output directory.
  - Use `--stats-dir` to store workflow logs/auxiliary files under `<output_file>_vcstats`.
  - Use `-` or `stdout` to stream output to stdout.
  - Any tooling or scripts that expected an output directory now need to pass a file path (e.g., `sample_vc.bcf`).
- CI workflow now runs on both tag pushes and main branch pushes (for codecov badge updates)
- Dockerfile optimizations: moved tests to test stage only, added pip cache mounts for faster builds

### Fixed
- Contig overlap is now reported on every cached annotation run and no-overlap failures are surfaced early.
- Fixed `AttributeError: 'Namespace' object has no attribute 'verbose'` when running `vcfcache compare` command
  - Added missing parent parser to compare command to inherit `--verbose`, `--quiet`, and `--debug` flags
- Enhanced `vcfcache compare` to show detailed timing breakdown similar to old benchmark mode
  - Now displays detailed step-by-step timings for each bcftools command in both cached and uncached runs
  - Shows comprehensive summary with total time, speedup factor, time saved (in seconds, minutes, and hours)
  - Improved workflow.log parsing to correctly extract "Workflow completed successfully in X.Xs" timing
  - Human-readable time formatting (e.g., "1h 23m 45.6s" or "12m 34.5s")
  - **NEW:** Displays genome_build and threads from params.snapshot.yaml for both runs
  - **NEW:** Counts and compares total variants in output files using `bcftools index -n`
  - **NEW:** Shows both variant count comparison and MD5 checksum verification
  - Removed warning about comparing runs in same mode (useful for comparing different caches)
  - Added comprehensive test suite for compare command (`tests/test_compare.py`) with 16 tests

## 0.4.1 (2026-01-06)

### Added
- Added formal YAML schema validation system (`vcfcache/utils/schemas.py`)
  - `ParamsYAMLSchema` class validates params.yaml structure with required/optional/forbidden fields
  - `AnnotationYAMLSchema` class validates annotation.yaml structure with required/optional/forbidden fields
  - Schema validation integrated into workflow loading with clear error messages
  - Detects when wrong YAML file type is provided (e.g., annotation.yaml passed to `-y` flag)
- Added comprehensive schema validation test suite (`tests/test_schemas.py`) with 23 tests covering all validation scenarios
- Added informative message when `vcfcache annotate` uses default params.snapshot.yaml from cache (when `-y` not provided)
- Added `threads` as required field in params.yaml (must be integer >= 1)
- Required fields in params.yaml: `annotation_tool_cmd`, `bcftools_cmd`, `temp_dir`, `threads`, `optional_checks`
- Required fields in annotation.yaml: `annotation_cmd`, `must_contain_info_tag`, `required_tool_version`, `optional_checks`

### Changed
- Moved `threads` configuration from CLI option (`-t/--threads`) to required field in params.yaml
  - Removed `-t/--threads` from `blueprint-init`, `blueprint-extend`, and `cache-build` commands
  - Auto-generated params files now include `threads: 1` as default
  - Ensures consistent configuration through params.yaml with schema validation
- Updated `-a` help text in `vcfcache annotate` to clarify it requires specific cache directory, not cache root
- Updated `-y` help text to explicitly mention default behavior (uses cache's params.snapshot.yaml if not provided)

### Fixed
- Fixed bug where accidentally providing annotation.yaml via `-y` flag (instead of params.yaml) would cause uncached annotation to run without proper validation
- Schema validation now provides helpful error messages when wrong YAML type is detected (e.g., "Did you accidentally provide annotation.yaml instead of params.yaml?")
- Fixed test fixtures in `test_contig_mismatches.py` to include all required schema fields for params.yaml and annotation.yaml

## 0.4.0 (2026-01-05)

### Added
- Added detailed timing breakdown to `blueprint-init`, `blueprint-extend`, `cache-build`, and `annotate` commands when `--debug` flag is provided (similar to `demo --smoke-test` output)
- Added comprehensive performance model documentation to WIKI.md explaining runtime efficiency, cache lookup mechanics, and expected speed-ups (2-10× for typical samples with 60-90% cache hits)
- Added brief performance summary to README.md with link to detailed model
- Added caveat section explaining when caching provides less benefit (small datasets, fast pipelines)
- Added troubleshooting note about VEP non-deterministic output in recent versions (ensembl-vep#1959)

### Changed
- Simplified `list` command interface:
  - Removed `vcfcache caches` and `vcfcache blueprints` convenience aliases (use `vcfcache list caches` and `vcfcache list blueprints` instead)
  - Made selector (`blueprints`/`caches`) a required positional argument
  - Simplified options: `--local [PATH]` now optionally takes a path directly
- Improved `demo` command help text to clearly distinguish smoke test mode vs benchmark mode
- `vcfcache demo` and `vcfcache list` now show compressed help when run without required arguments (full help available with `--help`)
- Updated all documentation to reflect new CLI interface and performance model

## 0.4.0b2 (2025-12-16)

### Added
- Added `vcfcache caches` / `vcfcache blueprints` shortcuts and `--local` listing (with optional `--path`) to discover locally available items.
- Added `vcfcache list --inspect <path-or-alias>` to show cache details and the exact `${params.*}` keys required by `annotation.yaml`.
- Improved Zenodo discovery output: formatted cache titles from alias keywords and ignored placeholder records smaller than 1MB.
- Added `-o/--output` to `vcfcache cache-build --doi ...` to control the download base directory per command.
- Added helper script to extract release notes from `CHANGELOG.md` (`scripts/release_notes_from_changelog.py`) and use it for GitHub Releases.

### Changed
- Overhauled user documentation (`WIKI.md`) into a walkthrough + complete CLI reference and clarified genome-agnostic support.
- Simplified GitHub Actions to a single tag-triggered CI workflow (no automated release/publish workflows).
- Updated Docker documentation and local build docs to reflect the current single runtime image (`docker/Dockerfile.vcfcache`).

### Fixed
- Made MD5 computation pure Python (no external `md5sum` dependency) for portability.
- Improved Zenodo sandbox integration test isolation (avoid interference from global `VCFCACHE_DIR`).
- Hardened CI builds: retry/fallback downloads for bcftools/htslib and cached the compiled toolchain to avoid rebuilding on every run.

## 0.4.0b2 (2025-12-16)

### Added
- Mainly adaptations in documentation and CI workflows to prepare for release.

## 0.4.0b1 (2025-12-15)

### Added
- Added tag-based GitHub Actions release workflow (`.github/workflows/release.yml`) to build/test, publish to (Test)PyPI, create GitHub Releases, and build/push Docker images.
- Updated `scripts/release.sh` to support PEP 440 pre-releases (e.g. `0.4.0b0`, `0.4.0rc1`) and mark GitHub releases as pre-release accordingly.
- Documented pre-release (beta/rc) guidance in `RELEASE.md`.

### Fixed
- Removed dependency on external `md5sum` binary by switching MD5 computation to pure Python (`hashlib`).
- Made Zenodo sandbox integration test robust to a globally set `VCFCACHE_DIR` by forcing an isolated cache dir under the test HOME.
- Fixed CI `PATH` handling so system tools remain available while still preferring the locally built `bcftools` binaries.
- Added retries/fallback URLs for downloading `bcftools`/`htslib` sources in GitHub Actions.
- Cached built `bcftools`/`htslib` in GitHub Actions to avoid recompiling on every run.
- Restricted CI to run only on `v*` tags (to avoid expensive builds on frequent pushes).

## 0.3.6 (2025-12-11)

### Fixed
- Fixed cache path resolution when using `vcfcache annotate -a <alias>` with downloaded caches
- Correctly locate `annotation.yaml` in extracted cache directory structure

## 0.3.5 (2025-12-11)

### Fixed
- Fixed packaging issue where `public_caches.yaml` was not included in PyPI wheel distribution
- Moved `public_caches.yaml` from project root to `vcfcache/` package directory to ensure proper installation
- `vcfcache list --public-caches` now works correctly when installed from PyPI

## 0.3.4 (2025-12-11)

### Added
- Added `--quiet/-q` flag to `vcfcache demo` command for minimal output
- Added automated release script (`scripts/release.sh`) with interactive prompts and smart version checking
- Added environment variable `VCFCACHE_BCFTOOLS` to override system bcftools location

### Changed
- Updated minimum bcftools requirement from 1.16 to 1.20 (for `--write-index` flag support)
- Improved bcftools error messages to mention override option

### Documentation
- Added bcftools configuration section to README.md, WIKI.md, and CLAUDE.md
- Updated RELEASE.md with automated release workflow

## 0.3.3 (2025-12-09)

### Added
- Extended `vcfcache demo` command with multiple modes:
  - `--smoke-test` flag for comprehensive testing of all 4 commands (blueprint-init, blueprint-extend, cache-build, annotate)
  - Benchmark mode (`-a <cache> --vcf <vcf> -y <params>`) for comparing cached vs uncached annotation performance
  - Benchmark mode shows detailed timing comparison, speedup calculation, and MD5 checksums of outputs
  - Demo command now shows help when invoked without arguments

### Changed
- Renamed internal `run_demo()` function to `run_smoke_test()` for clarity

### Fixed
- Fixed version lookup to handle uninstalled package (fallback to `__version__` from `__init__.py`)
- Fixed argparse conflict with parent parser arguments

## 0.3.2 (2025-12-09)

### Added
- Added demo 

## 0.3.1 (2025-12-09)

### Added
- Added zenodo as cache store
- Removed fat batteries included docker images in favor of relying on external annotation pipelines
- Removed nextflow

## 0.3.0 (2025-12-03)

### Added
- Comprehensive testing infrastructure with scenario-aware tests (vanilla, blueprint, annotated)
- Detailed testing documentation in `tests/README.md`
- Annotation workflow integration tests for all three scenarios
- Benchmark scripts for performance testing on BCF subsets (`tests/run_benchmarks_gnomad.sh`)
- Support for temporary output directories in benchmark scripts
- Dev dependencies (`pytest`, testing tools) to blueprint Docker images via `Dockerfile.cache`

### Changed
- None

### Deprecated
- None

### Removed
- None

### Fixed
- AF (allele frequency) threshold calculation in Docker build scripts - fixed octal interpretation bug where `af010` was incorrectly parsed as octal (8) instead of decimal (10)
- AF threshold display now correctly shows 0.01 for `af001` (1%) and 0.10 for `af010` (10%)
- Blueprint test path expectations - tests now correctly expect `cache_dir/cache/annotation_name` structure instead of incorrectly assuming `cache_dir/db/cache/annotation_name` (the `/db` layer only exists in Docker images, not in directly created caches)

### Security
- None

## 0.2.0 (2025-05-09)

### Added
- Complete rewrite including documentation and examples.
 
### Changed
- None

### Deprecated
- None

### Removed
- None

### Fixed
- None

### Security
- None

## 0.1.0 (2024-04-09)

### Added
- First release
 
### Changed
- None

### Deprecated
- None

### Removed
- None

### Fixed
- None

### Security
- None
