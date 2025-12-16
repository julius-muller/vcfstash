# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

## 0.4.0 (2025-12-16)

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
