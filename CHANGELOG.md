# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
- Detailed testing documentation in `tests/README.md` and `tests/README_TESTING.md`
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