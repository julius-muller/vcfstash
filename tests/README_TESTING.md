# VCFstash Testing Guide

This document describes the testing suite for VCFstash, including how tests adapt to different scenarios (vanilla, blueprint, annotated).

## Test Scenarios

VCFstash tests automatically detect and adapt to three scenarios:

### 1. **Vanilla** (No pre-built cache)
- Standard Python environment
- No `/cache` directory
- Tests cache creation from scratch
- Uses mock annotation (bcftools) for testing

### 2. **Blueprint** (Data-only Docker)
- Pre-built cache at `/cache`
- No VEP available
- Tests cache validation and usage
- Uses mock annotation (bcftools) for testing

### 3. **Annotated** (Full Docker with VEP)
- Pre-built cache at `/cache`
- VEP available
- Pre-annotated stash at `/cache/db/stash/vep_gnomad`
- Uses real VEP annotation for testing

## Test Files

### `test_blueprint.py`
Tests cache structure validation across scenarios:
- Cache directory structure
- Blueprint BCF file integrity
- Workflow files
- Chromosome naming conventions
- Python environment setup

**Applicable scenarios:** Blueprint, Annotated (skips Vanilla)

### `test_annotation.py`
Tests actual annotation workflows across all scenarios:

#### Annotated Scenario Tests:
- `test_annotated_stash_exists()` - Verifies pre-built VEP stash
- `test_annotate_with_cache()` - Tests annotation using pre-built cache
- `test_cache_hit_statistics()` - Verifies cache hit tracking
- `test_annotation_performance_with_cache()` - Performance measurement

#### Blueprint Scenario Tests:
- `test_blueprint_annotation_workflow()` - Full workflow: init → stash-annotate → annotate
- Uses mock annotation config

#### Vanilla Scenario Tests:
- `test_vanilla_annotation_workflow()` - Direct annotation without cache
- Simulates annotation behavior

#### Cross-Scenario Tests:
- `test_annotation_consistency_across_scenarios()` - Verifies MD5sum consistency
- Key test: annotations should be identical regardless of cache

**Applicable scenarios:** All (adapts behavior per scenario)

## Running Tests

### Run all tests (auto-detects scenario):
```bash
python -m pytest tests/
```

### Run specific test file:
```bash
python -m pytest tests/test_annotation.py -xvs
```

### Run specific test:
```bash
python -m pytest tests/test_annotation.py::test_annotate_with_cache -xvs
```

### Run with verbose output:
```bash
python -m pytest tests/ -vv
```

## Test Data

Test data is located in `tests/data/nodata/`:

- `gnomad_test.bcf` - Small gnomAD subset (cache hits in blueprint/annotated)
- `sample4.bcf` - Test sample variants (cache misses)
- `dbsnp_test.bcf` - Additional test data
- `crayz_db*.bcf` - Edge case test data

The fixture `test_sample_with_hits_and_misses()` in `test_annotation.py` combines these to create samples with both cache hits and misses for testing.

## Annotation Configs

Located in `tests/config/`:

### `test_annotation.config`
- Mock annotation using bcftools
- Adds `MOCK_ANNO` INFO tag
- Used for blueprint and vanilla scenarios
- Fast, no external dependencies

### `example_annotation.config`
- Real VEP annotation configuration
- Adds `CSQ` INFO tag with full annotations
- Used for annotated scenario only
- Requires VEP and cache directory

## Test Fixtures (conftest.py)

Key fixtures available to all tests:

- `test_scenario` - Auto-detects scenario ('vanilla', 'blueprint', 'annotated')
- `prebuilt_cache` - Path to `/cache` if available
- `annotation_config` - Appropriate config for scenario (VEP or mock)
- `test_output_dir` - Temporary directory for test outputs
- `test_sample_with_hits_and_misses` - Test VCF with cache hits and misses
- `annotation_stash_path` - Path to VEP stash in annotated scenario

## Key Testing Principles

### 1. Scenario Adaptation
Tests automatically skip when not applicable:
```python
if test_scenario != "annotated":
    pytest.skip("Only applicable to annotated scenario")
```

### 2. Cache Hits vs Misses
Tests verify both cached and non-cached variants are handled correctly:
- Cached variants (from gnomad_test.bcf) should use pre-computed annotations
- Non-cached variants (from sample4.bcf) should be annotated on-the-fly

### 3. Annotation Consistency
**Critical principle:** Annotation results must be MD5sum identical regardless of cache size or presence. The only difference should be annotation speed.

Tests verify:
- Same input → same output
- Cache only affects performance, not results
- No variants are missed or incorrectly annotated

### 4. Cleanup on Failure
Test output directories are preserved on failure for forensic analysis:
```
TEST FAILED! Temporary directory NOT removed for forensic analysis:
    /tmp/vcfstash_test_xxxxx
```

## Docker Testing

### Blueprint Docker:
```bash
docker run --rm ghcr.io/julius-muller/vcfstash-blueprint:gnomad-grch38-joint-af010 \
  python3 -m pytest tests/test_blueprint.py -xvs
```

### Annotated Docker:
```bash
docker run --rm --entrypoint='' ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af010-vep115 \
  bash -c 'cd /app && export PYTHONPATH=/app/venv/lib/python3.13/site-packages && python3 -m pytest tests/ -xvs'
```

Or use the default bash entrypoint:
```bash
docker run --rm ghcr.io/julius-muller/vcfstash-annotated:gnomad-v41-grch38-joint-af010-vep115 \
  -c 'cd /app && export PYTHONPATH=/app/venv/lib/python3.13/site-packages && python3 -m pytest tests/ -xvs'
```

## Expected Test Behavior by Scenario

| Test | Vanilla | Blueprint | Annotated |
|------|---------|-----------|-----------|
| `test_cache_directory_exists` | Skip | ✓ Pass | ✓ Pass |
| `test_annotated_stash_exists` | Skip | Skip | ✓ Pass |
| `test_annotate_with_cache` | Skip | Skip | ✓ Pass |
| `test_blueprint_annotation_workflow` | Skip | ✓ Pass | Skip |
| `test_vanilla_annotation_workflow` | ✓ Pass | Skip | Skip |
| `test_annotation_consistency` | ✓ Pass | ✓ Pass | ✓ Pass |

## Debugging Failed Tests

1. **Check test output directory** - preserved on failure
2. **Review stderr/stdout** - captured by pytest
3. **Verify scenario detection** - `print(test_scenario)` in test
4. **Check file paths** - use absolute paths for Docker
5. **Verify tools available** - bcftools, vep (if annotated)

## Adding New Tests

When adding tests:

1. Consider which scenarios they apply to
2. Use `pytest.skip()` for inapplicable scenarios
3. Add appropriate fixtures
4. Document expected behavior per scenario
5. Verify cleanup on both success and failure

Example:
```python
def test_new_feature(test_scenario, prebuilt_cache):
    """Test new feature."""
    if test_scenario == "vanilla":
        pytest.skip("Requires pre-built cache")

    # Test implementation
    assert cache_feature_works()
```

## Troubleshooting

### Tests fail to find bcftools
→ Check `setup_test_environment` fixture is running (autouse=True)

### Tests timeout in annotated scenario
→ VEP annotation is slow; increase timeout in subprocess.run()

### Cache hit/miss detection not working
→ Verify variants in test data overlap with cache variants

### MD5 sums don't match across scenarios
→ This is a **critical bug** - annotation should be deterministic
