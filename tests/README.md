# VCFstash Test Suite

## Test Organization

Tests are organized using pytest markers:

- **`@pytest.mark.blueprint`**: Tests that should pass in the vcfstash-blueprint Docker image (lean, data-only)
- **`@pytest.mark.unit`**: Unit tests with no external dependencies
- **`@pytest.mark.integration`**: Integration tests that require Nextflow and annotation tools

## Running Tests

### Run all tests
```bash
pytest tests/ -v
```

### Run only blueprint-compatible tests (for Docker image validation)
```bash
pytest -m blueprint tests/ -v
```

Alternatively, exclude integration tests:
```bash
pytest -m "not integration" tests/ -v
```

### Run only integration tests (requires Nextflow)
```bash
pytest -m integration tests/ -v
```

### Run only unit tests
```bash
pytest -m unit tests/ -v
```

## Docker Image Testing

For the **vcfstash-blueprint** Docker image (lean, data-only), use:

```bash
docker run --rm --entrypoint /bin/sh \
  ghcr.io/julius-muller/vcfstash-blueprint:TAG \
  -c 'cd /app && export PYTHONPATH=/app/venv/lib/python3.13/site-packages:$PYTHONPATH && python3 -m pytest -m blueprint tests/ -v'
```

Expected result: **5 tests should pass** (all blueprint-marked tests)

## Test Breakdown

### Blueprint Tests (5 tests - should pass in Docker)
- `test_annotate.py::test_sample_file_validity` - Validates BCF file integrity
- `test_core.py::test_error_handling` - Tests error conditions
- `test_core.py::test_file_validation` - Tests file validation and MD5 checks
- `test_core.py::test_vcf_reference_validation` - Tests VCF reference validation

### Integration Tests (4 tests - require Nextflow)
- `test_annotate.py::test_full_annotation_workflow` - Full pipeline from stash-init to annotate
- `test_annotate.py::test_cached_vs_uncached_annotation` - Compares cached vs uncached results
- `test_annotate.py::test_input_not_modified_during_annotation` - Verifies input preservation
- `test_annotate.py::test_normalization_flag` - Tests normalization functionality
