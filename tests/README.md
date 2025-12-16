# VCFcache Test Suite

This file is the canonical test documentation for VCFcache.

## Quick commands

Host (development):
```bash
python -m pytest tests/ -vv --color=yes --disable-warnings -rA
```

Docker (test stage from `docker/Dockerfile.vcfcache`):
```bash
docker build --target test -f docker/Dockerfile.vcfcache -t vcfcache:test .
docker run --rm vcfcache:test
```

---

## Philosophy: One suite, scenario-aware

## The Three Scenarios

### 1) Vanilla (dev environment)
**What it is:** Plain Python environment with no pre-built cache

**What's available:**
- ✅ Python package with vcfcache installed
- ✅ `bcftools` available on `PATH`
- ✅ Test data files (`tests/data/`)
- ❌ NO pre-built cache at `/cache`
- ❌ NO VEP

**How tests adapt:**
- Create cache from scratch using test data
- Use mock annotation (bcftools annotate)
- Full workflow testing from blueprint-init through annotate
- Cache validation tests are skipped (no `/cache`)

**Run (host):**
```bash
python -m pytest tests/ -vv --color=yes --disable-warnings -rA
```

---

### 2) Blueprint (data-only image)
**What it is:** Environment with a pre-built cache at `/cache`, but no VEP

**What's available:**
- ✅ Pre-built cache at `/cache` (gnomAD variants, AF-filtered)
- ✅ Test data files (copied during build)
- ❌ NO VEP

**How tests adapt:**
- Validate `/cache` structure and integrity
- Create mini test cache from test data for workflow testing
- Use mock annotation (bcftools annotate)
- Full workflow testing with mock annotation
- Tests verify cache is production-ready

**Run (example):**
Build a custom container/environment that has:
- VCFcache installed
- test data available
- a prebuilt cache mounted at `/cache`

Then run:
```bash
PYTEST_ADDOPTS="--color=yes --disable-warnings -rA --durations=10 -vv --maxfail=1" \
python -m pytest tests
```

---

### 3) Annotated (full-stack image)
**What it is:** Environment with VEP available and a pre-built cache at `/cache`

**What's available:**
- ✅ Pre-built cache at `/cache` (gnomAD + VEP annotations)
- ✅ VEP (Ensembl Variant Effect Predictor)
- ✅ Test data files (copied during build)

**How tests adapt:**
- Validate `/cache` structure and integrity
- Can use REAL VEP annotation (currently uses mock for testing)
- Full workflow testing with annotation
- Complete integration testing

**Run (example):**
In an environment where:
- `/cache` exists (prebuilt cache)
- `vep` is available on `PATH` (optional for some tests)

Run:
```bash
PYTEST_ADDOPTS="--color=yes --disable-warnings -rA --durations=10 -vv --maxfail=1" \
python -m pytest tests
```
Tip: add `-s` inside `PYTEST_ADDOPTS` if you want live stdout.

---

## How Scenario Detection Works (tl;dr)

Tests automatically detect which scenario they're running in:

```python
# In conftest.py
@pytest.fixture(scope="session")
def test_scenario():
    """Detect scenario based on environment."""
    cache_exists = Path("/cache").exists()
    vep_available = shutil.which("vep") is not None

    if cache_exists and vep_available:
        return "annotated"
    elif cache_exists:
        return "blueprint"
    else:
        return "vanilla"
```

Tests receive this as a fixture and adapt:

```python
def test_full_annotation_workflow(test_scenario, prebuilt_cache):
    """Adapts based on scenario."""
    if test_scenario == "vanilla":
        # Create cache from scratch
        run_blueprint_init(TEST_VCF, output_dir)
    else:
        # Use prebuilt cache
        print(f"Using cache at {prebuilt_cache}")

    # Rest of test proceeds...
```

---

## Test Categories (where to look)

### Pure Python Tests
Tests that work in all scenarios without adaptation:
- Module imports
- CLI help and version
- Error handling
- Utility functions (MD5, path resolution)
- Configuration validation

**Files:** `test_vanilla.py`

### Cache Validation Tests
Tests that validate cache structure (skip in vanilla scenario):
- Cache directory structure
- Blueprint BCF file integrity
- Cache metadata
- Workflow directory
- Cache query performance
- Chromosome naming conventions

**Files:** `test_blueprint.py`

**Behavior:**
- Vanilla: Skipped (no cache)
- Blueprint/Annotated: Validate `/cache`

### Workflow Integration Tests

### Performance / Cache Stats
- In annotated scenario only; long‑running (~30–55s each on Docker)
- `test_annotation.py::test_cache_hit_statistics`
- `test_annotation.py::test_annotation_performance_with_cache`
Tests that run full annotation workflows:
- blueprint-init, blueprint-extend, cache-build, annotate
- Cached vs uncached annotation comparison
- Input file preservation
- Normalization flag
- File validation

**Files:** `test_annotate.py`, `test_core.py`

**Behavior:**
- Vanilla: Create cache from test data, use mock annotation
- Blueprint: Create test cache from test data, use mock annotation, validate prebuilt cache
- Annotated: Create test cache from test data, can use VEP or mock annotation

---

## Running Tests

### Run all tests (auto-detects scenario)
```bash
# In development (vanilla)
pytest tests/ -v

# In Docker test stage
docker build --target test -f docker/Dockerfile.vcfcache -t vcfcache:test .
docker run --rm vcfcache:test
```

### Run with coverage
```bash
pytest tests/ --cov=vcfcache --cov-report=html -v
```

### Run specific test file
```bash
pytest tests/test_annotate.py -v
pytest tests/test_blueprint.py -v
pytest tests/test_vanilla.py -v
```

---

## Expected Test Results

| Scenario | Total Tests | Passed | Skipped | Failed |
|----------|------------|--------|---------|--------|
| **Vanilla** | 29 | ~18 | ~11 (cache tests) | 0 |
| **Blueprint** | 29 | ~29 | 0 | 0 |
| **Annotated** | 29 | ~29 | 0 | 0 |

---

## Test File Structure

```
tests/
├── README.md                    # This file
├── conftest.py                  # Scenario detection & fixtures
├── test_vanilla.py              # Pure Python tests (10 tests)
├── test_blueprint.py            # Cache validation tests (11 tests)
├── test_annotate.py             # Workflow tests (5 tests)
├── test_core.py                 # Core functionality tests (3 tests)
├── config/                      # Test configurations
│   ├── test_params.yaml         # Mock annotation params
│   └── test_annotation.config   # Mock annotation command
└── data/                        # Test data (small BCF files)
    ├── nodata/
    │   ├── crayz_db.bcf         # Test VCF 1
    │   ├── crayz_db2.bcf        # Test VCF 2
    │   └── sample4.bcf          # Test sample
    └── references/
        └── reference.fasta      # Mini reference genome
```

---

## Key Fixtures

### Scenario Detection
- `test_scenario`: Returns "vanilla", "blueprint", or "annotated"
- `prebuilt_cache`: Path to `/cache` if available, None otherwise
- `use_prebuilt_cache`: Boolean flag

### Test Data
- `params_file`: Temporary params file with paths resolved
- `test_output_dir`: Clean temporary directory for test outputs
- `annotation_config`: Path to mock annotation config

### Blueprint Helpers
- `mini_cache_dir`: Creates mini test cache from top 10 variants of blueprint (blueprint scenario only)

---

## Adding New Tests

All new tests should be scenario-aware:

```python
def test_my_new_feature(test_scenario, prebuilt_cache):
    """My new test that adapts to scenario."""
    print(f"Testing in {test_scenario} scenario")

    if test_scenario == "vanilla":
        # Vanilla-specific logic
        cache_dir = create_cache_from_scratch()
    else:
        # Blueprint/Annotated logic
        cache_dir = use_existing_cache(prebuilt_cache)

    # Common test logic
    result = run_test(cache_dir)
    assert result.success
```

### For cache validation tests:
```python
def test_cache_feature(test_scenario):
    """Test cache feature (skips in vanilla)."""
    if test_scenario == "vanilla":
        pytest.skip("No cache in vanilla scenario")

    # Cache validation logic
    assert Path("/cache").exists()
```

---

## Troubleshooting

### Tests fail with "Vanilla scenario has no pre-built cache"
**Expected behavior.** Cache validation tests skip in vanilla scenario since there's no `/cache` directory.

### All tests fail immediately
- Check Python environment: `python3 -m pytest --version`
- Verify installation: `pip show vcfcache`
- Check PYTHONPATH in Docker: `export PYTHONPATH=/app/venv/lib/python3.13/site-packages`

### "bcftools not found" errors
- Install `bcftools >= 1.20`, or set `VCFCACHE_BCFTOOLS=/path/to/bcftools`.
- In Docker, the runtime image installs `bcftools`; in the test stage image, tests run via `docker run --rm vcfcache:test`.

### Tests hang or timeout
- Increase timeout in pytest: `pytest --timeout=300`

---

## CI/CD Integration

### CI note

The repository CI workflow lives in `.github/workflows/ci.yml` and runs tests on version tags.
If you build your own “blueprint/annotated” images with `/cache` and/or VEP available, the same test suite will adapt automatically.

---

## Development Workflow

1. **Write tests:** Add scenario-aware tests to appropriate file
2. **Test locally (vanilla):** `pytest tests/ -v`
3. **Test in blueprint:** Build blueprint image, run tests
4. **Test in annotated:** Build annotated image, run tests
5. **Verify all scenarios pass:** All 29 tests should pass (with some skipped in vanilla)

---

## Philosophy Summary

**One codebase, three scenarios, 100% coverage everywhere.**

- Tests adapt to available resources
- Each scenario demonstrates full functionality
- Users can validate any scenario works end-to-end
- No duplication, maximum confidence
