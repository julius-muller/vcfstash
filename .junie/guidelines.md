# VCFstash Development Guidelines

## Project Overview
VCFstash is a tool to accelerate VCF annotations of large VCF files by maintaining a cache of frequently shared variants across human WGS samples. It manages a variant cache database and runs VCF annotations only on novel variants not present in the cache, significantly reducing annotation time.

## Build/Configuration Instructions

### Requirements
- Python 3.11+
- bcftools (version 1.20 shipped in directory `tools/bcftools`)
- Nextflow (version 24+, shipped in directory `workflow/.nextflow/framework/24.10.5/nextflow-24.10.5-one.jar`)

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/julius-muller/vcfstash.git
   cd vcfstash
   ```

2. Create and activate a virtual environment:
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

4. For Parquet support (optional):
   ```bash
   pip install -r requirements_parquet.txt
   ```

5. Install the package in development mode:
   ```bash
   pip install -e .
   ```

### Configuration
- Users typically have their annotation workflow (e.g. vep, annovar ...) setup as bash executable commands 
- The user copies the annotation workflow as bash script in variable annotation_cmd within the annotation.config file (test example in `tests/config/annotation.config`)
- After initialization of the database with a suitable vcf containing highly abundant variants, the normalized database containing only variants is stored as blueprint and cann be extended using `stash-add` command
- The user can run the annotation workflow on the blueprint using the `stash-annotate` command which instantiates a cache for vcf annotation.
- All parameters (paths, resources) in users annotation workflow that need to be configurable when running the workflow on vcf files (`annotate`), should be replaced by environment variables in annotation_cmd and variable definitions should be set in the YAML file
- The user can run the annotation workflow on a vcf file using the `annotate` command using the non-configurable annotation.config file from the `stash-annotate` step
- The application requires the VCFSTASH_ROOT environment which is set automatically in vcfstash.py

## Testing Information

### Test Structure
- Tests are located in the `tests/` directory
- Test files in directory `tests/data/expected_output` are created using `tests/update_reference.py` and should be used as gold standard for new tests
- The project uses pytest as the testing framework
- Test data is stored in `tests/data/nodata/` and `tests/data/expected_output/`
- Test configurations are in `tests/config/`

### Running Tests
1. Ensure you're in the project root directory with the virtual environment activated
2. Run all tests:
   ```bash
   python -m pytest
   ```
3. Run specific tests:
   ```bash
   python -m pytest tests/test_validation.py
   ```
4. Run tests with increased verbosity:
   ```bash
   python -m pytest -xvs tests/test_validation.py
   ```

### Adding New Tests
1. Create a new test file in the `tests/` directory with a name starting with `test_`
2. Import the necessary modules and functions
3. Write test functions with names starting with `test_`
4. Use pytest assertions to verify expected behavior

### Example Test
Here's a simple test for the `compute_md5` function:

```python
import os
import tempfile
import pytest
from pathlib import Path
from src.utils.validation import compute_md5

def test_compute_md5():
    """Test that compute_md5 correctly calculates MD5 hash of a file."""
    # Create a temporary file with known content
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file.write(b"test content for MD5 calculation")
        temp_file_path = temp_file.name

    try:
        # Calculate MD5 hash
        calculated_md5 = compute_md5(Path(temp_file_path))

        # Expected MD5 hash for "test content for MD5 calculation"
        expected_md5 = "11c1ea4414f9b160b0b9f98a3e53f3a2"

        # Assert that the calculated hash matches the expected hash
        assert calculated_md5 == expected_md5, f"Expected {expected_md5}, got {calculated_md5}"
    finally:
        # Clean up the temporary file
        os.unlink(temp_file_path)
```

## Additional Development Information

### Project Structure
- `src/`: Source code
  - `src/database/`: Database operations (initializer, updater, annotator)
  - `src/utils/`: Utility functions (validation, logging)
- `tests/`: Test files and data
- `workflow/`: Nextflow workflow files
- `resources/`: Resource files
- `tools/`: External tools
- `scripts/`: Utility scripts
- `vcfstash.py`: Main entry point

### Command-Line Interface
The application provides four main commands:
1. `stash-init`: Initialize a new database
2. `stash-add`: Add new VCF data to an existing database
3. `stash-annotate`: Run annotation workflow on the database
4. `annotate`: Annotate a VCF file using the cached database

### Dependencies
- Core: pysam, PyYAML
- Testing: pytest
- Optional: pandas, pyarrow (for Parquet support)

### Build System
- The project uses hatchling as the build backend
- The package is configured in pyproject.toml
- Entry points are defined in setup.py

### Environment Variables
- `VCFSTASH_ROOT`: Should be set to the project root directory (automatically set during tests)
