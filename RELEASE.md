# Releasing `vcfcache`

This document describes the steps to create a new `vcfcache` release and publish
it to PyPI. The flow is:

1. Bump version
2. Update changelog & docs
3. Build & upload to **TestPyPI**
4. Smoke-test from TestPyPI
5. Upload to **PyPI**
6. Tag and create a GitHub release

> Notes:
> - The project is built with `hatchling` via `pyproject.toml`.
> - Distributions are built with `python -m build`.
> - Uploads are done with `twine`.
> - Python version and metadata live in `pyproject.toml`. :contentReference[oaicite:0]{index=0}

---

## 0. Prerequisites

- Python >= 3.11 installed
- `uv` or `pip` to manage the local environment
- Accounts on:
  - https://pypi.org
  - https://test.pypi.org
- API tokens created on both:
  - Go to “Account settings → API tokens”
  - Scope: whole account is fine for this project
- Locally:
  ```bash
  uv venv .venv
  source .venv/bin/activate
  uv pip install build twine
  # or: pip install build twine
  ```

---

## 1. Bump Version

Update the version in `pyproject.toml`:

```toml
[project]
version = "0.3.2"  # or whatever the new version is
```

Update `vcfcache/__init__.py`:

```python
__version__ = "0.3.2"
```

---

## 2. Update Changelog & Docs

Update `CHANGELOG.md` with the new version and release notes.

---

## 3. Build Package

Clean old builds and create new distributions:

```bash
rm -rf dist/ build/ *.egg-info
python -m build
```

This creates:
- `dist/vcfcache-X.Y.Z.tar.gz` (source distribution)
- `dist/vcfcache-X.Y.Z-py3-none-any.whl` (wheel)

Verify the contents:

```bash
# Check what's in the wheel
unzip -l dist/vcfcache-*.whl

# Check what's in the sdist
tar -tzf dist/vcfcache-*.tar.gz
```

---

## 4. Upload to TestPyPI

Test the release process on TestPyPI first:

```bash
python -m twine upload --repository testpypi dist/*
```

You'll be prompted for your TestPyPI API token.

---

## 5. Smoke Test from TestPyPI

Create a fresh environment and install from TestPyPI:

```bash
# Create test environment
uv venv /tmp/vcfcache-test-env
source /tmp/vcfcache-test-env/bin/activate

# Install from TestPyPI (dependencies from PyPI)
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple/ \
            vcfcache

# Run comprehensive demo
vcfcache demo

# Test basic functionality
vcfcache --help
vcfcache --version

# Clean up
deactivate
rm -rf /tmp/vcfcache-test-env
```

**Expected demo results:**
- ✓ Step 1: blueprint-init - Creates initial cache
- ✓ Step 2: blueprint-extend - Adds more variants
- ✓ Step 3: cache-build - Annotates the blueprint
- ✓ Step 4: annotate - Uses cache to annotate sample
- Shows variant counts and validates annotation tags

**If bcftools is not installed**, the demo will fail with instructions:
```
Install bcftools >= 1.20:
  Ubuntu/Debian: sudo apt-get install bcftools
  macOS: brew install bcftools
  Conda: conda install -c bioconda bcftools
```

---

## 6. Upload to PyPI

If TestPyPI testing passes, upload to the real PyPI:

```bash
python -m twine upload dist/*
```

You'll be prompted for your PyPI API token.

---

## 7. Verify PyPI Installation

Test installation from the real PyPI:

```bash
# Fresh environment
uv venv /tmp/vcfcache-pypi-test
source /tmp/vcfcache-pypi-test/bin/activate

# Install from PyPI
pip install vcfcache

# Run comprehensive demo
vcfcache demo

# Clean up
deactivate
rm -rf /tmp/vcfcache-pypi-test
```

---

## 8. Tag and Create GitHub Release

After successful PyPI upload:

```bash
# Create and push git tag
git tag v0.3.2  # match the version
git push origin v0.3.2

# Create GitHub release
gh release create v0.3.2 \
  --title "v0.3.2" \
  --notes "See CHANGELOG.md for details" \
  dist/*
```

Or create the release manually on GitHub:
1. Go to https://github.com/julius-muller/vcfcache/releases/new
2. Choose the tag `v0.3.2`
3. Add release notes from CHANGELOG.md
4. Attach `dist/*.whl` and `dist/*.tar.gz` files
5. Publish release

---

## Quick Reference

**Full release workflow:**

```bash
# 1. Update version in pyproject.toml and vcfcache/__init__.py
# 2. Update CHANGELOG.md

# 3. Build
rm -rf dist/ build/ *.egg-info
python -m build

# 4. Test on TestPyPI
python -m twine upload --repository testpypi dist/*

# 5. Demo from TestPyPI
uv venv /tmp/test-env && source /tmp/test-env/bin/activate
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple/ vcfcache
vcfcache demo
deactivate && rm -rf /tmp/test-env

# 6. Upload to PyPI
python -m twine upload dist/*

# 7. Create git tag and GitHub release
git tag v0.3.2 && git push origin v0.3.2
gh release create v0.3.2 --title "v0.3.2" --notes-file CHANGELOG.md dist/*
```

---

## Troubleshooting

### Demo fails: "bcftools not found"

The demo requires bcftools >= 1.20 to be installed. This is a runtime dependency of vcfcache.

Users should install bcftools before using vcfcache:
- Ubuntu/Debian: `sudo apt-get install bcftools`
- macOS: `brew install bcftools`
- Conda: `conda install -c bioconda bcftools`

### Upload fails: "File already exists"

You cannot overwrite a version on PyPI once uploaded. You need to:
1. Bump the version number
2. Rebuild: `python -m build`
3. Upload the new version

### Wheel doesn't include expected files

Check `pyproject.toml` build configuration:
```toml
[tool.hatch.build.targets.wheel]
packages = ["vcfcache"]
include = ["vcfcache/recipes/**"]
```

Verify with: `unzip -l dist/vcfcache-*.whl`
