# Releasing `vcfcache`

Internal release checklist for publishing to PyPI.

---

## Model (maintainers)

- Develop on `main` (keep it green).
- Cut releases from immutable git tags (`v<version>`).
- Keep PyPI + Docker + GitHub Release artifacts aligned to the tag commit.

### Versioning (beta/rc/final)

Use PEP 440 versions:
- Beta: `0.4.0b0`, `0.4.0b1`, ...
- Release candidate: `0.4.0rc1`, `0.4.0rc2`, ...
- Final: `0.4.0`

Tag naming:
- `v0.4.0b0` / `v0.4.0rc1` / `v0.4.0`

---

## Quick Reference

**For a new release:**

```bash
# 1. Update version: pyproject.toml, __init__.py, CHANGELOG.md

# 2. Build & test locally
rm -rf dist/ build/ *.egg-info && python -m build
uv venv /tmp/t && source /tmp/t/bin/activate
uv pip install dist/vcfcache-*.whl
vcfcache demo --smoke-test && python -m pytest tests -q
deactivate && rm -rf /tmp/t

# 3. Upload to TestPyPI & verify
python -m twine upload --repository testpypi dist/*

# 4. Upload to PyPI
python -m twine upload dist/*

# 5. Build, tag & push Docker image
./scripts/local-build/build-and-push-final.sh --skip-push
docker tag ghcr.io/julius-muller/vcfcache:latest ghcr.io/julius-muller/vcfcache:v0.3.X
docker push ghcr.io/julius-muller/vcfcache:v0.3.X
docker push ghcr.io/julius-muller/vcfcache:latest

# 6. Tag & create GitHub release
git tag v0.3.X && git push origin v0.3.X
gh release create v0.3.X --title "vcfcache v0.3.X" --notes-file CHANGELOG.md dist/*
```

**For a patch bump without release:**
- Update version numbers, commit and push
- No PyPI upload, no Docker rebuild, no GitHub release

---

## Automated Release (Recommended)

Use the release script for streamlined, interactive releases:

```bash
./scripts/release.sh 0.3.4
```

### Pre-releases (beta/rc)

Use a PEP 440 pre-release version so PyPI/Docker/GitHub stay consistent:
- Beta: `0.4.0b0`, `0.4.0b1`, ...
- Release candidate: `0.4.0rc1`, `0.4.0rc2`, ...

The release script will mark GitHub Releases as pre-releases for `bN/rcN` versions and will avoid pushing the Docker `:latest` tag for pre-releases.
If you need to mark a GitHub release as pre-release without using a PEP 440 pre-release version, run:
```bash
./scripts/release.sh 0.4.0 --github-prerelease
```

The script will walk you through each step with prompts:
- **y** = Yes, proceed with this step
- **n** = No, cancel the release and exit
- **s** = Skip this step and continue

Steps:
1. Update version in pyproject.toml and __init__.py (interactive)
2. Prompt you to update CHANGELOG.md manually
3. Build and test locally - smoke test + full test suite (interactive)
4. Upload to TestPyPI (interactive)
5. Upload to PyPI (interactive)
6. Build and push Docker image (interactive)
7. Create GitHub release (interactive)

You can skip any step (e.g., skip Docker for a quick PyPI-only release).

---

## Prerequisites

- Python 3.11+
- `uv pip install build twine`
- API tokens for https://pypi.org and https://test.pypi.org
- Docker installed + logged in to GHCR: `echo $GITHUB_TOKEN | docker login ghcr.io -u USERNAME --password-stdin`
- GitHub CLI (`gh`) installed

---

## Detailed Release Workflow

### 1. Bump Version & Update Docs

Update version in **3 places**:
```bash
# pyproject.toml
version = "0.3.X"

# vcfcache/__init__.py
__version__ = "0.3.X"

# CHANGELOG.md
## [0.3.X] - 2025-MM-DD
### Added/Changed/Fixed
- ...
```

### 2. Build Package

```bash
rm -rf dist/ build/ *.egg-info
python -m build

# Verify contents
unzip -l dist/vcfcache-*.whl
tar -tzf dist/vcfcache-*.tar.gz
```

### 3. Upload to TestPyPI

```bash
python -m twine upload --repository testpypi dist/*
```

### 4. Test from TestPyPI

```bash
uv venv /tmp/test-env
source /tmp/test-env/bin/activate

uv pip install --python /tmp/test-env/bin/python \
  --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple/ \
  'vcfcache[dev]'

# Run smoke test
vcfcache demo --smoke-test

# Run test suite
python -m pytest tests -q

# Verify
vcfcache --version
vcfcache --help

# Cleanup
deactivate && rm -rf /tmp/test-env
```

**Expected results:**
- ✓ Demo completes all 4 steps
- ✓ All tests pass
- ✓ Version matches release

### 5. Upload to PyPI

```bash
python -m twine upload dist/*
```

### 6. Verify from PyPI

```bash
uv venv /tmp/pypi-test && source /tmp/pypi-test/bin/activate
pip install 'vcfcache[dev]'
vcfcache demo --smoke-test
python -m pytest tests -q
deactivate && rm -rf /tmp/pypi-test
```

### 7. Build & Push Docker Image

**Use the build script:**
```bash
# Build and test (but don't push yet)
./scripts/local-build/build-and-push-final.sh --skip-push
```

This script:
- Builds `docker/Dockerfile.vcfcache`
- Tags as `ghcr.io/julius-muller/vcfcache:latest`
- Runs tests inside the container
- Smoke tests with `vcfcache --version`

**Tag with version and push:**
```bash
# Tag with release version
docker tag ghcr.io/julius-muller/vcfcache:latest \
  ghcr.io/julius-muller/vcfcache:v0.3.X

# Push both tags
docker push ghcr.io/julius-muller/vcfcache:v0.3.X
docker push ghcr.io/julius-muller/vcfcache:latest
```

**Script options:**
- `--skip-tests`: Skip container tests (not recommended)
- `--skip-push`: Don't push after building (use for manual tagging)
- `--force`: Force rebuild even if image exists

### 8. Tag & Create GitHub Release

```bash
# Create and push tag
git tag v0.3.X
git push origin v0.3.X

# Create GitHub release with artifacts
gh release create v0.3.X \
  --title "vcfcache v0.3.X" \
  --notes-file CHANGELOG.md \
  dist/*
```

This creates a release with:
- Release notes from CHANGELOG.md
- Wheel and source distribution as downloadable artifacts

---

## Optional: Blueprint & Annotated Docker Images

**These are only needed for major releases with significant changes.**

### Blueprint Image (lightweight demo, no VEP)

```bash
./scripts/local-build/03-build-blueprint.sh \
  path/to/gnomad_subset.bcf \
  --push

# Tags and pushes to:
# ghcr.io/julius-muller/vcfcache-blueprint:gnomad-grch38-joint-chry-af010
```

### Annotated Image (full VEP + pre-annotated cache)

```bash
# Step 1: Build base image with VEP
./scripts/local-build/04a-build-base-image.sh \
  path/to/gnomad.bcf \
  --yes

# Step 2: Run annotation and commit (requires VEP cache mounted)
./scripts/local-build/04b-annotate-and-commit.sh \
  --base-image <image-from-step-1> \
  --vep-cache-dir /path/to/vep/cache \
  --yes
```

**When to rebuild these:**
- Major version changes (1.0.0, 2.0.0)
- Significant workflow changes affecting caching
- Updated gnomAD version
- **Skip for patch releases**

---

## Troubleshooting

**Demo fails: "bcftools not found"**
- Install bcftools >= 1.20: `apt install bcftools` or `brew install bcftools`

**Tests fail: "No module named pytest"**
- Install dev dependencies: `pip install 'vcfcache[dev]'`

**Upload fails: "File already exists"**
- Cannot overwrite PyPI versions. Bump version, rebuild, re-upload.

**Docker push fails: "authentication required"**
- Login: `echo $GITHUB_TOKEN | docker login ghcr.io -u USERNAME --password-stdin`
- Token needs `write:packages` scope

**Docker build fails**
- Try with `--force` flag: `./scripts/local-build/build-and-push-final.sh --force --skip-push`
- Check Dockerfile syntax: `docker build -f docker/Dockerfile.vcfcache .`
