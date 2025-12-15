# Releases (maintainers)

This file is for maintainers/dev workflows. User docs live in `README.md` and `WIKI.md`.

## Recommended model

- Develop on `main` (keep it green).
- Cut releases from immutable git tags.
- Keep PyPI + Docker + GitHub Releases aligned to the tag commit.

## Versioning (beta/rc/final)

Use PEP 440 versions:
- Beta: `0.4.0b0`, `0.4.0b1`, ...
- Release candidate: `0.4.0rc1`, `0.4.0rc2`, ...
- Final: `0.4.0`

Tag naming:
- `v0.4.0b0` / `v0.4.0rc1` / `v0.4.0`

## Release workflows

### 1) Automated (recommended): GitHub Actions on tags

Workflow: `.github/workflows/release.yml`

Behavior:
- Trigger: `push` of tags matching `v*`.
- Validates that tag `v<version>` matches `pyproject.toml`’s version.
- Runs tests and builds `dist/*`.
- Creates a GitHub Release (marked as pre-release for `bN/rcN` tags).
- Publishes:
  - pre-releases → TestPyPI
  - finals → PyPI
- Builds/pushes Docker:
  - always pushes `ghcr.io/<owner>/vcfcache:v<version>`
  - pushes `ghcr.io/<owner>/vcfcache:latest` only for final releases

Required secrets:
- `PYPI_API_TOKEN`
- `TEST_PYPI_API_TOKEN`

### 2) Manual (local): `scripts/release.sh`

The script is an interactive checklist that can:
- bump versions in `pyproject.toml` and `vcfcache/__init__.py`
- build + smoke test + run pytest
- upload to (Test)PyPI
- build/tag/push Docker
- create a GitHub Release (pre-release for `bN/rcN`), with optional `--github-prerelease`

See `RELEASE.md` for the detailed step-by-step checklist.

## Practical checklist

1. Update code on `main`; ensure CI green.
2. Pick next version:
   - beta: `0.4.0b0`
   - final: `0.4.0`
3. Update `pyproject.toml`, `vcfcache/__init__.py`, `CHANGELOG.md`.
4. Tag and push:
   - `git tag v0.4.0b0 && git push origin v0.4.0b0`
5. Let `.github/workflows/release.yml` publish artifacts.

