#!/usr/bin/env bash
set -euo pipefail

# Release automation script for vcfcache
# Usage: ./scripts/release.sh <version>

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Parse arguments
VERSION=""
FORCE_GH_PRERELEASE=0

show_help() {
  cat << EOF
Usage: $0 <version>

Release automation script for vcfcache.
After each step, you'll be prompted to continue, skip, or cancel.

Arguments:
  <version>        Version to release (e.g., 0.4.0, 0.4.0b0, 0.4.0rc1)

Options:
  -h, --help       Show this help message
  --github-prerelease  Mark the GitHub Release as pre-release (even if version is not a PEP 440 pre-release)

Examples:
  $0 0.3.4
  $0 0.4.0b0
  $0 0.4.0 --github-prerelease

Interactive prompts:
  y - Yes, proceed with this step
  n - No, cancel the release and exit
  s - Skip this step and continue
EOF
}

while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
      show_help
      exit 0
      ;;
    --github-prerelease)
      FORCE_GH_PRERELEASE=1
      shift
      ;;
    -*)
      echo "Unknown option: $1"
      show_help
      exit 1
      ;;
    *)
      if [[ -z "$VERSION" ]]; then
        VERSION="$1"
      else
        echo "Error: Unexpected argument: $1"
        show_help
        exit 1
      fi
      shift
      ;;
  esac
done

if [[ -z "$VERSION" ]]; then
  echo "Error: Version is required"
  show_help
  exit 1
fi

log() { echo "[$(date +%H:%M:%S)] $*"; }

pip_install() {
  # Prefer uv for speed/caching; fall back to pip.
  local pybin
  pybin="$(command -v python || true)"
  if command -v uv >/dev/null 2>&1; then
    uv pip install --python "$pybin" "$@" >/dev/null
  else
    python -m pip install -U "$@" >/dev/null
  fi
}

pyproject_build_requires() {
  # Print PEP517 build-system.requires entries (space-separated), or empty on failure.
  python - <<'PY' 2>/dev/null || true
import tomllib
from pathlib import Path
data = tomllib.loads(Path("pyproject.toml").read_text(encoding="utf-8"))
reqs = data.get("build-system", {}).get("requires", [])
print(" ".join(reqs))
PY
}

is_pep440_prerelease() {
  # Prefer a real PEP 440 check (packaging); fall back to a heuristic.
  local v="$1"
  if python -c "from packaging.version import Version; import sys; print('1' if Version(sys.argv[1]).is_prerelease else '0')" "$v" >/dev/null 2>&1; then
    python -c "from packaging.version import Version; import sys; print('1' if Version(sys.argv[1]).is_prerelease else '0')" "$v"
    return 0
  fi
  if [[ "$v" =~ (a|b|rc)[0-9]+$ ]]; then
    echo "1"
  else
    echo "0"
  fi
}

ask_yn_skip() {
  local prompt="$1"
  while true; do
    read -p "$prompt (y=yes, n=cancel, s=skip) " -n 1 -r
    echo
    case $REPLY in
      [Yy])
        return 0  # yes
        ;;
      [Nn])
        log "Cancelled by user"
        exit 1
        ;;
      [Ss])
        return 1  # skip
        ;;
      *)
        echo "Invalid input. Please enter y (yes), n (cancel), or s (skip)"
        ;;
    esac
  done
}

ensure_clean_git() {
  if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    return 0
  fi
  if [[ -n "$(git status --porcelain)" ]]; then
    log "  ✗ Git working tree is not clean."
    git status --porcelain
    log "    Commit/stash changes or run the release from a clean state."
    return 1
  fi
  return 0
}

commit_version_bump() {
  # Commit version-related file changes if present.
  # Only stages known files to avoid sweeping unrelated changes.
  local files=()
  [[ -f pyproject.toml ]] && files+=("pyproject.toml")
  [[ -f vcfcache/__init__.py ]] && files+=("vcfcache/__init__.py")
  [[ -f CHANGELOG.md ]] && files+=("CHANGELOG.md")

  if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    log "  ✗ Not a git repository; cannot commit/tag."
    return 1
  fi

  if [[ ${#files[@]} -eq 0 ]]; then
    log "  ✗ No version files found to commit."
    return 1
  fi

  # Stage only if they have changes.
  local any=false
  for f in "${files[@]}"; do
    if ! git diff --quiet -- "$f" >/dev/null 2>&1; then
      git add "$f"
      any=true
    fi
  done

  if ! $any; then
    log "  ✓ No version file changes to commit"
    return 0
  fi

  git commit -m "Bump version to $VERSION" >/dev/null
  log "  ✓ Committed version bump"
  return 0
}

push_current_branch() {
  local branch
  branch="$(git rev-parse --abbrev-ref HEAD)"
  if [[ -z "$branch" ]] || [[ "$branch" == "HEAD" ]]; then
    log "  ✗ Detached HEAD; cannot push branch."
    return 1
  fi
  git push origin "$branch" >/dev/null
  log "  ✓ Pushed branch to origin ($branch)"
  return 0
}

create_and_push_tag() {
  local tag="v$VERSION"
  if git rev-parse "$tag" >/dev/null 2>&1; then
    log "  ✓ Git tag $tag already exists locally"
  else
    git tag -a "$tag" -m "Release $tag"
    log "  ✓ Created git tag $tag"
  fi
  git push origin "$tag" >/dev/null
  log "  ✓ Pushed git tag $tag"
  return 0
}

docker_login_ghcr() {
  # Non-interactive GHCR login using a token.
  # Prefer GHCR_TOKEN (PAT with write:packages). Fallback to gh auth token if available.
  local registry="ghcr.io"

  local user="${GHCR_USERNAME:-}"
  if [[ -z "$user" ]] && command -v gh >/dev/null 2>&1; then
    user="$(gh api user -q .login 2>/dev/null || true)"
  fi
  if [[ -z "$user" ]]; then
    log "  ✗ Could not determine GH username. Set GHCR_USERNAME (e.g. julius-muller)."
    return 1
  fi

  local token="${GHCR_TOKEN:-}"
  if [[ -z "$token" ]] && command -v gh >/dev/null 2>&1; then
    token="$(gh auth token 2>/dev/null || true)"
  fi
  if [[ -z "$token" ]]; then
    log "  ✗ No GHCR token available."
    log "    Export GHCR_TOKEN (PAT with write:packages) then run this script again, or run:"
    log "      echo \"\$GHCR_TOKEN\" | docker login ghcr.io -u $user --password-stdin"
    return 1
  fi

  echo "$token" | docker login "$registry" -u "$user" --password-stdin >/dev/null
  return 0
}

derive_wiki_repo_url() {
  # Derive https URL for the GitHub wiki repo from origin remote.
  # Returns: prints URL or empty string.
  # Can be overridden with WIKI_REPO_URL env var.
  if [[ -n "${WIKI_REPO_URL:-}" ]]; then
    echo "$WIKI_REPO_URL"
    return 0
  fi

  local origin
  origin="$(git remote get-url origin 2>/dev/null || true)"
  if [[ -z "$origin" ]]; then
    echo ""
    return 0
  fi

  # Normalize origin to owner/repo
  local slug=""
  if [[ "$origin" =~ ^git@github\.com:(.+)\.git$ ]]; then
    slug="${BASH_REMATCH[1]}"
  elif [[ "$origin" =~ ^https?://github\.com/(.+)\.git$ ]]; then
    slug="${BASH_REMATCH[1]}"
  elif [[ "$origin" =~ ^https?://github\.com/(.+)$ ]]; then
    slug="${BASH_REMATCH[1]}"
  fi

  if [[ -z "$slug" ]]; then
    echo ""
    return 0
  fi

  echo "https://github.com/${slug}.wiki.git"
}

sync_github_wiki_from_file() {
  # Sync ./WIKI.md into the GitHub wiki repo as Home.md.
  # Requires git push permissions to <repo>.wiki.git
  local wiki_url
  wiki_url="$(derive_wiki_repo_url)"
  if [[ -z "$wiki_url" ]]; then
    log "  ✗ Could not derive wiki repo URL from origin; set WIKI_REPO_URL manually."
    return 1
  fi
  if [[ ! -f "WIKI.md" ]]; then
    log "  ✗ WIKI.md not found in repo root."
    return 1
  fi

  local tmp
  tmp="$(mktemp -d -t vcfcache-wiki.XXXXXX)"
  log "  → Cloning wiki repo: $wiki_url"
  if ! git clone --depth 1 "$wiki_url" "$tmp" >/dev/null 2>&1; then
    log "  ✗ Failed to clone wiki repo."
    log "    If the wiki is disabled, enable it in GitHub repo settings (Features → Wikis), or skip this step."
    log "    If auth fails, ensure you can push to the wiki repo: $wiki_url"
    rm -rf "$tmp"
    return 1
  fi

  cp WIKI.md "$tmp/Home.md"
  (
    cd "$tmp"
    if git diff --quiet -- Home.md; then
      log "  ✓ Wiki already up-to-date (no changes)"
      return 0
    fi
    git add Home.md
    git commit -m "Sync wiki from WIKI.md (v$VERSION)" >/dev/null
    git push >/dev/null
    log "  ✓ Updated GitHub Wiki Home.md"
  )
  rm -rf "$tmp"
  return 0
}

cd "$PROJECT_ROOT"

# Require a clean working tree for predictable release tags.
if ! ensure_clean_git; then
  exit 1
fi

# Get current versions from files
CURRENT_VERSION_PYPROJECT=$(grep '^version = ' pyproject.toml | sed 's/version = "\(.*\)"/\1/')
CURRENT_VERSION_INIT=$(grep '^__version__ = ' vcfcache/__init__.py | sed 's/__version__ = "\(.*\)"/\1/')
IS_PRERELEASE="$(is_pep440_prerelease "$VERSION")"
GH_PRERELEASE_FLAG=""
if [[ "$IS_PRERELEASE" == "1" ]] || [[ "$FORCE_GH_PRERELEASE" == "1" ]]; then
  GH_PRERELEASE_FLAG="--prerelease"
fi

release_notes_file() {
  # Create a temp release notes file containing only the relevant CHANGELOG section.
  # Falls back to CHANGELOG.md if extraction fails.
  local tmp
  tmp="$(mktemp -t "vcfcache-release-notes.XXXXXX.md")"
  if .venv/bin/python scripts/release_notes_from_changelog.py "$VERSION" >"$tmp" 2>/dev/null; then
    echo "$tmp"
    return 0
  fi
  if python scripts/release_notes_from_changelog.py "$VERSION" >"$tmp" 2>/dev/null; then
    echo "$tmp"
    return 0
  fi
  rm -f "$tmp"
  echo "CHANGELOG.md"
}

log "Starting release process for version $VERSION"
log "Current versions: pyproject.toml=$CURRENT_VERSION_PYPROJECT, __init__.py=$CURRENT_VERSION_INIT"
if [[ "$IS_PRERELEASE" == "1" ]]; then
  log "Detected pre-release version (PEP 440): $VERSION"
fi
if [[ "$FORCE_GH_PRERELEASE" == "1" ]] && [[ "$IS_PRERELEASE" != "1" ]]; then
  log "GitHub Release will be marked as pre-release (--github-prerelease)"
fi
echo ""

# Step 1: Update version + commit
log "Step 1: Update version (and commit)"

# Check if version is already updated
if [[ "$CURRENT_VERSION_PYPROJECT" == "$VERSION" ]] && [[ "$CURRENT_VERSION_INIT" == "$VERSION" ]]; then
  log "  ✓ Version already set to $VERSION, skipping"
else
  if ask_yn_skip "Proceed with version update? ($CURRENT_VERSION_PYPROJECT -> $VERSION)"; then

    # Update pyproject.toml
    if [[ -f "pyproject.toml" ]]; then
      sed -i "s/^version = .*/version = \"$VERSION\"/" pyproject.toml
      log "  ✓ Updated pyproject.toml"
    else
      log "  ✗ pyproject.toml not found"
      exit 1
    fi

    # Update __init__.py
    if [[ -f "vcfcache/__init__.py" ]]; then
      sed -i "s/^__version__ = .*/__version__ = \"$VERSION\"/" vcfcache/__init__.py
      log "  ✓ Updated vcfcache/__init__.py"
    else
      log "  ✗ vcfcache/__init__.py not found"
      exit 1
    fi

    # Reminder to update CHANGELOG.md manually
    log "  ⚠ Please manually update CHANGELOG.md with:"
    log "    ## $VERSION ($(date +%Y-%m-%d))"
    log "    ### Added/Changed/Fixed"
    log "    - ..."
    echo ""
    read -p "Press Enter once CHANGELOG.md is updated..."
    log "  ✓ Version update complete"
  else
    log "  ⊘ Skipped version update"
  fi
fi
echo ""

# Step 1b: Commit + push version bump
log "Step 1b: Commit & push version bump"
if ask_yn_skip "Commit and push version bump to origin?"; then
  commit_version_bump
  push_current_branch
else
  log "  ⊘ Skipped commit/push"
fi
echo ""

# Step 2: Build and test locally
log "Step 2: Build and test locally"

# Check if package is already built
if [[ -f "dist/vcfcache-${VERSION}-py3-none-any.whl" ]] && [[ -f "dist/vcfcache-${VERSION}.tar.gz" ]]; then
  log "  ✓ Package v$VERSION already built, skipping"
else
  if ask_yn_skip "Build package v$VERSION and run tests?"; then
    rm -rf dist/ build/ *.egg-info
    log "  → Building distributions..."
    # Prefer a fast, non-isolated build (avoids pip downloading build backends in an isolated env).
    # Fall back to isolated build if the local environment is missing build requirements.
    pip_install build
    build_reqs="$(pyproject_build_requires)"
    if [[ -n "$build_reqs" ]]; then
      # shellcheck disable=SC2086
      pip_install $build_reqs
    else
      pip_install hatchling editables
    fi

    if python -m build --no-isolation; then
      log "  ✓ Built package (no isolation)"
    else
      log "  ⚠ No-isolation build failed; retrying with isolated build (may download build deps)..."
      if command -v uv >/dev/null 2>&1; then
        python -m build --installer uv
      else
        python -m build
      fi
      log "  ✓ Built package (isolated)"
    fi

    log "  → Installing in temporary venv and running tests..."
    uv venv /tmp/vcfcache-release-test
    source /tmp/vcfcache-release-test/bin/activate
    uv pip install "dist/vcfcache-${VERSION}-py3-none-any.whl[dev]"

    log "  → Running smoke test..."
    vcfcache demo --smoke-test --quiet

    log "  → Running test suite..."
    python -m pytest tests -q

    deactivate
    rm -rf /tmp/vcfcache-release-test

    log "  ✓ Local tests passed"
  else
    log "  ⊘ Skipped build and test"
  fi
fi
echo ""

# Step 3: Upload to TestPyPI
log "Step 3: Upload to TestPyPI"

if [[ "$IS_PRERELEASE" == "1" ]]; then
  log "  → Pre-release detected: TestPyPI is recommended for betas/RCs"
fi
if ask_yn_skip "Upload v$VERSION to TestPyPI?"; then
  python -m twine upload --repository testpypi dist/*
  log "  ✓ Uploaded to TestPyPI"
  log "  → Verify at: https://test.pypi.org/project/vcfcache/$VERSION/"
  echo ""
  read -p "Press Enter once TestPyPI verification is complete..."
else
  log "  ⊘ Skipped TestPyPI upload"
fi
echo ""

# Step 4: Upload to PyPI
log "Step 4: Upload to PyPI"

if [[ "$IS_PRERELEASE" == "1" ]]; then
  log "  ⚠ Pre-release detected: skipping PyPI upload is usually recommended"
fi
if ask_yn_skip "Upload v$VERSION to PyPI?"; then
  python -m twine upload dist/*
  log "  ✓ Uploaded to PyPI"
  log "  → Live at: https://pypi.org/project/vcfcache/$VERSION/"
else
  log "  ⊘ Skipped PyPI upload"
fi
echo ""

# Step 5: Docker build and push
log "Step 5: Build and push Docker image"

# Ask user if they want to build/push Docker image
if ! ask_yn_skip "Build Docker image for v$VERSION?"; then
  log "  ⊘ Skipped Docker build"
else
  # Build image
  log "  → Building Docker image..."
  ./scripts/local-build/build-and-push-final.sh --skip-push --force

  # Tag with version
  log "  → Tagging as v$VERSION..."
  docker tag ghcr.io/julius-muller/vcfcache:latest "ghcr.io/julius-muller/vcfcache:v$VERSION"

  # Push
  if ask_yn_skip "Push Docker images to GHCR?"; then
    log "  → Logging in to GHCR (ghcr.io)..."
    if ! docker_login_ghcr; then
      log "  ✗ GHCR login failed."
      log "    If you previously logged into Docker Hub, make sure to login specifically to GHCR:"
      log "      docker logout ghcr.io || true"
      log "      echo \"\$GHCR_TOKEN\" | docker login ghcr.io -u <github-username> --password-stdin"
      log "    Token must include: write:packages (and read:packages)."
      exit 1
    fi

    docker push "ghcr.io/julius-muller/vcfcache:v$VERSION"
    if [[ "$IS_PRERELEASE" != "1" ]]; then
      docker push ghcr.io/julius-muller/vcfcache:latest
    else
      log "  → Pre-release detected: not pushing :latest"
    fi
    log "  ✓ Pushed Docker images"

    # Verify
    log "  → Verifying pushed image..."
    docker pull "ghcr.io/julius-muller/vcfcache:v$VERSION" --quiet
    docker run --rm "ghcr.io/julius-muller/vcfcache:v$VERSION" demo --smoke-test --quiet
    log "  ✓ Docker image verified"
  else
    log "  ⊘ Skipped Docker push (images tagged locally)"
  fi
fi
echo ""

# Step 6: Tag + GitHub release
log "Step 6: Tag and create GitHub release"

# Check if release already exists
if gh release view "v$VERSION" >/dev/null 2>&1; then
  log "  ✓ GitHub release v$VERSION already exists, skipping"
  exit 0
fi

# Ensure tag exists and is pushed (tag should always point at the version bump commit)
if ask_yn_skip "Create and push git tag v$VERSION now?"; then
  create_and_push_tag
else
  log "  ⊘ Skipped tag creation"
fi

# Ask user if they want to create the release
if ! ask_yn_skip "Create GitHub release v$VERSION with artifacts?"; then
  log "  ⊘ Skipped GitHub release creation"
  exit 0
fi

if ! git rev-parse "v$VERSION" >/dev/null 2>&1; then
  log "  ✗ Git tag v$VERSION does not exist. Create it first."
  exit 1
fi

# Create GitHub release
log "  → Creating GitHub release..."
NOTES_FILE="$(release_notes_file)"
gh release create "v$VERSION" \
  --title "vcfcache v$VERSION" \
  $GH_PRERELEASE_FLAG \
  --notes-file "$NOTES_FILE" \
  dist/*
log "  ✓ GitHub release created"
log "  → View at: https://github.com/julius-muller/vcfcache/releases/tag/v$VERSION"

echo ""

# Step 7: Sync GitHub Wiki (optional)
log "Step 7: Sync GitHub Wiki (optional)"
log "  → This updates the GitHub Wiki homepage (Home.md) from repo WIKI.md."
if ask_yn_skip "Sync GitHub Wiki from WIKI.md now?"; then
  sync_github_wiki_from_file || log "  ⚠ Wiki sync failed (non-fatal)."
else
  log "  ⊘ Skipped wiki sync"
fi

echo ""
log "================================================"
log "Release v$VERSION complete! 🎉"
log "================================================"
