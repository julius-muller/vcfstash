#!/usr/bin/env bash
set -euo pipefail

#########################
# Configurable settings #
#########################

# bcftools version (can be overridden via env: BCFTOOLS_VERSION=1.22 ./update-bcftools.sh)
BCFTOOLS_VERSION="${BCFTOOLS_VERSION:-1.22}"

# Download base: GitHub releases
# Final URL:
#   https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2
BCFTOOLS_URL_BASE="${BCFTOOLS_URL_BASE:-https://github.com/samtools/bcftools/releases/download}"

# Number of make jobs (override via MAKE_JOBS)
MAKE_JOBS="${MAKE_JOBS:-$(nproc)}"

#########################
# Derived paths         #
#########################

# Project root = directory containing this script, then go up one (assuming scripts/…)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

TOOLS_DIR="$PROJECT_ROOT/tools"
TOOLS_LIB_DIR="$TOOLS_DIR/lib"

BCFTOOLS_WRAPPER="$TOOLS_DIR/bcftools"
BCFTOOLS_REAL="$TOOLS_DIR/bcftools.real"

#########################
# Sanity checks         #
#########################

if [[ ! -d "$TOOLS_DIR" ]]; then
    echo "ERROR: tools/ directory not found under $PROJECT_ROOT" >&2
    exit 1
fi

mkdir -p "$TOOLS_LIB_DIR"

if [[ ! -x "$BCFTOOLS_WRAPPER" ]]; then
    echo "WARNING: $BCFTOOLS_WRAPPER is not executable."
    echo "         Make sure you have created the wrapper script as discussed."
fi

#########################
# Build in temp dir     #
#########################

TMPDIR="$(mktemp -d -t bcftools-build-XXXXXX)"
trap 'rm -rf "$TMPDIR"' EXIT

echo "Using temp build dir: $TMPDIR"

TARBALL="bcftools-${BCFTOOLS_VERSION}.tar.bz2"
TARBALL_URL="${BCFTOOLS_URL_BASE}/${BCFTOOLS_VERSION}/${TARBALL}"

echo "Downloading bcftools ${BCFTOOLS_VERSION} from:"
echo "  ${TARBALL_URL}"

cd "$TMPDIR"

# Prefer curl, fall back to wget
if command -v curl >/dev/null 2>&1; then
    curl -L -o "$TARBALL" "$TARBALL_URL"
elif command -v wget >/dev/null 2>&1; then
    wget -O "$TARBALL" "$TARBALL_URL"
else
    echo "ERROR: Neither curl nor wget found; please install one of them." >&2
    exit 1
fi

# Quick sanity check: should be bzip2 (first 3 bytes "BZh")
if ! head -c 3 "$TARBALL" | grep -q 'BZh'; then
    echo "ERROR: $TARBALL does not look like a bzip2 file (wrong URL / HTML redirect?)." >&2
    file "$TARBALL" || true
    exit 1
fi

echo "Extracting ${TARBALL}..."
tar xjf "$TARBALL"

SRC_DIR="$TMPDIR/bcftools-${BCFTOOLS_VERSION}"
if [[ ! -d "$SRC_DIR" ]]; then
    echo "ERROR: Expected source dir $SRC_DIR not found after extraction." >&2
    exit 1
fi

cd "$SRC_DIR"

PREFIX="$TMPDIR/prefix"
mkdir -p "$PREFIX"

echo "Configuring bcftools with prefix: $PREFIX"
./configure --prefix="$PREFIX"

echo "Building bcftools (make -j${MAKE_JOBS})..."
make -j"${MAKE_JOBS}"

echo "Installing into prefix..."
make install

#########################
# Copy bcftools.real    #
#########################

if [[ ! -x "$PREFIX/bin/bcftools" ]]; then
    echo "ERROR: Installed bcftools not found at $PREFIX/bin/bcftools" >&2
    exit 1
fi

echo "Updating $BCFTOOLS_REAL..."
cp "$PREFIX/bin/bcftools" "$BCFTOOLS_REAL"
chmod +x "$BCFTOOLS_REAL"

#########################
# Copy dependent libs   #
#########################

# Helper: copy one dependency found by ldd (if present)
copy_dep() {
    local pattern="$1"
    local name="$2"

    local line
    line="$(ldd "$BCFTOOLS_REAL" | grep "$pattern" || true)"

    if [[ -z "$line" ]]; then
        echo "NOTE: $name ($pattern) not found in ldd output (maybe linked statically)."
        return 0
    fi

    # Try to extract path after '=>'
    local path
    path="$(awk -v pat="$pattern" '
        $0 ~ pat {
            for (i=1; i<=NF; i++) {
                if ($i == "=>") { print $(i+1); exit }
            }
        }
    ' <<< "$line")"

    # Some systems print "libX.so => not found" so guard that
    if [[ -z "$path" || "$path" == "not" ]]; then
        echo "WARNING: Could not extract path for $name from: $line" >&2
        return 1
    fi

    echo "Copying $name from $path -> $TOOLS_LIB_DIR/"
    cp "$path" "$TOOLS_LIB_DIR"/
}

echo "Cleaning old bundled libs..."
rm -f "$TOOLS_LIB_DIR"/libhts.so.* \
      "$TOOLS_LIB_DIR"/libdeflate.so.* \
      "$TOOLS_LIB_DIR"/libhtscodecs.so.*

echo "Inspecting dynamic dependencies with ldd on $BCFTOOLS_REAL..."
ldd "$BCFTOOLS_REAL" || true

# Try to copy the usual suspects (if they are dynamically linked)
copy_dep "libhts\.so"      "libhts"
copy_dep "libdeflate\.so"  "libdeflate"
copy_dep "libhtscodecs\.so" "libhtscodecs"

#########################
# Final check           #
#########################

echo "Verifying that bcftools (via wrapper) runs..."
set +e
"$BCFTOOLS_WRAPPER" --version || {
    echo "ERROR: bcftools wrapper failed to run." >&2
}
set -e

echo "Done. bcftools ${BCFTOOLS_VERSION} has been built and installed into:"
echo "  $BCFTOOLS_REAL"
echo "Any dynamic libs found for hts/deflate/htscodecs have been copied into:"
echo "  $TOOLS_LIB_DIR"

