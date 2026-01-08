#!/usr/bin/env bash
set -euo pipefail

: "${ZENODO_TOKEN:?Please export ZENODO_TOKEN (prod token)}"
ZENODO_SANDBOX=${ZENODO_SANDBOX:-0}

ALIAS="GRCh38-af010-vep115.2_basic"
work=/tmp/vcfcache_dummy_upload
rm -rf "$work"
mkdir -p "$work"

CACHE_DIR="$work/cache_$ALIAS"
TAR_PATH="$work/${ALIAS}.tar.gz"
DL_DIR="$work/download"

echo "[1/6] Create dummy cache"
python - "$CACHE_DIR" "$ALIAS" <<'PY'
import os, sys
from pathlib import Path
from tests.test_cli_alias_and_pull import make_dummy_cache

cache_dir = Path(sys.argv[1])
alias = sys.argv[2]
tmp = cache_dir.parent
tmp.mkdir(parents=True, exist_ok=True)
cache = make_dummy_cache(tmp, alias)
print(cache)
PY

echo "[2/6] Tar cache"
python - "$CACHE_DIR" "$TAR_PATH" <<'PY'
import sys
from pathlib import Path
from vcfcache.utils.archive import tar_cache
tar_cache(Path(sys.argv[1]), Path(sys.argv[2]))
PY

echo "[3/6] Upload to Zenodo (prod)"
DOI=$(python - "$ZENODO_TOKEN" "$ZENODO_SANDBOX" "$TAR_PATH" "$ALIAS" <<'PY'
import os, sys, requests
from pathlib import Path
from vcfcache.integrations import zenodo

token = sys.argv[1]
sandbox = sys.argv[2] == "1"
tar_path = Path(sys.argv[3])
alias = sys.argv[4]

blueprints = {
    "title": f"VCFcache dummy cache {alias}",
    "upload_type": "dataset",
    "description": "Dummy cache generated for automated testing of VCFcache Zenodo upload pipeline.",
    "creators": [{"name": "VCFcache Bot"}],
}

api_base = zenodo.ZENODO_SANDBOX_API if sandbox else zenodo.ZENODO_API

dep = zenodo.create_deposit(token, sandbox=sandbox)
requests.put(
    f"{api_base}/deposit/depositions/{dep['id']}",
    params={"access_token": token},
    json={"metadata": blueprints},
    timeout=30,
).raise_for_status()
zenodo.upload_file(dep, tar_path, token, sandbox=sandbox)
if not sandbox:
    dep = zenodo.publish_deposit(dep, token, sandbox=sandbox)
print(dep.get("doi", "draft"))
PY
)
echo "DOI: $DOI"

md5=$(python - "$TAR_PATH" <<'PY'
import sys
from pathlib import Path
from vcfcache.utils.archive import file_md5
print(file_md5(Path(sys.argv[1])))
PY
)

echo "[4/6] Download to verify"
mkdir -p "$DL_DIR"
python - "$DOI" "$DL_DIR" <<'PY'
import sys
from pathlib import Path
from vcfcache.integrations.zenodo import download_doi
from vcfcache.utils.archive import extract_cache

doi = sys.argv[1]
dest_dir = Path(sys.argv[2])
tar_dest = dest_dir / "cache.tar.gz"
download_doi(doi, tar_dest, sandbox=sandbox)
extracted = extract_cache(tar_dest, dest_dir)
print("Extracted:", extracted)
PY

echo "[5/6] Manifest entry to add:"
cat <<EOF
- alias: $ALIAS
  doi: $DOI
  version: "0.3.0"
  genome: GRCh38
  af: "0.10"
  tool: vep115.2
  updated_at: $(date +%Y-%m-%d)
  md5: $md5
  annotation_yaml_md5: placeholder
EOF

echo "[6/6] Done. Work dir: $work"
