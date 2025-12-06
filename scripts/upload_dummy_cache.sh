#!/usr/bin/env bash
set -euo pipefail

: "${ZENODO_TOKEN:?Please export ZENODO_TOKEN (prod token)}"
ZENODO_SANDBOX=${ZENODO_SANDBOX:-0}

export ALIAS="GRCh38-af010-vep115.2_basic"
work=/tmp/vcfstash_dummy_upload
rm -rf "$work"
mkdir -p "$work"
export CACHE_DIR="$work/cache_$ALIAS"
export TAR_PATH="$work/${ALIAS}.tar.gz"

echo "[1/6] Create dummy cache"
python - <<'PY'
import os
from pathlib import Path
from tests.test_cli_alias_and_pull import make_dummy_cache

tmp = Path(os.environ["CACHE_DIR"]).parent
tmp.mkdir(parents=True, exist_ok=True)
cache = make_dummy_cache(tmp, os.environ["ALIAS"])
print(cache)
PY

echo "[2/6] Tar cache"
python - <<'PY'
import os
from pathlib import Path
from vcfstash.utils.archive import tar_cache

tar_cache(Path(os.environ["CACHE_DIR"]), Path(os.environ["TAR_PATH"]))
PY

echo "[3/6] Upload to Zenodo (prod)"
DOI=$(python - <<'PY'
import os
from pathlib import Path
import requests
from vcfstash.integrations import zenodo

alias = os.environ["ALIAS"]
metadata = {
    "title": f"VCFstash dummy cache {alias}",
    "upload_type": "dataset",
    "description": "Dummy cache generated for automated testing of VCFstash Zenodo upload pipeline.",
    "creators": [{"name": "VCFstash Bot"}],
}

token = os.environ["ZENODO_TOKEN"]
sandbox = os.environ.get("ZENODO_SANDBOX", "0") == "1"
api_base = zenodo.ZENODO_SANDBOX_API if sandbox else zenodo.ZENODO_API

tar_path = Path(os.environ["TAR_PATH"])

dep = zenodo.create_deposit(token, sandbox=sandbox)
requests.put(
    f"{api_base}/deposit/depositions/{dep['id']}",
    params={"access_token": token},
    json={"metadata": metadata},
    timeout=30,
).raise_for_status()
zenodo.upload_file(dep, tar_path, token, sandbox=sandbox)
if not sandbox:
    dep = zenodo.publish_deposit(dep, token, sandbox=sandbox)
print(dep.get("doi", "draft"))
PY)

echo "DOI: $DOI"

md5=$(python - <<'PY'
from pathlib import Path
from vcfstash.utils.archive import file_md5
print(file_md5(Path("$TAR_PATH")))
PY)

echo "[4/6] Download to verify"
dl_dir="$work/download"
mkdir -p "$dl_dir"
python - <<'PY'
from pathlib import Path
from vcfstash.integrations.zenodo import download_doi
from vcfstash.utils.archive import extract_cache
doi = "$DOI"
tar_dest = Path("$dl_dir") / "cache.tar.gz"
download_doi(doi, tar_dest)
extracted = extract_cache(tar_dest, Path("$dl_dir"))
print("Extracted:", extracted)
PY

echo "[5/6] Manifest entry to add:"
cat <<EOF2
- alias: $ALIAS
  doi: $DOI
  genome: GRCh38
  af: "0.10"
  tool: vep115.2
  image_tag: vcfstash:vep115.2_basic
  updated_at: $(date +%Y-%m-%d)
  md5: $md5
  annotation_yaml_md5: placeholder
EOF2

echo "[6/6] Done. Work dir: $work"
