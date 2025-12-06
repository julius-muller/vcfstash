#!/usr/bin/env bash
set -euo pipefail

: "${ZENODO_TOKEN:?Please export ZENODO_TOKEN (prod token)}"
ZENODO_SANDBOX=${ZENODO_SANDBOX:-0}

alias="GRCh38-af010-vep115.2_basic"
work=/tmp/vcfstash_dummy_upload
rm -rf "$work"
mkdir -p "$work"

echo "[1/6] Create dummy cache"
python - <<'PY'
from pathlib import Path
from tests.test_cli_alias_and_pull import make_dummy_cache
tmp = Path("$work")
tmp.mkdir(parents=True, exist_ok=True)
cache = make_dummy_cache(tmp, "$alias")
print(cache)
PY

cache_dir="$work/cache_$alias"

echo "[2/6] Tar cache"
tar_path="$work/${alias}.tar.gz"
python - <<'PY'
from pathlib import Path
from vcfstash.utils.archive import tar_cache
tar_cache(Path("$cache_dir"), Path("$tar_path"))
PY

echo "[3/6] Upload to Zenodo (prod)"
python - <<'PY'
import os
from pathlib import Path
from vcfstash.integrations import zenodo

token = os.environ["ZENODO_TOKEN"]
sandbox = os.environ.get("ZENODO_SANDBOX", "0") == "1"
tar_path = Path("$tar_path")

dep = zenodo.create_deposit(token, sandbox=sandbox)
zenodo.upload_file(dep, tar_path, token, sandbox=sandbox)
if not sandbox:
    dep = zenodo.publish_deposit(dep, token, sandbox=sandbox)
print(dep.get("doi", "draft"), dep.get("id"))
PY

doi=$(python - <<'PY'
import os
from pathlib import Path
from vcfstash.integrations import zenodo

token = os.environ["ZENODO_TOKEN"]
sandbox = os.environ.get("ZENODO_SANDBOX", "0") == "1"
tar_path = Path("$tar_path")

dep = zenodo.create_deposit(token, sandbox=sandbox)
zenodo.upload_file(dep, tar_path, token, sandbox=sandbox)
if not sandbox:
    dep = zenodo.publish_deposit(dep, token, sandbox=sandbox)
print(dep.get("doi", "draft"))
PY)

md5=$(python - <<'PY'
from pathlib import Path
from vcfstash.utils.archive import file_md5
print(file_md5(Path("$tar_path")))
PY)

echo "[4/6] Download to verify"
dl_dir="$work/download"
mkdir -p "$dl_dir"
python - <<'PY'
from pathlib import Path
from vcfstash.integrations.zenodo import download_doi
from vcfstash.utils.archive import extract_cache
doi = "$doi"
tar_dest = Path("$dl_dir") / "cache.tar.gz"
download_doi(doi, tar_dest)
extracted = extract_cache(tar_dest, Path("$dl_dir"))
print("Extracted:", extracted)
PY

echo "[5/6] Manifest entry to add:"
cat <<EOF2
- alias: $alias
  doi: $doi
  genome: GRCh38
  af: "0.10"
  tool: vep115.2
  image_tag: vcfstash:vep115.2_basic
  updated_at: $(date +%Y-%m-%d)
  md5: $md5
  annotation_yaml_md5: placeholder
EOF2

echo "[6/6] Done. Work dir: $work"
