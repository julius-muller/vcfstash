from pathlib import Path

from vcfcache.manifest import find_alias, load_manifest


def test_manifest_load_and_find_alias(tmp_path):
    manifest_file = tmp_path / "m.yaml"
    manifest_file.write_text(
        """
- alias: GRCh38-af010-vep115.2_basic
  doi: 10.5281/zenodo.12345
  genome: GRCh38
  af: "0.10"
  tool: vep115.2
  image_tag: vcfcache:vep115.2_basic
  updated_at: 2025-01-01
  md5: abc
"""
    )

    entries = load_manifest(str(manifest_file))
    hit = find_alias(entries, "GRCh38-af010-vep115.2_basic")
    assert hit is not None
    assert hit.doi == "10.5281/zenodo.12345"

