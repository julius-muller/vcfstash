import json
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import sys

import vcfcache.cli as cli
from vcfcache.utils.archive import tar_cache


def make_dummy_cache(tmp_path: Path, alias: str) -> Path:
    cache_root = tmp_path / f"cache_{alias}"
    cache_dir = cache_root / "cache" / alias
    workflow_dir = cache_root / "workflow"
    cache_dir.mkdir(parents=True, exist_ok=True)
    workflow_dir.mkdir(parents=True, exist_ok=True)

    # Minimal required files
    (cache_dir / "vcfcache_annotated.bcf").write_bytes(b"dummy")
    (workflow_dir / "init.yaml").write_text("bcftools_cmd: echo\n")
    (cache_dir / "annotation.yaml").write_text("annotation_tool_cmd: echo annotate\n")
    return cache_root


def test_cli_pull_downloads_and_extracts(tmp_path, monkeypatch, capsys):
    alias = "GRCh38-af010-vep115.2_basic"
    cache_root = make_dummy_cache(tmp_path, alias)
    tar_path = tmp_path / "dummy.tar.gz"
    tar_cache(cache_root, tar_path)

    def fake_download(doi, dest):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(tar_path.read_bytes())
        return dest

    monkeypatch.setattr(cli, "download_doi", fake_download)

    args = [
        "pull",
        "--doi",
        "10.5281/zenodo.fake",
        "--dest",
        str(tmp_path / "out"),
    ]

    monkeypatch.setattr(cli, "sys", sys)
    monkeypatch.setattr(sys, "argv", ["vcfcache"] + args)
    cli.main()

    extracted = tmp_path / "out" / cache_root.name
    assert extracted.exists()
    assert (extracted / "cache" / alias / "vcfcache_annotated.bcf").exists()


def test_cli_annotate_alias_resolves_and_prints_command(tmp_path, monkeypatch, capsys):
    alias = "GRCh38-af010-vep115.2_basic"
    cache_root = make_dummy_cache(tmp_path, alias)
    tar_path = tmp_path / "dummy.tar.gz"
    tar_cache(cache_root, tar_path)

    def fake_download(doi, dest):
        dest.write_bytes(tar_path.read_bytes())
        return dest

    monkeypatch.setattr(cli, "download_doi", fake_download)

    manifest = tmp_path / "manifest.yaml"
    manifest.write_text(
        json.dumps(
            [
                {
                    "alias": alias,
                    "doi": "10.5281/zenodo.fake",
                    "genome": "GRCh38",
                    "af": "0.10",
                    "tool": "vep115.2",
                    "image_tag": "vcfcache:vep115.2_basic",
                    "updated_at": "2025-01-01",
                }
            ]
        )
    )

    args = [
        "annotate",
        "-a",
        alias,
        "--vcf",
        str(cache_root / "blueprint" / "dummy.bcf"),  # not used when --show-command
        "--output",
        str(tmp_path / "out"),
        "--show-command",
        "--manifest",
        str(manifest),
    ]

    monkeypatch.setattr(cli.Path, "home", lambda: tmp_path)
    monkeypatch.setattr(cli, "sys", sys)
    monkeypatch.setattr(sys, "argv", ["vcfcache"] + args)
    cli.main()

    captured = capsys.readouterr()
    assert "echo annotate" in captured.out


def test_cli_list_manifest(tmp_path, monkeypatch, capsys):
    manifest = tmp_path / "manifest.yaml"
    manifest.write_text(
        """
- alias: GRCh38-af010-vep115.2_basic
  doi: 10.5281/zenodo.fake
  genome: GRCh38
  af: "0.10"
  tool: vep115.2
  image_tag: vcfcache:vep115.2_basic
  updated_at: 2025-01-01
"""
    )

    args = ["list", "--public-caches", "--manifest", str(manifest)]
    monkeypatch.setattr(cli, "sys", sys)
    monkeypatch.setattr(sys, "argv", ["vcfcache"] + args)
    cli.main()
    captured = capsys.readouterr()
    assert "alias" in captured.out
    assert "GRCh38-af010-vep115.2_basic" in captured.out
