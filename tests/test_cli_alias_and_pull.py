import json
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest
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
    (workflow_dir / "init.yaml").write_text(
        "annotation_tool_cmd: echo\n"
        "bcftools_cmd: echo\n"
        "temp_dir: /tmp\n"
        "threads: 1\n"
        "genome_build: GRCh38\n"
    )
    (cache_dir / "annotation.yaml").write_text(
        "annotation_cmd: echo annotate\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )
    return cache_root


@pytest.mark.skip(
    reason="Annotate command no longer resolves Zenodo aliases. "
    "Use 'vcfcache cache-build --doi' to download caches instead."
)
def test_cli_annotate_alias_resolves_and_prints_command(tmp_path, monkeypatch, capsys):
    alias = "cache-hg38-gnomad-4.1wgs-AF0100-vep-115.2-basic"
    cache_root = make_dummy_cache(tmp_path, alias)
    tar_path = tmp_path / "dummy.tar.gz"
    tar_cache(cache_root, tar_path)

    def fake_download(doi, dest, sandbox=False):
        dest.write_bytes(tar_path.read_bytes())
        return dest

    monkeypatch.setattr(cli, "download_doi", fake_download)
    monkeypatch.setattr(cli, "resolve_zenodo_alias", lambda a, item_type, sandbox: ("10.5281/zenodo.fake", alias))

    args = [
        "annotate",
        "-a",
        alias,
        "--vcf",
        str(cache_root / "blueprint" / "dummy.bcf"),  # not used when --requirements
        "--output",
        str(tmp_path / "out"),
        "--requirements",
    ]

    monkeypatch.setattr(cli.Path, "home", lambda: tmp_path)
    monkeypatch.setattr(cli, "sys", sys)
    monkeypatch.setattr(sys, "argv", ["vcfcache"] + args)
    cli.main()

    captured = capsys.readouterr()
    assert "echo annotate" in captured.out


def test_cli_list_queries_zenodo(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(
        cli,
        "search_zenodo_records",
        lambda item_type, genome=None, source=None, sandbox=False, min_size_mb=1.0: [
            {"title": "VCFcache cache: x", "doi": "10.5281/zenodo.fake", "created": "2025-01-01", "size_mb": 0.01}
        ],
    )

    args = ["list", "caches"]
    monkeypatch.setattr(cli, "sys", sys)
    monkeypatch.setattr(sys, "argv", ["vcfcache"] + args)
    cli.main()
    captured = capsys.readouterr()
    assert "Available vcfcache caches on Zenodo" in captured.out
    assert "10.5281/zenodo.fake" in captured.out
