# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

import os
import subprocess
import sys
from pathlib import Path

import pytest

import vcfcache.cli as cli


def test_load_dotenv_sets_missing(monkeypatch, tmp_path):
    home = tmp_path / "home"
    cwd = tmp_path / "cwd"
    home.mkdir()
    cwd.mkdir()
    (home / ".env").write_text("A=1\nB=2\n")
    (cwd / ".env").write_text("B=3\nC=4\n")

    monkeypatch.setattr(cli.Path, "home", lambda: home)
    monkeypatch.setattr(cli.Path, "cwd", lambda: cwd)
    monkeypatch.delenv("A", raising=False)
    monkeypatch.delenv("B", raising=False)
    monkeypatch.delenv("C", raising=False)

    cli._load_dotenv()
    assert os.environ["A"] == "1"
    assert os.environ["B"] == "2"  # cwd does not override existing values
    assert os.environ["C"] == "4"


def test_load_dotenv_preserves_existing(monkeypatch, tmp_path):
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    (cwd / ".env").write_text("A=1\n")
    monkeypatch.setattr(cli.Path, "home", lambda: tmp_path / "nohome")
    monkeypatch.setattr(cli.Path, "cwd", lambda: cwd)
    monkeypatch.setenv("A", "existing")

    cli._load_dotenv()
    assert os.environ["A"] == "existing"


def test_show_detailed_timings(capsys, tmp_path):
    log = tmp_path / "workflow.log"
    log.write_text(
        "[12:00:00] Command completed in 1.234s: bcftools view\n"
        "[12:00:01] Command completed in 0.500s: bcftools annotate\n"
    )
    cli._show_detailed_timings(log)
    out = capsys.readouterr().out
    assert "Detailed timing" in out
    assert "bcftools view" in out
    assert "Total" in out


def test_find_cache_dir_variants(tmp_path):
    root = tmp_path / "root"
    cache = root / "cache"
    cache.mkdir(parents=True)

    assert cli._find_cache_dir(root) == cache
    assert cli._find_cache_dir(cache) == cache

    anno = cache / "anno1"
    anno.mkdir()
    (anno / "vcfcache_annotated.bcf").write_text("x")
    assert cli._find_cache_dir(anno) == cache

    with pytest.raises(FileNotFoundError):
        cli._find_cache_dir(tmp_path / "missing")


def test_list_annotation_caches_marks_incomplete(tmp_path):
    cache = tmp_path / "cache"
    cache.mkdir()
    complete = cache / "complete"
    incomplete = cache / "incomplete"
    complete.mkdir()
    incomplete.mkdir()
    (complete / "vcfcache_annotated.bcf.csi").write_text("idx")

    names = cli._list_annotation_caches(cache)
    assert "complete" in names
    assert "incomplete (incomplete)" in names


def test_print_annotation_command_single_cache(tmp_path, capsys):
    cache_root = tmp_path / "root"
    cache_dir = cache_root / "cache" / "anno1"
    cache_dir.mkdir(parents=True)
    (cache_dir / "annotation.yaml").write_text(
        "annotation_cmd: echo hello\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )

    cli._print_annotation_command(cache_root)
    out = capsys.readouterr().out
    assert "Annotation command (with params substituted)" in out
    assert "echo hello" in out


def test_print_annotation_command_multiple_caches(tmp_path, capsys):
    cache_root = tmp_path / "root"
    cache_dir = cache_root / "cache"
    (cache_dir / "a").mkdir(parents=True)
    (cache_dir / "b").mkdir(parents=True)
    (cache_dir / "a" / "vcfcache_annotated.bcf").write_text("x")

    cli._print_annotation_command(cache_root)
    out = capsys.readouterr().out
    assert "Multiple caches found" in out


def test_load_dotenv_quotes_and_comments(monkeypatch, tmp_path):
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    (cwd / ".env").write_text(
        "# comment\n"
        "A='one two'\n"
        "B=\"three\"\n"
        "C=plain\n"
        "INVALID_LINE\n"
    )
    monkeypatch.setattr(cli.Path, "home", lambda: tmp_path / "nohome")
    monkeypatch.setattr(cli.Path, "cwd", lambda: cwd)
    monkeypatch.delenv("A", raising=False)
    monkeypatch.delenv("B", raising=False)
    monkeypatch.delenv("C", raising=False)

    cli._load_dotenv()
    assert os.environ["A"] == "one two"
    assert os.environ["B"] == "three"
    assert os.environ["C"] == "plain"


def test_main_version(monkeypatch, capsys):
    monkeypatch.setattr(cli, "pkg_version", lambda *_args, **_kwargs: "1.2.3")
    monkeypatch.setattr(cli.sys, "argv", ["vcfcache", "--version"])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 0
    out = capsys.readouterr().out.strip()
    assert out == "1.2.3"


def test_main_blueprint_init_local(monkeypatch, tmp_path):
    called = {}

    class _Init:
        def __init__(self, **kwargs):
            called["init_kwargs"] = kwargs

        def initialize(self):
            called["initialized"] = True

    monkeypatch.setattr(cli, "DatabaseInitializer", _Init)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)

    out_dir = tmp_path / "out"
    in_bcf = tmp_path / "in.bcf"
    in_bcf.write_text("bcf")
    monkeypatch.setattr(cli.sys, "argv", ["vcfcache", "blueprint-init", "-i", str(in_bcf), "-o", str(out_dir)])
    cli.main()

    assert called.get("initialized") is True
    assert called["init_kwargs"]["input_file"] == in_bcf
    assert called["init_kwargs"]["output_dir"] == out_dir


def test_main_annotate_requires_args(monkeypatch):
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")
    monkeypatch.setattr(cli.sys, "argv", ["vcfcache", "annotate"])
    with pytest.raises(SystemExit):
        cli.main()


def test_main_cache_build_local(monkeypatch, tmp_path):
    called = {}

    class _Annot:
        def __init__(self, **kwargs):
            called["annot_kwargs"] = kwargs

        def annotate(self, *args, **kwargs):
            called["annotated"] = True

    monkeypatch.setattr(cli, "DatabaseAnnotator", _Annot)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)

    db_dir = tmp_path / "cache_root"
    (db_dir / "blueprint").mkdir(parents=True)
    (db_dir / "blueprint" / "vcfcache.bcf").write_text("bcf")
    anno = tmp_path / "anno.yaml"
    anno.write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "cache-build",
            "-d",
            str(db_dir),
            "-n",
            "test_cache",
            "-a",
            str(anno),
        ],
    )
    cli.main()

    assert called.get("annotated") is True
    assert called["annot_kwargs"]["annotation_name"] == "test_cache"


def test_main_annotate_runs(monkeypatch, tmp_path):
    called = {}

    class _Annot:
        def __init__(self, **kwargs):
            called["annot_kwargs"] = kwargs

        def annotate(self, **kwargs):
            called["ran"] = kwargs

    monkeypatch.setattr(cli, "VCFAnnotator", _Annot)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)

    ann_db = tmp_path / "cache" / "anno"
    ann_db.mkdir(parents=True)
    (ann_db / "vcfcache_annotated.bcf").write_text("x")
    (ann_db / "annotation.yaml").write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )
    (ann_db / "params.snapshot.yaml").write_text(
        "annotation_tool_cmd: bcftools\n"
        "bcftools_cmd: bcftools\n"
        "temp_dir: /tmp\n"
        "threads: 1\n"
        "genome_build: GRCh38\n"
    )

    vcf = tmp_path / "in.bcf"
    vcf.write_text("bcf")
    (tmp_path / "in.bcf.csi").write_text("idx")

    out_dir = tmp_path / "out"
    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "annotate",
            "-a",
            str(ann_db),
            "-i",
            str(vcf),
            "-o",
            str(out_dir),
            "--preserve-unannotated",
        ],
    )
    cli.main()
    assert called.get("ran") is not None
    assert called["ran"]["preserve_unannotated"] is True


def test_main_cache_build_doi_blueprint(monkeypatch, tmp_path):
    called = {}

    class _Annot:
        def __init__(self, **kwargs):
            called["annot_kwargs"] = kwargs

        def annotate(self):
            called["annotated"] = True

    monkeypatch.setattr(cli, "DatabaseAnnotator", _Annot)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(cli, "download_doi", lambda *_args, **_kwargs: None)

    cache_base = tmp_path / "cache_base"
    extract_root = cache_base / "temp"
    extract_root.mkdir(parents=True)
    extracted = extract_root / "bp1"
    (extracted / "blueprint").mkdir(parents=True)
    (extracted / "blueprint" / "vcfcache.bcf").write_text("bcf")

    monkeypatch.setattr(cli, "extract_cache", lambda *_args, **_kwargs: extracted)

    anno = tmp_path / "anno.yaml"
    anno.write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "cache-build",
            "--doi",
            "10.5281/zenodo.123",
            "-o",
            str(cache_base),
            "-n",
            "test_cache",
            "-a",
            str(anno),
        ],
    )
    cli.main()
    assert called.get("annotated") is True
    assert called["annot_kwargs"]["annotation_name"] == "test_cache"


def test_main_annotate_list(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(cli, "_list_annotation_caches", lambda *_args, **_kwargs: ["a", "b"])
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")

    monkeypatch.setattr(cli.sys, "argv", ["vcfcache", "annotate", "--list", "-a", str(tmp_path)])
    cli.main()
    out = capsys.readouterr().out
    assert "Available cached annotations" in out


def test_main_annotate_requirements(monkeypatch, tmp_path):
    called = {}
    monkeypatch.setattr(cli, "_print_annotation_command", lambda *_args, **_kwargs: called.__setitem__("ok", True))
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")

    monkeypatch.setattr(cli.sys, "argv", ["vcfcache", "annotate", "--requirements", "-a", str(tmp_path)])
    cli.main()
    assert called.get("ok") is True


def test_requirements_shows_snapshot_and_values(tmp_path):
    annotation_dir = tmp_path / "cache" / "anno1"
    annotation_dir.mkdir(parents=True)
    (annotation_dir / "annotation.yaml").write_text(
        "annotation_cmd: \"echo ${params.foo} ${params.bar}\"\n"
        "must_contain_info_tag: CSQ\n"
        "required_tool_version: \"1.0\"\n"
        "genome_build: GRCh38\n"
    )
    (annotation_dir / "params.snapshot.yaml").write_text(
        "foo: hello\n"
        "bcftools_cmd: bcftools\n"
    )

    result = subprocess.run(
        [sys.executable, "-m", "vcfcache.cli", "annotate", "--requirements", "-a", str(annotation_dir)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "params.yaml: " in result.stdout
    assert "from cache snapshot" in result.stdout
    assert "foo: hello" in result.stdout
    assert "bar: " in result.stdout
    assert "<missing>" in result.stdout


def test_main_cache_build_doi_prebuilt_cache(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(cli, "download_doi", lambda *_args, **_kwargs: None)

    def _should_not_call(*_args, **_kwargs):
        raise AssertionError("DatabaseAnnotator should not be called for prebuilt cache")

    monkeypatch.setattr(cli, "DatabaseAnnotator", _should_not_call)

    cache_base = tmp_path / "cache_base"
    extract_root = cache_base / "temp"
    extract_root.mkdir(parents=True)
    extracted = extract_root / "cache1"
    (extracted / "blueprint").mkdir(parents=True)
    (extracted / "blueprint" / "vcfcache.bcf").write_text("bcf")
    anno_dir = extracted / "cache" / "anno1"
    anno_dir.mkdir(parents=True)
    (anno_dir / "vcfcache_annotated.bcf").write_text("x")
    (anno_dir / "annotation.yaml").write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )

    monkeypatch.setattr(cli, "extract_cache", lambda *_args, **_kwargs: extracted)

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "cache-build",
            "--doi",
            "10.5281/zenodo.999",
            "-o",
            str(cache_base),
        ],
    )
    cli.main()


def test_main_list_inspect_blueprint_only(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(cli, "search_zenodo_records", lambda *_args, **_kwargs: [])

    root = tmp_path / "bp_only"
    (root / "blueprint").mkdir(parents=True)
    (root / "blueprint" / "vcfcache.bcf").write_text("bcf")

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "list",
            "blueprints",
            "--inspect",
            str(root),
        ],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "VCFcache inspect" in out
    assert "blueprint-only" in out.lower()


def test_main_list_local_blueprints(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.setattr("vcfcache.utils.validation.find_bcftools", lambda: "/usr/bin/bcftools")

    def _run(*_args, **_kwargs):
        return type("R", (), {"stdout": "5\n", "returncode": 0})()

    monkeypatch.setattr("subprocess.run", _run)

    base = tmp_path / "vcfcache"
    blueprints = base / "blueprints" / "bp1"
    (blueprints / "blueprint").mkdir(parents=True)
    (blueprints / "blueprint" / "vcfcache.bcf").write_bytes(b"x" * (1024 * 1024 + 10))
    (blueprints / "blueprint" / "sources.info").write_text("{}")
    (blueprints / ".vcfcache_complete").write_text("completed: true\nmode: blueprint-init\n")

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "list",
            "blueprints",
            "--local",
            str(base),
        ],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "Local vcfcache blueprints" in out
    assert "Variants:" in out


def test_main_list_local_caches(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.setattr("vcfcache.utils.validation.find_bcftools", lambda: "/usr/bin/bcftools")

    def _run(*_args, **_kwargs):
        return type("R", (), {"stdout": "7\n", "returncode": 0})()

    monkeypatch.setattr("subprocess.run", _run)

    base = tmp_path / "vcfcache"
    cache_root = base / "caches" / "cache1"
    (cache_root / "blueprint").mkdir(parents=True)
    (cache_root / "blueprint" / "vcfcache.bcf").write_bytes(b"x" * (1024 * 1024 + 10))
    (cache_root / "blueprint" / "sources.info").write_text("{}")
    anno = cache_root / "cache" / "anno1"
    anno.mkdir(parents=True)
    (anno / "vcfcache_annotated.bcf").write_bytes(b"x" * (1024 * 1024 + 10))
    (anno / ".vcfcache_complete").write_text("completed: true\nmode: cache-build\n")
    (anno / "annotation.yaml").write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )
    (anno / "params.snapshot.yaml").write_text(
        "annotation_tool_cmd: bcftools\n"
        "bcftools_cmd: bcftools\n"
        "temp_dir: /tmp\n"
        "threads: 1\n"
        "genome_build: GRCh38\n"
    )

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "list",
            "caches",
            "--local",
            str(base),
        ],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "Local vcfcache caches" in out


def test_main_list_remote_no_records(monkeypatch, capsys):
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(cli, "search_zenodo_records", lambda *_args, **_kwargs: [])

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "list",
            "caches",
        ],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "No caches found on Zenodo" in out


def test_main_list_inspect_annotation_dir(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)

    root = tmp_path / "cache_root"
    anno = root / "cache" / "anno1"
    anno.mkdir(parents=True)
    (root / "blueprint").mkdir(parents=True)
    (root / "blueprint" / "vcfcache.bcf").write_text("bcf")
    (anno / "vcfcache_annotated.bcf").write_text("x")
    (anno / "annotation.yaml").write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )
    (anno / "params.snapshot.yaml").write_text(
        "annotation_tool_cmd: bcftools\n"
        "bcftools_cmd: bcftools\n"
        "temp_dir: /tmp\n"
        "threads: 1\n"
        "genome_build: GRCh38\n"
    )

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "list",
            "caches",
            "--inspect",
            str(anno),
        ],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "Annotation dir" in out


def test_main_list_remote_records(monkeypatch, capsys):
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)

    records = [
        {
            "title": "cache-foo",
            "keywords": ["cache-hg38-gnomad-v4-AF0100-vep-110-basic"],
            "doi": "10.1234/zenodo.1",
            "created": "2024-01-01",
            "size_mb": 10.5,
        }
    ]
    monkeypatch.setattr(cli, "search_zenodo_records", lambda *_args, **_kwargs: records)

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "list",
            "caches",
        ],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "Available vcfcache caches on Zenodo" in out
    assert "DOI:" in out


def test_main_cache_build_local_missing_args(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")

    db_dir = tmp_path / "cache_root"
    (db_dir / "blueprint").mkdir(parents=True)
    (db_dir / "blueprint" / "vcfcache.bcf").write_text("bcf")

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "cache-build",
            "-d",
            str(db_dir),
        ],
    )
    with pytest.raises(ValueError):
        cli.main()


def test_main_annotate_list_empty(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(cli, "_list_annotation_caches", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")

    monkeypatch.setattr(cli.sys, "argv", ["vcfcache", "annotate", "--list", "-a", str(tmp_path)])
    cli.main()
    out = capsys.readouterr().out
    assert "No cached annotations found" in out


def test_main_push_invalid_cache_dir(monkeypatch, tmp_path):
    monkeypatch.setenv("ZENODO_TOKEN", "token")
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)

    bad_dir = tmp_path / "bad"
    bad_dir.mkdir()

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "push",
            "--cache-dir",
            str(bad_dir),
        ],
    )
    with pytest.raises(ValueError):
        cli.main()


def test_main_demo_smoke(monkeypatch):
    """Test that demo command runs smoke test successfully."""
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)

    demo_mod = type(
        "D",
        (),
        {"run_smoke_test": lambda *a, **_k: 0},  # Only smoke test now
    )()
    monkeypatch.setitem(sys.modules, "vcfcache.demo", demo_mod)

    # Demo now always runs smoke test (no --smoke-test flag needed)
    monkeypatch.setattr(cli.sys, "argv", ["vcfcache", "demo"])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 0


def test_main_cache_build_doi_missing_name(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")
    monkeypatch.setattr(cli, "download_doi", lambda *_a, **_k: None)

    cache_base = tmp_path / "base"
    extract_root = cache_base / "temp"
    extract_root.mkdir(parents=True)
    extracted = extract_root / "bp1"
    (extracted / "blueprint").mkdir(parents=True)
    (extracted / "blueprint" / "vcfcache.bcf").write_text("bcf")
    monkeypatch.setattr(cli, "extract_cache", lambda *_a, **_k: extracted)

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "cache-build",
            "--doi",
            "10.5281/zenodo.1",
            "-o",
            str(cache_base),
            "-a",
            str(tmp_path / "anno.yaml"),
        ],
    )
    with pytest.raises(ValueError):
        cli.main()


def test_main_cache_build_doi_cache_with_anno(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")
    monkeypatch.setattr(cli, "download_doi", lambda *_a, **_k: None)

    cache_base = tmp_path / "base"
    extract_root = cache_base / "temp"
    extract_root.mkdir(parents=True)
    extracted = extract_root / "cache1"
    (extracted / "blueprint").mkdir(parents=True)
    (extracted / "blueprint" / "vcfcache.bcf").write_text("bcf")
    anno_dir = extracted / "cache" / "anno1"
    anno_dir.mkdir(parents=True)
    (anno_dir / "vcfcache_annotated.bcf").write_text("x")
    (anno_dir / "annotation.yaml").write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )
    monkeypatch.setattr(cli, "extract_cache", lambda *_a, **_k: extracted)

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "cache-build",
            "--doi",
            "10.5281/zenodo.2",
            "-o",
            str(cache_base),
            "-a",
            str(tmp_path / "anno.yaml"),
        ],
    )
    with pytest.raises(ValueError):
        cli.main()


def test_main_push_metadata_error(monkeypatch, tmp_path):
    monkeypatch.setenv("ZENODO_TOKEN", "token")
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)

    cache_dir = tmp_path / "cache"
    (cache_dir / "blueprint").mkdir(parents=True)
    (cache_dir / "blueprint" / "vcfcache.bcf").write_text("bcf")
    (cache_dir / ".vcfcache_complete").write_text("completed: true\nmode: blueprint-init\n")

    import types
    import vcfcache.utils.archive as archive_mod
    import vcfcache.integrations as integrations_mod

    zenodo_mod = types.ModuleType("zenodo")
    zenodo_mod.create_deposit = lambda *_args, **_kwargs: {"id": 1, "doi": "10.1234/zenodo.1"}
    zenodo_mod.upload_file = lambda *_args, **_kwargs: None
    zenodo_mod.publish_deposit = lambda dep, *_args, **_kwargs: dep
    zenodo_mod._api_base = lambda *_args, **_kwargs: "https://example.org"

    monkeypatch.setattr(integrations_mod, "zenodo", zenodo_mod)
    monkeypatch.setitem(sys.modules, "vcfcache.integrations.zenodo", zenodo_mod)

    def _tar_cache(*_args, **_kwargs):
        tar_path = _args[1] if len(_args) > 1 else _kwargs.get("tar_path")
        Path(tar_path).write_text("tar")
        return Path(tar_path)

    monkeypatch.setattr(cli, "tar_cache_subset", _tar_cache)
    monkeypatch.setattr(archive_mod, "file_md5", lambda *_args, **_kwargs: "md5")

    class _Resp:
        ok = False
        status_code = 400
        reason = "Bad"

        def json(self):
            return {"error": "bad"}

        def raise_for_status(self):
            pass

    monkeypatch.setattr(cli.requests, "put", lambda *_a, **_k: _Resp())

    meta = tmp_path / "meta.json"
    meta.write_text('{"title": "x"}')

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "push",
            "--cache-dir",
            str(cache_dir),
            "--metadata",
            str(meta),
            "--yes",
        ],
    )
    with pytest.raises(RuntimeError):
        cli.main()


def test_main_cache_build_doi_blueprint_success(monkeypatch, tmp_path):
    called = {}
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")
    monkeypatch.setattr(cli, "download_doi", lambda *_a, **_k: None)

    class _Annot:
        def __init__(self, **kwargs):
            called["ok"] = True

        def annotate(self):
            called["annotated"] = True

    monkeypatch.setattr(cli, "DatabaseAnnotator", _Annot)

    cache_base = tmp_path / "base"
    extract_root = cache_base / "temp"
    extract_root.mkdir(parents=True)
    extracted = extract_root / "bp1"
    (extracted / "blueprint").mkdir(parents=True)
    (extracted / "blueprint" / "vcfcache.bcf").write_text("bcf")

    monkeypatch.setattr(cli, "extract_cache", lambda *_a, **_k: extracted)

    anno = tmp_path / "anno.yaml"
    anno.write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "cache-build",
            "--doi",
            "10.5281/zenodo.1",
            "-o",
            str(cache_base),
            "-n",
            "test_cache",
            "-a",
            str(anno),
        ],
    )
    cli.main()
    assert called.get("annotated") is True


def test_main_list_inspect_cache_root_multiple_annos(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)

    root = tmp_path / "cache_root"
    (root / "blueprint").mkdir(parents=True)
    (root / "blueprint" / "vcfcache.bcf").write_text("bcf")
    anno1 = root / "cache" / "a"
    anno2 = root / "cache" / "b"
    anno1.mkdir(parents=True)
    anno2.mkdir(parents=True)
    (anno1 / "vcfcache_annotated.bcf").write_text("x")
    (anno1 / "annotation.yaml").write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )
    (anno2 / "vcfcache_annotated.bcf").write_text("x")
    (anno2 / "annotation.yaml").write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "list",
            "caches",
            "--inspect",
            str(root),
        ],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "Cache annotations" in out


def test_main_push_publish_autometadata(monkeypatch, tmp_path):
    monkeypatch.setenv("ZENODO_TOKEN", "token")
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)

    cache_dir = tmp_path / "cache"
    (cache_dir / "blueprint").mkdir(parents=True)
    (cache_dir / "blueprint" / "vcfcache.bcf").write_text("bcf")
    (cache_dir / ".vcfcache_complete").write_text("completed: true\nmode: blueprint-init\n")

    import types
    import vcfcache.utils.archive as archive_mod
    import vcfcache.integrations as integrations_mod

    zenodo_mod = types.ModuleType("zenodo")
    zenodo_mod.create_deposit = lambda *_args, **_kwargs: {"id": 1, "doi": "10.1234/zenodo.1"}
    zenodo_mod.upload_file = lambda *_args, **_kwargs: None
    zenodo_mod.publish_deposit = lambda dep, *_args, **_kwargs: dep
    zenodo_mod._api_base = lambda *_args, **_kwargs: "https://example.org"

    monkeypatch.setattr(integrations_mod, "zenodo", zenodo_mod)
    monkeypatch.setitem(sys.modules, "vcfcache.integrations.zenodo", zenodo_mod)

    def _tar_cache(*_args, **_kwargs):
        tar_path = _args[1] if len(_args) > 1 else _kwargs.get("tar_path")
        Path(tar_path).write_text("tar")
        return Path(tar_path)

    monkeypatch.setattr(cli, "tar_cache_subset", _tar_cache)
    monkeypatch.setattr(archive_mod, "file_md5", lambda *_args, **_kwargs: "md5")
    monkeypatch.setattr(cli.requests, "put", lambda *_a, **_k: type("Resp", (), {"ok": True, "raise_for_status": lambda *_: None})())

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "push",
            "--cache-dir",
            str(cache_dir),
            "--publish",
            "--yes",
        ],
    )
    cli.main()


def test_main_push_test_token(monkeypatch, tmp_path):
    monkeypatch.setenv("ZENODO_SANDBOX_TOKEN", "token")
    monkeypatch.delenv("ZENODO_TOKEN", raising=False)
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)

    cache_dir = tmp_path / "cache"
    (cache_dir / "blueprint").mkdir(parents=True)
    (cache_dir / "blueprint" / "vcfcache.bcf").write_text("bcf")
    (cache_dir / ".vcfcache_complete").write_text("completed: true\nmode: blueprint-init\n")

    import types
    import vcfcache.utils.archive as archive_mod
    import vcfcache.integrations as integrations_mod

    zenodo_mod = types.ModuleType("zenodo")
    zenodo_mod.create_deposit = lambda *_args, **_kwargs: {"id": 1, "doi": "10.1234/zenodo.1"}
    zenodo_mod.upload_file = lambda *_args, **_kwargs: None
    zenodo_mod.publish_deposit = lambda dep, *_args, **_kwargs: dep
    zenodo_mod._api_base = lambda *_args, **_kwargs: "https://example.org"

    monkeypatch.setattr(integrations_mod, "zenodo", zenodo_mod)
    monkeypatch.setitem(sys.modules, "vcfcache.integrations.zenodo", zenodo_mod)

    def _tar_cache(*_args, **_kwargs):
        tar_path = _args[1] if len(_args) > 1 else _kwargs.get("tar_path")
        Path(tar_path).write_text("tar")
        return Path(tar_path)

    monkeypatch.setattr(cli, "tar_cache_subset", _tar_cache)
    monkeypatch.setattr(archive_mod, "file_md5", lambda *_args, **_kwargs: "md5")
    monkeypatch.setattr(cli.requests, "put", lambda *_a, **_k: type("Resp", (), {"ok": True, "raise_for_status": lambda *_: None})())

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "push",
            "--cache-dir",
            str(cache_dir),
            "--test",
            "--yes",
        ],
    )
    cli.main()


def test_main_list_inspect_invalid_path(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)
    monkeypatch.setattr(cli, "search_zenodo_records", lambda *_a, **_k: [])

    bad = tmp_path / "nope"
    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "list",
            "caches",
            "--inspect",
            str(bad),
        ],
    )
    with pytest.raises(FileNotFoundError):
        cli.main()


def test_main_annotate_missing_args_errors(monkeypatch):
    monkeypatch.setattr(cli.sys, "argv", ["vcfcache", "annotate", "-a", "cache"])
    with pytest.raises(SystemExit):
        cli.main()


def test_main_blueprint_init_doi_downloads(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")

    tar_holder = {}

    def _download_doi(_doi, tar_path, sandbox=False):
        tar_holder["path"] = Path(tar_path)
        Path(tar_path).write_text("tar")

    def _extract_cache(_tar_path, output_dir):
        extracted = Path(output_dir) / "bp1"
        extracted.mkdir(parents=True, exist_ok=True)
        return extracted

    monkeypatch.setattr(cli, "download_doi", _download_doi)
    monkeypatch.setattr(cli, "extract_cache", _extract_cache)

    out_dir = tmp_path / "out"
    monkeypatch.setattr(
        cli.sys,
        "argv",
        ["vcfcache", "blueprint-init", "--doi", "10.5281/zenodo.1", "-o", str(out_dir)],
    )
    cli.main()

    assert tar_holder["path"].exists() is False
    assert (out_dir / "bp1").exists()


def test_main_blueprint_init_debug_shows_timings(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")

    called = {"timings": None}

    def _show_timing(path):
        called["timings"] = Path(path)

    class _Init:
        def __init__(self, **_kwargs):
            pass

        def initialize(self):
            pass

    monkeypatch.setattr(cli, "_show_detailed_timings", _show_timing)
    monkeypatch.setattr(cli, "DatabaseInitializer", _Init)

    out_dir = tmp_path / "out"
    (out_dir / "blueprint").mkdir(parents=True, exist_ok=True)
    (out_dir / "blueprint" / "workflow.log").write_text("Command completed in 1.0s: step\n")
    vcf = tmp_path / "input.vcf"
    vcf.write_text("vcf")
    (tmp_path / "input.vcf.tbi").write_text("idx")

    monkeypatch.setattr(
        cli.sys,
        "argv",
        ["vcfcache", "blueprint-init", "-i", str(vcf), "-o", str(out_dir), "--debug"],
    )
    cli.main()
    assert called["timings"] == out_dir / "blueprint" / "workflow.log"


def test_main_blueprint_extend_debug_shows_timings(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")

    called = {"timings": None}

    def _show_timing(path):
        called["timings"] = Path(path)

    class _Upd:
        def __init__(self, **_kwargs):
            pass

        def add(self):
            pass

    monkeypatch.setattr(cli, "_show_detailed_timings", _show_timing)
    monkeypatch.setattr(cli, "DatabaseUpdater", _Upd)

    db_dir = tmp_path / "db"
    (db_dir / "blueprint").mkdir(parents=True, exist_ok=True)
    (db_dir / "blueprint" / "workflow.log").write_text("Command completed in 1.0s: step\n")

    vcf = tmp_path / "input.vcf"
    vcf.write_text("vcf")
    (tmp_path / "input.vcf.tbi").write_text("idx")

    monkeypatch.setattr(
        cli.sys,
        "argv",
        ["vcfcache", "blueprint-extend", "-d", str(db_dir), "-i", str(vcf), "--debug"],
    )
    cli.main()
    assert called["timings"] == db_dir / "blueprint" / "workflow.log"


def test_show_detailed_timings_output(tmp_path, capsys):
    log = tmp_path / "workflow.log"
    log.write_text(
        "Command completed in 75.5s: step1\nCommand completed in 10.0s: step2\n"
    )
    cli._show_detailed_timings(log)
    out = capsys.readouterr().out
    assert "Detailed timing" in out
    assert "Total" in out


def test_show_detailed_timings_missing_file(tmp_path, capsys):
    cli._show_detailed_timings(tmp_path / "missing.log")
    assert capsys.readouterr().out == ""


def test_print_annotation_command_single_cache(tmp_path, capsys):
    cache_root = tmp_path / "cache_root"
    cache_dir = cache_root / "cache"
    cache_dir.mkdir(parents=True)
    anno_dir = cache_dir / "anno1"
    anno_dir.mkdir()
    (anno_dir / "annotation.yaml").write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )

    cli._print_annotation_command(cache_root)
    out = capsys.readouterr().out
    assert "Annotation command (with params substituted)" in out
    assert "echo ok" in out


def test_print_annotation_command_multiple_caches(tmp_path, capsys):
    cache_root = tmp_path / "cache_root"
    cache_dir = cache_root / "cache"
    cache_dir.mkdir(parents=True)
    (cache_dir / "a").mkdir()
    (cache_dir / "b").mkdir()
    (cache_dir / "a" / "vcfcache_annotated.bcf").write_text("bcf")

    cli._print_annotation_command(cache_root)
    out = capsys.readouterr().out
    assert "Multiple caches found" in out
    assert " (incomplete)" in out


def test_print_annotation_command_no_caches(tmp_path):
    cache_root = tmp_path / "cache_root"
    (cache_root / "cache").mkdir(parents=True)
    with pytest.raises(FileNotFoundError):
        cli._print_annotation_command(cache_root)


def test_list_annotation_caches_marks_incomplete(tmp_path):
    cache_root = tmp_path / "cache_root"
    cache_dir = cache_root / "cache"
    cache_dir.mkdir(parents=True)
    (cache_dir / "complete").mkdir()
    (cache_dir / "complete" / "vcfcache_annotated.bcf.csi").write_text("idx")
    (cache_dir / "incomplete").mkdir()
    (cache_dir / "note.txt").write_text("skip")

    names = cli._list_annotation_caches(cache_root)
    assert names == ["complete", "incomplete (incomplete)"]


def test_main_list_inspect_blueprint_only(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)

    root = tmp_path / "bp"
    (root / "blueprint").mkdir(parents=True)
    (root / "blueprint" / "vcfcache.bcf").write_text("bcf")

    monkeypatch.setattr(
        cli.sys,
        "argv",
        ["vcfcache", "list", "caches", "--inspect", str(root)],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "blueprint-only" in out


def test_main_list_inspect_invalid_structure(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)

    bad = tmp_path / "bad"
    bad.mkdir()
    monkeypatch.setattr(
        cli.sys,
        "argv",
        ["vcfcache", "list", "caches", "--inspect", str(bad)],
    )
    with pytest.raises(ValueError):
        cli.main()


def test_main_list_inspect_missing_params_snapshot(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)

    root = tmp_path / "cache_root"
    (root / "blueprint").mkdir(parents=True)
    (root / "blueprint" / "vcfcache.bcf").write_text("bcf")
    anno_dir = root / "cache" / "anno1"
    anno_dir.mkdir(parents=True)
    (anno_dir / "vcfcache_annotated.bcf").write_text("bcf")
    (anno_dir / "annotation.yaml").write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
        "# ${params.foo}\n"
    )

    monkeypatch.setattr(
        cli.sys,
        "argv",
        ["vcfcache", "list", "caches", "--inspect", str(anno_dir)],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "Params snapshot: (missing) params.snapshot.yaml" in out
    assert "MISSING" in out


def test_main_list_local_caches_formats_alias(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)

    base = tmp_path / "local"
    caches_dir = base / "caches"
    caches_dir.mkdir(parents=True)
    alias = "cache-GRCh38-gnomad-3.1-AF0100-vep-110-default"
    root = caches_dir / alias
    (root / "blueprint").mkdir(parents=True)
    (root / "cache").mkdir()
    anno_dir = root / "cache" / "anno1"
    anno_dir.mkdir(parents=True)
    (anno_dir / "annotation.yaml").write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )
    (anno_dir / "vcfcache_annotated.bcf").write_bytes(b"x" * (2 * 1024 * 1024))
    (anno_dir / "vcfcache_annotated.bcf.csi").write_text("idx")
    (anno_dir / ".vcfcache_complete").write_text("completed: true\nmode: cache-build\n")
    (root / "blueprint" / "vcfcache.bcf").write_text("bcf")
    (root / "blueprint" / "vcfcache.bcf.csi").write_text("idx")
    (root / "blueprint" / "sources.info").write_text("{}")
    (anno_dir / "params.snapshot.yaml").write_text(
        "annotation_tool_cmd: bcftools\n"
        "bcftools_cmd: bcftools\n"
        "temp_dir: /tmp\n"
        "threads: 1\n"
        "genome_build: GRCh38\n"
    )

    monkeypatch.setattr(
        cli.sys,
        "argv",
        ["vcfcache", "list", "caches", "--local", str(base)],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "gnomAD" in out
    assert "VEP" in out
    assert "AF ≥ 0.10" in out


def test_main_list_local_missing_dir(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)

    missing = tmp_path / "missing"
    monkeypatch.setattr(
        cli.sys,
        "argv",
        ["vcfcache", "list", "caches", "--local", str(missing)],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "No local caches found" in out


def test_push_requires_completion_flags(monkeypatch, tmp_path):
    monkeypatch.setattr(
        cli,
        "setup_logging",
        lambda *_args, **_kwargs: type(
            "L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None}
        )(),
    )
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(cli.os, "environ", {**cli.os.environ, "ZENODO_TOKEN": "token"})
    monkeypatch.setattr(cli, "tar_cache_subset", lambda *_args, **_kwargs: None)
    monkeypatch.setattr("vcfcache.utils.archive.file_md5", lambda *_args, **_kwargs: "md5")
    monkeypatch.setattr("vcfcache.integrations.zenodo.create_deposit", lambda *_args, **_kwargs: {"id": 1})
    monkeypatch.setattr("vcfcache.integrations.zenodo.upload_file", lambda *_args, **_kwargs: None)

    # Missing blueprint completion flag
    base = tmp_path / "bp_root"
    (base / "blueprint").mkdir(parents=True)
    (base / "blueprint" / "vcfcache.bcf").write_text("bcf")

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "push",
            "--cache-dir",
            str(base),
            "--yes",
        ],
    )
    with pytest.raises(ValueError, match="Missing .vcfcache_complete"):
        cli.main()

    # Cache with missing per-annotation completion flag (select cache explicitly)
    cache_root = tmp_path / "cache_root"
    (cache_root / "blueprint").mkdir(parents=True)
    (cache_root / "blueprint" / "vcfcache.bcf").write_text("bcf")
    (cache_root / "cache" / "anno1").mkdir(parents=True)
    (cache_root / "cache" / "anno1" / "vcfcache_annotated.bcf").write_text("bcf")
    (cache_root / ".vcfcache_complete").write_text("completed: true\nmode: blueprint-init\n")

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "push",
            "--cache-dir",
            str(cache_root / "cache" / "anno1"),
            "--yes",
        ],
    )
    with pytest.raises(ValueError, match="Missing .vcfcache_complete"):
        cli.main()


def test_main_version_fallback_uses_package_version(monkeypatch):
    monkeypatch.setattr(cli, "pkg_version", lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("nope")))
    monkeypatch.setattr(cli.sys, "argv", ["vcfcache", "--version"])
    with pytest.raises(SystemExit):
        cli.main()


def test_main_annotate_requirements_and_list_conflict(monkeypatch):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "annotate",
            "--requirements",
            "--list",
            "-a",
            "cache",
        ],
    )
    with pytest.raises(SystemExit):
        cli.main()


def test_main_annotate_uncached_print(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(cli, "setup_logging", lambda *_a, **_k: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_a, **_k: None)
    monkeypatch.setattr(cli, "check_bcftools_installed", lambda: "/usr/bin/bcftools")

    class _Annot:
        def __init__(self, **kwargs):
            pass

        def annotate(self, **kwargs):
            return None

    monkeypatch.setattr(cli, "VCFAnnotator", _Annot)

    ann_db = tmp_path / "cache" / "anno"
    ann_db.mkdir(parents=True)
    (ann_db / "vcfcache_annotated.bcf").write_text("x")
    (ann_db / "annotation.yaml").write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )
    (ann_db / "params.snapshot.yaml").write_text(
        "annotation_tool_cmd: bcftools\n"
        "bcftools_cmd: bcftools\n"
        "temp_dir: /tmp\n"
        "threads: 1\n"
        "genome_build: GRCh38\n"
    )

    vcf = tmp_path / "in.bcf"
    vcf.write_text("bcf")
    (tmp_path / "in.bcf.csi").write_text("idx")

    out_dir = tmp_path / "out"
    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "annotate",
            "-a",
            str(ann_db),
            "-i",
            str(vcf),
            "-o",
            str(out_dir),
            "--uncached",
        ],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "uncached mode" in out


def test_main_push_requires_token(monkeypatch, tmp_path):
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)
    monkeypatch.delenv("ZENODO_TOKEN", raising=False)
    monkeypatch.delenv("ZENODO_SANDBOX_TOKEN", raising=False)
    monkeypatch.setenv("ZENODO_TOKEN", "")

    cache_dir = tmp_path / "cache"
    (cache_dir / "blueprint").mkdir(parents=True)
    (cache_dir / "blueprint" / "vcfcache.bcf").write_text("bcf")
    (cache_dir / ".vcfcache_complete").write_text("completed: true\nmode: blueprint-init\n")

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "push",
            "--cache-dir",
            str(cache_dir),
        ],
    )
    with pytest.raises(RuntimeError):
        cli.main()


def test_main_push_blueprint_success(monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("ZENODO_TOKEN", "token")
    monkeypatch.setattr(cli, "setup_logging", lambda *_args, **_kwargs: type("L", (), {"debug": lambda *_: None, "info": lambda *_: None, "error": lambda *_: None})())
    monkeypatch.setattr(cli, "log_command", lambda *_args, **_kwargs: None)

    cache_dir = tmp_path / "cache"
    (cache_dir / "blueprint").mkdir(parents=True)
    (cache_dir / "blueprint" / "vcfcache.bcf").write_text("bcf")
    (cache_dir / ".vcfcache_complete").write_text("completed: true\nmode: blueprint-init\n")

    import types
    import vcfcache.utils.archive as archive_mod
    import vcfcache.integrations as integrations_mod

    zenodo_mod = types.ModuleType("zenodo")
    zenodo_mod.create_deposit = lambda *_args, **_kwargs: {"id": 1, "doi": "10.1234/zenodo.1"}
    zenodo_mod.upload_file = lambda *_args, **_kwargs: None
    zenodo_mod.publish_deposit = lambda dep, *_args, **_kwargs: dep

    monkeypatch.setattr(integrations_mod, "zenodo", zenodo_mod)
    monkeypatch.setitem(sys.modules, "vcfcache.integrations.zenodo", zenodo_mod)

    def _tar_cache(*_args, **_kwargs):
        tar_path = _args[1] if len(_args) > 1 else _kwargs.get("tar_path")
        Path(tar_path).write_text("tar")
        return Path(tar_path)

    monkeypatch.setattr(cli, "tar_cache_subset", _tar_cache)
    monkeypatch.setattr(archive_mod, "file_md5", lambda *_args, **_kwargs: "md5")

    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "vcfcache",
            "push",
            "--cache-dir",
            str(cache_dir),
            "--yes",
        ],
    )
    cli.main()
    out = capsys.readouterr().out
    assert "Upload complete" in out
