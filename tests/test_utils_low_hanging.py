# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

import logging
import runpy
import sys
from pathlib import Path

import pytest

from vcfcache.utils.archive import dir_md5
from vcfcache.utils.logging import log_command, setup_logging
from vcfcache.utils.naming import CacheName
from vcfcache.utils.paths import get_project_root, get_resource_path, get_vcfcache_root
from vcfcache.utils.validation import (
    check_duplicate_md5,
    compare_versions,
    find_bcftools,
    generate_test_command,
    parse_bcftools_version,
    validate_vcf_format,
)


def test_cache_name_blueprint_alias_roundtrip():
    name = CacheName(genome="hg38", source="gnomad", release="v4", filt="AF0100")
    assert name.is_cache is False
    assert name.blueprint_alias == "bp-hg38-gnomad-v4-AF0100"
    with pytest.raises(ValueError):
        _ = name.cache_alias

    parsed = CacheName.parse(name.blueprint_alias)
    assert parsed == name


def test_cache_name_cache_alias_roundtrip():
    name = CacheName(
        genome="hg38",
        source="gnomad",
        release="v4",
        filt="AF0100",
        tool="vep",
        tool_version="110",
        preset="my-preset",
    )
    assert name.is_cache is True
    assert (
        name.cache_alias
        == "cache-hg38-gnomad-v4-AF0100-vep-110-my-preset"
    )

    parsed = CacheName.parse(name.cache_alias)
    assert parsed == name


def test_cache_name_parse_invalid():
    with pytest.raises(ValueError):
        CacheName.parse("cache-too-few-parts")
    with pytest.raises(ValueError):
        CacheName.parse("unknown-hg38-gnomad")


def test_setup_logging_levels_and_file(tmp_path):
    log_path = tmp_path / "logs" / "vcfcache.log"
    logger = setup_logging(verbosity=2, log_file=log_path)

    # Ensure handlers are configured as expected.
    assert any(isinstance(h, logging.FileHandler) for h in logger.handlers)
    assert log_path.exists()

    # Console handler should be DEBUG in verbosity=2
    console = next(h for h in logger.handlers if not isinstance(h, logging.FileHandler))
    assert console.level == logging.DEBUG


def test_log_command_records_message(monkeypatch, caplog):
    logger = setup_logging(verbosity=2)
    monkeypatch.setattr(sys, "argv", ["vcfcache", "demo", "--quiet"])
    with caplog.at_level(logging.DEBUG, logger="vcfcache"):
        log_command(logger)
    assert "Script command:" in caplog.text


def test_paths_env_override(monkeypatch, tmp_path):
    monkeypatch.setenv("VCFCACHE_ROOT", str(tmp_path))
    assert get_project_root() == tmp_path
    assert get_vcfcache_root() == tmp_path
    assert get_resource_path(Path("resources/x.txt")) == tmp_path / "resources/x.txt"


def test_dir_md5_changes_with_content(tmp_path):
    a = tmp_path / "a.txt"
    b = tmp_path / "b.txt"
    a.write_bytes(b"one")
    b.write_bytes(b"two")
    first = dir_md5([a, b])
    b.write_bytes(b"two!")
    second = dir_md5([a, b])
    assert first != second


def test_parse_and_compare_versions():
    assert parse_bcftools_version("1.20") == (1, 20, 0)
    assert parse_bcftools_version("1.22+htslib-1.22") == (1, 22, 0)
    assert parse_bcftools_version("1.2.3") == (1, 2, 3)
    with pytest.raises(ValueError):
        parse_bcftools_version("not-a-version")

    assert compare_versions("1.20", "1.21") == -1
    assert compare_versions("1.21", "1.21") == 0
    assert compare_versions("1.22", "1.21") == 1


def test_find_bcftools_env_preferred(monkeypatch, tmp_path):
    fake = tmp_path / "bcftools"
    fake.write_text("#!/bin/sh\necho ok\n")
    fake.chmod(0o755)
    monkeypatch.setenv("VCFCACHE_BCFTOOLS", str(fake))
    assert find_bcftools() == str(fake.resolve())


def test_find_bcftools_env_invalid_falls_back(monkeypatch):
    monkeypatch.setenv("VCFCACHE_BCFTOOLS", "/no/such/bcftools")
    monkeypatch.setattr("vcfcache.utils.validation.shutil.which", lambda _: "/usr/bin/bcftools")
    assert find_bcftools() == "/usr/bin/bcftools"


def test_check_duplicate_md5():
    db_info = {"input_files": [{"md5": "abc"}, {"md5": "def"}]}
    assert check_duplicate_md5(db_info, "abc") is True
    assert check_duplicate_md5(db_info, "xyz") is False
    assert check_duplicate_md5({}, "abc") is False


def test_generate_test_command_output(capsys):
    cmd = generate_test_command(
        vcfcache_path="vcfcache",
        vcf_path="input.vcf",
        output_dir="/tmp/out",
        yaml_path="params.yaml",
        annotation_config="ann.yaml",
        add_vcf_path="add.vcf",
        input_vcf_path="query.vcf",
        annotate_name="demo",
        annotation_db="/tmp/db",
        annotation_output="/tmp/out.vcf",
        force=False,
    )
    captured = capsys.readouterr().out
    assert "blueprint-init" in captured
    assert "cache-build" in cmd


def test_validate_vcf_format_ok_and_missing(tmp_path):
    good_vcf = tmp_path / "good.vcf"
    good_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "1\t1\t.\tA\tT\t.\tPASS\t.\tGT:AD\t0/1:3,4\n"
    )
    ok, err = validate_vcf_format(good_vcf)
    assert ok is True
    assert err is None

    bad_vcf = tmp_path / "bad.vcf"
    bad_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "1\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
    )
    ok, err = validate_vcf_format(bad_vcf)
    assert ok is False
    assert err is not None


def test_main_entrypoint_invokes_cli(monkeypatch):
    import vcfcache.cli as cli

    called = {"ok": False}

    def _main():
        called["ok"] = True

    monkeypatch.setattr(cli, "main", _main)
    runpy.run_module("vcfcache.__main__", run_name="__main__")
    assert called["ok"] is True


def test_main_module_import_no_exec():
    import importlib

    mod = importlib.import_module("vcfcache.__main__")
    assert hasattr(mod, "main")
