# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

import logging
import warnings
from pathlib import Path

import pytest

from vcfcache.database.outputs import (
    AnnotatedCacheOutput,
    AnnotatedUserOutput,
    CacheOutput,
    MAXCHAR_CACHENAME,
)
from vcfcache.utils.logging import log_command, setup_logging
from vcfcache.utils import validation as v


def test_setup_logging_levels_and_no_file():
    logger = setup_logging(verbosity=0, log_file=None)
    assert logger.level == logging.DEBUG
    assert not any(isinstance(h, logging.FileHandler) for h in logger.handlers)
    console = logger.handlers[0]
    assert console.level == logging.WARNING

    logger = setup_logging(verbosity=1, log_file=None)
    console = logger.handlers[0]
    assert console.level == logging.INFO


def test_setup_logging_with_file(tmp_path):
    log_path = tmp_path / "logs" / "vcfcache.log"
    logger = setup_logging(verbosity=0, log_file=log_path)
    assert any(isinstance(h, logging.FileHandler) for h in logger.handlers)
    file_handler = next(h for h in logger.handlers if isinstance(h, logging.FileHandler))
    assert file_handler.level == logging.INFO
    assert log_path.exists()


def test_log_command_info(monkeypatch, caplog):
    logger = setup_logging(verbosity=1)
    monkeypatch.setattr("sys.argv", ["vcfcache", "demo"])
    with caplog.at_level(logging.INFO, logger="vcfcache"):
        log_command(logger, info=True)
    assert "Script command:" in caplog.text


def test_get_bcf_stats_requires_bcftools():
    with pytest.raises(ValueError):
        v.get_bcf_stats(Path("x.bcf"), bcftools_path=None)


def test_get_bcf_stats_parses_sn(monkeypatch):
    class _Res:
        stdout = "SN\t0\tnumber of records:\t42\nXX\tfoo\n"

    monkeypatch.setattr(v.subprocess, "run", lambda *a, **k: _Res())
    stats = v.get_bcf_stats(Path("x.bcf"), bcftools_path=Path("/usr/bin/bcftools"))
    assert stats["number of records"] == "42"


def test_get_bcf_stats_error(monkeypatch):
    def _run(*_args, **_kwargs):
        raise v.subprocess.CalledProcessError(1, "bcftools")

    monkeypatch.setattr(v.subprocess, "run", _run)
    stats = v.get_bcf_stats(Path("x.bcf"), bcftools_path=Path("/usr/bin/bcftools"))
    assert "error" in stats


def test_validate_bcf_header_requires_bcftools():
    with pytest.raises(ValueError):
        v.validate_bcf_header(Path("x.bcf"), bcftools_path=None)


def test_validate_bcf_header_missing_norm(monkeypatch):
    class _Res:
        stdout = "##fileformat=VCFv4.2\n##contig=<ID=1>\n"

    monkeypatch.setattr(v.subprocess, "run", lambda *a, **k: _Res())
    ok, msg = v.validate_bcf_header(Path("x.bcf"), norm=True, bcftools_path=Path("/usr/bin/bcftools"))
    assert ok is False
    assert "Missing bcftools_normCommand" in msg


def test_validate_bcf_header_missing_options(monkeypatch):
    class _Res:
        stdout = "##bcftools_normCommand=norm -m-\n##contig=<ID=1>\n"

    monkeypatch.setattr(v.subprocess, "run", lambda *a, **k: _Res())
    ok, msg = v.validate_bcf_header(Path("x.bcf"), norm=True, bcftools_path=Path("/usr/bin/bcftools"))
    assert ok is False
    assert "Missing required normalization options" in msg


def test_validate_bcf_header_missing_contigs(monkeypatch):
    class _Res:
        stdout = "##bcftools_normCommand=norm -c x -m-\n"

    monkeypatch.setattr(v.subprocess, "run", lambda *a, **k: _Res())
    ok, msg = v.validate_bcf_header(Path("x.bcf"), norm=True, bcftools_path=Path("/usr/bin/bcftools"))
    assert ok is False
    assert "No contig lines found" in msg


def test_validate_vcf_format_empty(tmp_path):
    vcf = tmp_path / "empty.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    )
    ok, msg = v.validate_vcf_format(vcf)
    assert ok is False
    msg_lower = msg.lower()
    assert "empty" in msg_lower or "error reading vcf file" in msg_lower


def test_validate_vcf_format_missing_file(tmp_path):
    missing = tmp_path / "missing.vcf"
    ok, msg = v.validate_vcf_format(missing)
    assert ok is False
    assert "Error reading VCF file" in msg


def test_compute_md5_missing(tmp_path):
    missing = tmp_path / "nope.bin"
    with pytest.raises(FileNotFoundError):
        v.compute_md5(missing)


def test_outputs_validate_label():
    with pytest.raises(ValueError):
        CacheOutput.validate_label("x" * (MAXCHAR_CACHENAME + 1))
    with pytest.raises(ValueError):
        CacheOutput.validate_label("bad name")
    with pytest.raises(ValueError):
        CacheOutput.validate_label("bad*name")


def test_cache_output_structure(tmp_path):
    out = CacheOutput(str(tmp_path / "cache_root"))
    out.create_structure()
    assert out.validate_structure() is True

    (out.cache_root_dir / "cache").rmdir()
    with warnings.catch_warnings(record=True) as w:
        assert out.validate_structure() is False
        assert any("Missing required path" in str(x.message) for x in w)


def test_annotated_cache_output_structure(tmp_path):
    cache_root = tmp_path / "root"
    cache_out = CacheOutput(str(cache_root))
    cache_out.create_structure()
    init_yaml = cache_out.workflow_dir / "init.yaml"
    init_yaml.parent.mkdir(parents=True, exist_ok=True)
    init_yaml.write_text("params: {}")

    annotation_dir = cache_root / "cache" / "anno1"
    ann = AnnotatedCacheOutput(str(annotation_dir))
    ann.create_structure()
    assert ann.validate_structure() is True

    init_yaml.unlink()
    with warnings.catch_warnings(record=True) as w:
        assert ann.validate_structure() is False
        assert any("Missing required path" in str(x.message) for x in w)


def test_annotated_user_output_structure(tmp_path):
    output_dir = tmp_path / "out_user"
    ann = AnnotatedUserOutput(str(output_dir))
    ann.create_structure()
    assert ann.validate_structure() is True

    (output_dir / "workflow").rmdir()
    with warnings.catch_warnings(record=True) as w:
        assert ann.validate_structure() is False
        assert any("Missing required path" in str(x.message) for x in w)
