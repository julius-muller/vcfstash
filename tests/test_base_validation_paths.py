# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

import importlib
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from vcfcache.database import base as db_base
from vcfcache.database.workflow_base import WorkflowBase
from vcfcache.utils import validation as v


class _DummyWorkflow(WorkflowBase):
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def run(self, *args, **kwargs):  # pragma: no cover - not used in unit tests
        raise NotImplementedError

    def cleanup_work_dir(self):  # pragma: no cover - not used in unit tests
        raise NotImplementedError


def test_create_workflow_inserts_stub(monkeypatch, tmp_path):
    monkeypatch.setitem(
        sys.modules,
        "vcfcache.database.workflow_manager",
        SimpleNamespace(WorkflowManager=_DummyWorkflow),
    )
    wf = db_base.create_workflow(
        input_file=tmp_path / "input.bcf",
        output_dir=tmp_path,
        name="demo",
        verbosity=0,
    )
    assert isinstance(wf, _DummyWorkflow)
    assert wf.kwargs["workflow"] == tmp_path / "workflow.stub"


def test_create_workflow_requires_output_dir(monkeypatch):
    monkeypatch.setitem(
        sys.modules,
        "vcfcache.database.workflow_manager",
        SimpleNamespace(WorkflowManager=_DummyWorkflow),
    )
    with pytest.raises(ValueError):
        db_base.create_workflow(input_file=Path("input.bcf"))


def test_parse_vcf_info_casting():
    db = db_base.VCFDatabase.__new__(db_base.VCFDatabase)
    db.TRANSCRIPT_KEYS = ["SYMBOL", "PICK", "DISTANCE"]
    data = [
        {"SYMBOL": "GENE1", "PICK": "1", "DISTANCE": "10"},
        {"SYMBOL": "GENE2", "PICK": "0", "DISTANCE": "1.5"},
    ]
    parsed = db.parse_vcf_info(data)
    assert parsed[0]["PICK"] is True
    assert parsed[1]["PICK"] is False
    assert parsed[0]["DISTANCE"] == 10
    assert parsed[1]["DISTANCE"] == 1.5


def test_parse_vcf_info_missing_key():
    db = db_base.VCFDatabase.__new__(db_base.VCFDatabase)
    db.TRANSCRIPT_KEYS = ["SYMBOL", "PICK"]
    with pytest.raises(ValueError):
        db.parse_vcf_info([{"SYMBOL": "GENE"}])


def test_ensure_indexed(tmp_path):
    db = db_base.VCFDatabase.__new__(db_base.VCFDatabase)
    db.logger = None

    bcf = tmp_path / "x.bcf"
    with pytest.raises(FileNotFoundError):
        db.ensure_indexed(bcf)

    bcf.write_bytes(b"bcf")
    with pytest.raises(RuntimeError):
        db.ensure_indexed(bcf)

    (tmp_path / "x.bcf.csi").write_text("idx")
    db.ensure_indexed(bcf)


def test_validate_bcf_header_wrapper(monkeypatch):
    db = db_base.VCFDatabase.__new__(db_base.VCFDatabase)
    db.logger = None
    db.bcftools_path = Path("/usr/bin/bcftools")

    monkeypatch.setattr(
        "vcfcache.database.base.validate_bcf_header",
        lambda *_args, **_kwargs: (False, "bad header"),
    )
    ok, msg = db.validate_bcf_header(Path("x.bcf"), norm=False)
    assert ok is False
    assert msg == "bad header"


def test_check_bcftools_version_ok(monkeypatch):
    class _Result:
        stdout = "1.22\n"

    monkeypatch.setattr(v.subprocess, "run", lambda *a, **k: _Result())
    assert v.check_bcftools_version("/usr/bin/bcftools") == "1.22"


def test_check_bcftools_version_timeout(monkeypatch):
    def _run(*_args, **_kwargs):
        raise v.subprocess.TimeoutExpired(cmd="bcftools", timeout=5)

    monkeypatch.setattr(v.subprocess, "run", _run)
    with pytest.raises(RuntimeError):
        v.check_bcftools_version("/usr/bin/bcftools")


def test_check_bcftools_installed_missing(monkeypatch):
    monkeypatch.setattr(v, "find_bcftools", lambda: None)
    with pytest.raises(FileNotFoundError):
        v.check_bcftools_installed()


def test_check_bcftools_installed_too_old(monkeypatch):
    monkeypatch.setattr(v, "find_bcftools", lambda: "/usr/bin/bcftools")
    monkeypatch.setattr(v, "check_bcftools_version", lambda *_: "1.0")
    with pytest.raises(RuntimeError):
        v.check_bcftools_installed(min_version="1.20")


def test_paths_pkg_root_fallback(monkeypatch):
    import pathlib
    import vcfcache.utils.paths as paths_mod

    orig_exists = pathlib.Path.exists

    def _exists(self):
        if self.name == "pyproject.toml":
            return False
        return orig_exists(self)

    monkeypatch.delenv("VCFCACHE_ROOT", raising=False)
    monkeypatch.setattr(pathlib.Path, "exists", _exists)

    reloaded = importlib.reload(paths_mod)
    root = reloaded.get_project_root()
    assert root == Path(reloaded.__file__).resolve().parent.parent

    monkeypatch.delenv("VCFCACHE_ROOT", raising=False)
    assert reloaded.get_vcfcache_root() == root
