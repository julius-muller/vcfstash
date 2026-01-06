# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

from pathlib import Path
import subprocess

from vcfcache.database.workflow_base import WorkflowBase


class _DummyWorkflow(WorkflowBase):
    def run(
        self,
        db_mode: str,
        trace: bool = False,
        db_bcf: Path | None = None,
        dag: bool = False,
        timeline: bool = False,
        report: bool = False,
        temp: Path | str = "/tmp",
    ) -> subprocess.CompletedProcess:
        return subprocess.CompletedProcess(args=["dummy"], returncode=0)

    def cleanup_work_dir(self) -> None:
        return None


class _Logger:
    def __init__(self):
        self.warnings: list[str] = []

    def warning(self, msg: str) -> None:
        self.warnings.append(msg)


def test_init_uses_placeholder_workflow(tmp_path):
    wf = _DummyWorkflow(
        input_file=tmp_path / "in.bcf",
        output_dir=tmp_path,
        name="demo",
        workflow=None,
    )
    assert wf.workflow_file == tmp_path / "workflow.stub"
    assert wf.workflow_dir == tmp_path


def test_init_uses_given_workflow_path(tmp_path):
    wf_path = tmp_path / "wf.nf"
    wf = _DummyWorkflow(
        input_file=tmp_path / "in.bcf",
        output_dir=tmp_path / "out",
        name="demo",
        workflow=wf_path,
    )
    assert wf.workflow_file == wf_path
    assert wf.workflow_dir == wf_path.parent


def test_get_temp_files_none(tmp_path):
    wf = _DummyWorkflow(
        input_file=tmp_path / "in.bcf",
        output_dir=tmp_path,
        name="demo",
    )
    wf.work_dir = None
    assert wf._get_temp_files() == []


def test_get_temp_files_present(tmp_path):
    wf = _DummyWorkflow(
        input_file=tmp_path / "in.bcf",
        output_dir=tmp_path,
        name="demo",
    )
    wf.work_dir = tmp_path / "work"
    assert wf._get_temp_files() == [wf.work_dir]


def test_warn_temp_files_no_logger(tmp_path):
    wf = _DummyWorkflow(
        input_file=tmp_path / "in.bcf",
        output_dir=tmp_path,
        name="demo",
    )
    wf.work_dir = tmp_path / "work"
    wf.warn_temp_files()


def test_warn_temp_files_empty(tmp_path):
    wf = _DummyWorkflow(
        input_file=tmp_path / "in.bcf",
        output_dir=tmp_path,
        name="demo",
    )
    wf.logger = _Logger()
    wf.work_dir = None
    wf.warn_temp_files()
    assert wf.logger.warnings == []


def test_warn_temp_files_with_logger_existing_path(tmp_path):
    wf = _DummyWorkflow(
        input_file=tmp_path / "in.bcf",
        output_dir=tmp_path,
        name="demo",
    )
    wf.logger = _Logger()
    wf.work_dir = tmp_path / "work"
    wf.work_dir.mkdir()
    wf.warn_temp_files()
    assert "Temporary files from failed run exist:" in wf.logger.warnings[0]
    assert any(str(wf.work_dir) in msg for msg in wf.logger.warnings)
    assert any("remove these files manually" in msg for msg in wf.logger.warnings)


def test_warn_temp_files_with_logger_missing_path(tmp_path):
    wf = _DummyWorkflow(
        input_file=tmp_path / "in.bcf",
        output_dir=tmp_path,
        name="demo",
    )
    wf.logger = _Logger()
    wf.work_dir = tmp_path / "missing"
    wf.warn_temp_files()
    assert "Temporary files from failed run exist:" in wf.logger.warnings[0]
    assert not any(str(wf.work_dir) in msg for msg in wf.logger.warnings)
    assert any("remove these files manually" in msg for msg in wf.logger.warnings)


def test_call_base_abstract_methods_for_coverage(tmp_path):
    wf = _DummyWorkflow(
        input_file=tmp_path / "in.bcf",
        output_dir=tmp_path,
        name="demo",
    )
    assert WorkflowBase.run(wf, db_mode="annotate") is None
    assert WorkflowBase.cleanup_work_dir(wf) is None
