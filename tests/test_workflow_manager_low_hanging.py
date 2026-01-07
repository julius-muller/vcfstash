# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

from pathlib import Path

import pytest

from vcfcache.database import workflow_manager as wf


class _DummyLogger:
    def __init__(self):
        self.debugs = []
        self.infos = []
        self.warnings = []
        self.errors = []

    def debug(self, msg):
        self.debugs.append(msg)

    def info(self, msg):
        self.infos.append(msg)

    def warning(self, msg):
        self.warnings.append(msg)

    def error(self, msg):
        self.errors.append(msg)


def _make_manager(tmp_path):
    mgr = wf.WorkflowManager.__new__(wf.WorkflowManager)
    mgr.output_dir = tmp_path
    mgr.output_file = None
    mgr.name = "demo"
    mgr.logger = _DummyLogger()
    mgr.params_file_content = {
        "bcftools_cmd": "bcftools",
        "threads": 2,
        "genome_build": "GRCh38",
    }
    mgr.nfa_config_content = {
        "annotation_cmd": "echo ${INPUT_BCF} > ${OUTPUT_BCF}",
        "must_contain_info_tag": "MOCK",
        "required_tool_version": "1.0",
        "genome_build": "GRCh38",
    }
    mgr.work_dir = None
    return mgr


def test_substitute_variables_env_params_and_extra(monkeypatch, tmp_path):
    mgr = _make_manager(tmp_path)
    mgr.params_file_content = {
        "bcftools_cmd": "/usr/bin/bcftools",
        "vep_cache": "/opt/vep",
        "genome_build": "GRCh38",
    }
    monkeypatch.setenv("VCFCACHE_ROOT", "/data/vcfcache")

    text = "root=${VCFCACHE_ROOT} cache=${params.vep_cache} in=${INPUT_BCF} out=\\${OUTPUT_BCF}"
    out = mgr._substitute_variables(
        text,
        extra_vars={"INPUT_BCF": "in.bcf", "OUTPUT_BCF": "out.bcf"},
    )
    assert "root=/data/vcfcache" in out
    assert "cache=/opt/vep" in out
    assert "in=in.bcf" in out
    assert "out=out.bcf" in out


def test_substitute_variables_skip_var(tmp_path):
    mgr = _make_manager(tmp_path)
    text = "in=${INPUT_BCF} out=${OUTPUT_BCF}"
    out = mgr._substitute_variables(
        text,
        extra_vars={"INPUT_BCF": "in.bcf", "OUTPUT_BCF": "out.bcf"},
        skip_vars=["OUTPUT_BCF"],
    )
    assert "in=in.bcf" in out
    assert "${OUTPUT_BCF}" in out


def test_write_trace_file(tmp_path):
    mgr = _make_manager(tmp_path)
    mgr._write_trace_file("annotate", start_time=wf.datetime.datetime.now(), end_time=wf.datetime.datetime.now())
    trace = tmp_path / "demo_trace.txt"
    assert trace.exists()
    assert "annotate" in trace.read_text()


def test_create_work_dir(tmp_path):
    mgr = _make_manager(tmp_path)
    mgr._create_work_dir(tmp_path, dirname="work")
    assert mgr.work_dir == tmp_path / "work"
    assert mgr.work_dir.exists()

    mgr._create_work_dir(tmp_path, dirname="work")
    assert mgr.work_dir.exists()


def test_cleanup_work_dir(tmp_path):
    mgr = _make_manager(tmp_path)
    mgr._create_work_dir(tmp_path, dirname="work")
    (mgr.work_dir / "x.txt").write_text("x")
    mgr.cleanup_work_dir()
    assert mgr.work_dir is None


def test_validate_info_tag_success(monkeypatch, tmp_path):
    mgr = _make_manager(tmp_path)

    class _Res:
        returncode = 0

    monkeypatch.setattr(wf.subprocess, "run", lambda *a, **k: _Res())
    mgr._validate_info_tag(Path("x.bcf"), "MOCK")


def test_validate_info_tag_failure(monkeypatch, tmp_path):
    mgr = _make_manager(tmp_path)

    class _Res:
        returncode = 1

    monkeypatch.setattr(wf.subprocess, "run", lambda *a, **k: _Res())
    with pytest.raises(RuntimeError):
        mgr._validate_info_tag(Path("x.bcf"), "MOCK")


def test_run_invalid_mode(tmp_path):
    mgr = _make_manager(tmp_path)
    with pytest.raises(ValueError):
        mgr.run("nope")


def test_workflow_manager_genome_mismatch(tmp_path):
    params = tmp_path / "params.yaml"
    params.write_text(
        "annotation_tool_cmd: bcftools\n"
        "bcftools_cmd: bcftools\n"
        "temp_dir: /tmp\n"
        "threads: 1\n"
        "genome_build: GRCh38\n"
    )
    anno = tmp_path / "annotation.yaml"
    anno.write_text(
        "annotation_cmd: echo ${INPUT_BCF} > ${OUTPUT_BCF}\n"
        "must_contain_info_tag: MOCK\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh37\n"
    )

    with pytest.raises(ValueError):
        wf.WorkflowManager(
            input_file=tmp_path / "in.bcf",
            output_dir=tmp_path / "out",
            name="demo",
            params_file=params,
            anno_config_file=anno,
        )


def test_bcftools_command_run_writes_outputs(monkeypatch, tmp_path):
    logger = _DummyLogger()
    cmd = "bcftools view -h input.bcf"

    class _Res:
        returncode = 0
        stdout = "ok"
        stderr = ""

    monkeypatch.setattr(wf.subprocess, "run", lambda *a, **k: _Res())
    bc = wf.BcftoolsCommand(cmd, logger, tmp_path)
    result = bc.run()

    assert result.returncode == 0
    assert (tmp_path / "command.sh").exists()
    assert (tmp_path / "stdout.txt").read_text() == "ok"
    assert (tmp_path / "stderr.txt").read_text() == ""
    assert (tmp_path / "timing.txt").exists()


def test_bcftools_command_run_failure(monkeypatch, tmp_path):
    logger = _DummyLogger()
    cmd = "bcftools view -h input.bcf"

    class _Res:
        returncode = 1
        stdout = ""
        stderr = "boom"

    monkeypatch.setattr(wf.subprocess, "run", lambda *a, **k: _Res())
    bc = wf.BcftoolsCommand(cmd, logger, tmp_path)
    with pytest.raises(wf.subprocess.CalledProcessError):
        bc.run(check=True)


def test_run_blueprint_init_normalize_flag(tmp_path):
    mgr = _make_manager(tmp_path)
    called = {"normalize": None}

    def _run_blueprint_init(normalize=False):
        called["normalize"] = normalize
        return "ok"

    mgr._run_blueprint_init = _run_blueprint_init  # type: ignore[method-assign]
    result = mgr.run("blueprint-init", nextflow_args=["--normalize"])
    assert result == "ok"
    assert called["normalize"] is True


def test_run_requires_db_bcf(tmp_path):
    mgr = _make_manager(tmp_path)
    with pytest.raises(ValueError):
        mgr.run("blueprint-extend")
    with pytest.raises(ValueError):
        mgr.run("cache-build")
    with pytest.raises(ValueError):
        mgr.run("annotate")


def test_run_annotate_nocache_calls_method(tmp_path):
    mgr = _make_manager(tmp_path)
    called = {"ok": False}

    def _run_annotate_nocache(skip_split_multiallelic=False):
        called["ok"] = True
        return "done"

    mgr._run_annotate_nocache = _run_annotate_nocache  # type: ignore[method-assign]
    result = mgr.run("annotate-nocache")
    assert called["ok"] is True
    assert result == "done"


def test_run_blueprint_extend_dispatches(tmp_path):
    mgr = _make_manager(tmp_path)
    called = {"db": None}

    def _run_blueprint_extend(db_bcf, normalize=False):
        called["db"] = db_bcf
        return "ok"

    mgr._run_blueprint_extend = _run_blueprint_extend  # type: ignore[method-assign]
    result = mgr.run("blueprint-extend", db_bcf=Path("db.bcf"))
    assert result == "ok"
    assert called["db"] == Path("db.bcf")


def test_run_cache_build_dispatches(tmp_path):
    mgr = _make_manager(tmp_path)
    called = {"db": None}

    def _run_cache_build(db_bcf):
        called["db"] = db_bcf
        return "ok"

    mgr._run_cache_build = _run_cache_build  # type: ignore[method-assign]
    result = mgr.run("cache-build", db_bcf=Path("db.bcf"))
    assert result == "ok"
    assert called["db"] == Path("db.bcf")


def test_run_annotate_dispatches(tmp_path):
    mgr = _make_manager(tmp_path)
    called = {"db": None, "preserve": None}

    def _run_annotate(db_bcf, preserve_unannotated=False, skip_split_multiallelic=False):
        called["db"] = db_bcf
        called["preserve"] = preserve_unannotated
        return "ok"

    mgr._run_annotate = _run_annotate  # type: ignore[method-assign]
    result = mgr.run("annotate", db_bcf=Path("cache.bcf"), preserve_unannotated=True)
    assert result == "ok"
    assert called["db"] == Path("cache.bcf")
    assert called["preserve"] is True


def test_run_blueprint_init_command_builds(tmp_path, monkeypatch):
    mgr = _make_manager(tmp_path)
    mgr.input_file = tmp_path / "input.bcf"
    mgr.input_file.write_text("bcf")
    mgr.output_dir = tmp_path / "out"
    mgr.output_dir.mkdir()
    mgr.work_dir = tmp_path / "work"
    mgr.work_dir.mkdir()

    calls = {"cmds": []}

    class _Cmd:
        def __init__(self, cmd, logger, work_dir):
            calls["cmds"].append(cmd)

        def run(self, check=True):
            return "ok"

    monkeypatch.setattr(wf, "BcftoolsCommand", _Cmd)

    mgr._run_blueprint_init(normalize=False)
    assert "norm -m-" not in calls["cmds"][0]
    calls["cmds"].clear()
    mgr._run_blueprint_init(normalize=True)
    assert "norm -m-" in calls["cmds"][0]


def test_run_blueprint_extend_commands(tmp_path, monkeypatch):
    mgr = _make_manager(tmp_path)
    mgr.input_file = tmp_path / "input.bcf"
    mgr.input_file.write_text("bcf")
    mgr.output_dir = tmp_path / "out"
    mgr.output_dir.mkdir()
    mgr.work_dir = tmp_path / "work"
    mgr.work_dir.mkdir()

    calls = {"cmds": []}

    class _Cmd:
        def __init__(self, cmd, logger, work_dir):
            calls["cmds"].append(cmd)

        def run(self, check=True):
            return "ok"

    monkeypatch.setattr(wf, "BcftoolsCommand", _Cmd)
    mgr._run_blueprint_extend(Path("db.bcf"), normalize=False)
    assert len(calls["cmds"]) == 2
    assert "merge" in calls["cmds"][1]


def test_run_cache_build_writes_logs(tmp_path, monkeypatch):
    mgr = _make_manager(tmp_path)
    mgr.output_dir = tmp_path / "out"
    mgr.output_dir.mkdir()
    mgr.work_dir = tmp_path / "work"
    mgr.work_dir.mkdir()

    class _Cmd:
        def __init__(self, cmd, logger, work_dir):
            self.work_dir = work_dir

        def run(self, check=True):
            (self.work_dir / "stdout.txt").write_text("anno ok")
            (self.work_dir / "stderr.txt").write_text("")
            return "ok"

    monkeypatch.setattr(wf, "BcftoolsCommand", _Cmd)
    monkeypatch.setattr(mgr, "_validate_info_tag", lambda *_args, **_kwargs: None)

    mgr.nfa_config_content = {
        "annotation_cmd": "echo ${INPUT_BCF} > ${OUTPUT_BCF}",
        "must_contain_info_tag": "MOCK",
        "required_tool_version": "1.0",
        "genome_build": "GRCh38",
    }
    mgr.params_file_content = {"bcftools_cmd": "bcftools", "threads": 1, "genome_build": "GRCh38"}

    mgr._run_cache_build(Path("db.bcf"))
    assert (mgr.output_dir / "annotation_tool.log").exists()


def test_run_annotate_missing_variants(monkeypatch, tmp_path):
    mgr = _make_manager(tmp_path)
    mgr.input_file = tmp_path / "sample.bcf"
    mgr.input_file.write_text("bcf")
    mgr.output_dir = tmp_path / "out"
    mgr.output_dir.mkdir()
    mgr.work_dir = tmp_path / "work"
    mgr.work_dir.mkdir()
    mgr.nfa_config_content = {
        "annotation_cmd": "echo ${INPUT_BCF} > ${OUTPUT_BCF}",
        "must_contain_info_tag": "MOCK",
        "required_tool_version": "1.0",
        "genome_build": "GRCh38",
    }
    mgr.params_file_content = {"bcftools_cmd": "bcftools", "threads": 1, "genome_build": "GRCh38"}

    def _run(cmd, *args, **kwargs):
        class _Res:
            returncode = 0
            stdout = ""
            stderr = ""
        if "index -n" in cmd:
            _Res.stdout = "2\n"
        if cmd.strip().startswith("echo"):
            _Res.stdout = "anno\n"
        return _Res()

    monkeypatch.setattr(wf.subprocess, "run", _run)
    result = mgr._run_annotate(Path("cache.bcf"), preserve_unannotated=True)
    assert result is not None
    assert (mgr.output_dir / "annotation_tool.log").exists()


def test_run_annotate_no_missing(monkeypatch, tmp_path):
    mgr = _make_manager(tmp_path)
    mgr.input_file = tmp_path / "sample.bcf"
    mgr.input_file.write_text("bcf")
    mgr.output_dir = tmp_path / "out"
    mgr.output_dir.mkdir()
    mgr.work_dir = tmp_path / "work"
    mgr.work_dir.mkdir()
    mgr.nfa_config_content = {
        "annotation_cmd": "echo ${INPUT_BCF} > ${OUTPUT_BCF}",
        "must_contain_info_tag": "MOCK",
        "required_tool_version": "1.0",
        "genome_build": "GRCh38",
    }
    mgr.params_file_content = {"bcftools_cmd": "bcftools", "threads": 1, "genome_build": "GRCh38"}

    def _run(cmd, *args, **kwargs):
        class _Res:
            returncode = 0
            stdout = "0\n" if "index -n" in cmd else ""
            stderr = ""
        return _Res()

    monkeypatch.setattr(wf.subprocess, "run", _run)
    result = mgr._run_annotate(Path("cache.bcf"), preserve_unannotated=False)
    assert result is not None
