# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

import json
from pathlib import Path

import pytest
import yaml

from vcfcache.database import base as db_base
from vcfcache.database import initializer as init_mod
from vcfcache.database import updater as upd_mod


class _DummyWorkflow:
    def __init__(self):
        self.ran = False
        self.cleaned = False

    def run(self, **_kwargs):
        self.ran = True

    def cleanup_work_dir(self):
        self.cleaned = True


class _DummyCacheOutput:
    def __init__(self, root_dir: Path, valid: bool = True):
        self.root_dir = root_dir
        self._valid = valid
        self.created = False

    def validate_structure(self) -> bool:
        return self._valid

    def create_structure(self) -> None:
        self.root_dir.mkdir(parents=True, exist_ok=True)
        self.created = True


def test_copy_workflow_srcfiles_skip_config(tmp_path):
    source = tmp_path / "src"
    source.mkdir()
    (source / "a.txt").write_text("x")
    (source / "b.config").write_text("y")

    dest = tmp_path / "dest"
    db_base.VCFDatabase._copy_workflow_srcfiles(source, dest, skip_config=True)

    assert (dest / "a.txt").exists()
    assert not (dest / "b.config").exists()


def test_copy_workflow_srcfiles_source_missing(tmp_path):
    source = tmp_path / "missing"
    dest = tmp_path / "dest"
    db_base.VCFDatabase._copy_workflow_srcfiles(source, dest, skip_config=True)
    assert dest.exists()


def test_initializer_setup_cache_force(tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    root = tmp_path / "cache"
    root.mkdir()
    init.cached_output = _DummyCacheOutput(root, valid=True)

    init._setup_cache(force=True)
    assert init.cached_output.created is True


def test_initializer_setup_cache_existing_no_force(tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    root = tmp_path / "cache"
    root.mkdir()
    init.cached_output = _DummyCacheOutput(root, valid=True)

    with pytest.raises(FileExistsError):
        init._setup_cache(force=False)


def test_initializer_setup_cache_invalid_structure(tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    root = tmp_path / "cache"
    root.mkdir()
    init.cached_output = _DummyCacheOutput(root, valid=False)

    with pytest.raises(FileExistsError):
        init._setup_cache(force=False)


def test_initializer_log_contigs_success(monkeypatch, tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.input_file = tmp_path / "in.bcf"
    init.bcftools_path = Path("/usr/bin/bcftools")
    init.logger = None

    class _Res:
        stdout = "1\t100\n2\t200\nchrM\t30\n"

    monkeypatch.setattr(init_mod.subprocess, "run", lambda *a, **k: _Res())
    init._log_contigs()


def test_initializer_log_contigs_success_with_logger(monkeypatch, tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.input_file = tmp_path / "in.bcf"
    init.bcftools_path = Path("/usr/bin/bcftools")

    class _Logger:
        def __init__(self):
            self.infos = []

        def info(self, msg):
            self.infos.append(msg)

    init.logger = _Logger()

    class _Res:
        stdout = "1\t100\n2\t200\n"

    monkeypatch.setattr(init_mod.subprocess, "run", lambda *a, **k: _Res())
    init._log_contigs()
    assert any("contigs" in msg for msg in init.logger.infos)


def test_initializer_log_contigs_error(monkeypatch, tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.input_file = tmp_path / "in.bcf"
    init.bcftools_path = Path("/usr/bin/bcftools")
    init.logger = None

    def _run(*_args, **_kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(init_mod.subprocess, "run", _run)
    init._log_contigs()


def test_initializer_log_contigs_error_with_logger(monkeypatch, tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.input_file = tmp_path / "in.bcf"
    init.bcftools_path = Path("/usr/bin/bcftools")

    class _Logger:
        def __init__(self):
            self.warnings = []

        def warning(self, msg):
            self.warnings.append(msg)

    init.logger = _Logger()

    def _run(*_args, **_kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(init_mod.subprocess, "run", _run)
    init._log_contigs()
    assert any("Could not list contigs" in msg for msg in init.logger.warnings)


def test_initializer_create_database_success(monkeypatch, tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.logger = None
    init.cached_output = _DummyCacheOutput(tmp_path / "cache", valid=True)
    init.cached_output.create_structure()
    init.cache_name = "demo"
    init.input_file = tmp_path / "input.bcf"
    init.input_file.write_bytes(b"bcf")
    init.info_file = tmp_path / "blueprint" / "sources.info"
    init.info_file.parent.mkdir(parents=True, exist_ok=True)
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.nx_workflow = _DummyWorkflow()
    init.debug = False

    monkeypatch.setattr(init_mod, "compute_md5", lambda *_: "abc123")
    monkeypatch.setattr(init, "_log_contigs", lambda: None)

    init._create_database()
    assert init.info_file.exists()
    data = json.loads(init.info_file.read_text())
    assert data["input_files"][0]["md5"] == "abc123"
    assert init.nx_workflow.ran is True
    assert init.nx_workflow.cleaned is True


def test_initializer_create_database_success_with_logger(monkeypatch, tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)

    class _Logger:
        def __init__(self):
            self.infos = []

        def info(self, msg):
            self.infos.append(msg)

        def error(self, msg):
            self.infos.append(msg)

    init.logger = _Logger()
    init.cached_output = _DummyCacheOutput(tmp_path / "cache", valid=True)
    init.cached_output.create_structure()
    init.cache_name = "demo"
    init.input_file = tmp_path / "input.bcf"
    init.input_file.write_bytes(b"bcf")
    init.info_file = tmp_path / "blueprint" / "sources.info"
    init.info_file.parent.mkdir(parents=True, exist_ok=True)
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.nx_workflow = _DummyWorkflow()
    init.debug = False

    monkeypatch.setattr(init_mod, "compute_md5", lambda *_: "abc123")
    monkeypatch.setattr(init, "_log_contigs", lambda: None)

    init._create_database()
    assert any("Database creation completed successfully" in msg for msg in init.logger.infos)


def test_initializer_create_database_md5_error(monkeypatch, tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.logger = None
    init.cached_output = _DummyCacheOutput(tmp_path / "cache", valid=True)
    init.cached_output.create_structure()
    init.cache_name = "demo"
    init.input_file = tmp_path / "input.bcf"
    init.input_file.write_bytes(b"bcf")
    init.info_file = tmp_path / "blueprint" / "sources.info"
    init.info_file.parent.mkdir(parents=True, exist_ok=True)
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.nx_workflow = _DummyWorkflow()
    init.debug = True

    def _raise(*_args, **_kwargs):
        raise RuntimeError("md5 failed")

    monkeypatch.setattr(init_mod, "compute_md5", _raise)
    monkeypatch.setattr(init, "_log_contigs", lambda: None)

    with pytest.raises(RuntimeError):
        init._create_database()


def test_initializer_create_database_missing_input_file(tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.logger = None
    init.cached_output = _DummyCacheOutput(tmp_path / "cache", valid=True)
    init.cached_output.create_structure()
    init.cache_name = "demo"
    init.input_file = None
    init.info_file = tmp_path / "blueprint" / "sources.info"
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.nx_workflow = _DummyWorkflow()
    init.debug = True

    with pytest.raises(RuntimeError):
        init._create_database()


def test_initializer_initialize_missing_input_logs(tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)

    class _Logger:
        def __init__(self):
            self.errors = []

        def error(self, msg):
            self.errors.append(msg)

        def debug(self, *_args, **_kwargs):
            pass

        def debug(self, *_args, **_kwargs):
            pass

        def debug(self, *_args, **_kwargs):
            pass

    init.logger = _Logger()
    init.input_file = tmp_path / "missing.bcf"
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init._create_database = lambda: None

    with pytest.raises(FileNotFoundError):
        init.initialize()
    assert any("does not exist" in msg for msg in init.logger.errors)


def test_initializer_initialize_existing_output_logs(tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)

    class _Logger:
        def __init__(self):
            self.errors = []

        def error(self, msg):
            self.errors.append(msg)

        def debug(self, *_args, **_kwargs):
            pass

    init.logger = _Logger()
    init.input_file = tmp_path / "input.bcf"
    init.input_file.write_text("bcf")
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.blueprint_bcf.parent.mkdir(parents=True, exist_ok=True)
    init.blueprint_bcf.write_text("bcf")
    init._create_database = lambda: None

    with pytest.raises(FileExistsError):
        init.initialize()
    assert any("Output database already exists" in msg for msg in init.logger.errors)


def test_initializer_validate_inputs_existing_output_logs(tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)

    class _Logger:
        def __init__(self):
            self.errors = []

        def error(self, msg):
            self.errors.append(msg)

        def debug(self, *_args, **_kwargs):
            pass

    init.logger = _Logger()
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.blueprint_bcf.parent.mkdir(parents=True, exist_ok=True)
    init.blueprint_bcf.write_text("bcf")
    init.input_file = tmp_path / "input.bcf"
    init.ensure_indexed = lambda *_a, **_k: None

    with pytest.raises(FileExistsError):
        init._validate_inputs()
    assert any("already exists" in msg for msg in init.logger.errors)


def test_updater_validate_input_files_missing_db(tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = None
    upd.blueprint_bcf = tmp_path / "missing.bcf"
    upd.input_file = tmp_path / "input.bcf"
    with pytest.raises(FileNotFoundError):
        upd._validate_input_files()


def test_updater_validate_inputs_missing_input_logs(tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)

    class _Logger:
        def __init__(self):
            self.errors = []

        def error(self, msg):
            self.errors.append(msg)

        def debug(self, *_args, **_kwargs):
            pass

    upd.logger = _Logger()
    upd.input_file = tmp_path / "missing.bcf"

    with pytest.raises(FileNotFoundError):
        upd._validate_inputs()
    assert any("not found" in msg for msg in upd.logger.errors)


def test_updater_validate_input_files_success_logs(tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)

    class _Logger:
        def __init__(self):
            self.debugs = []

        def debug(self, msg):
            self.debugs.append(msg)

        def error(self, *_args, **_kwargs):
            pass

    upd.logger = _Logger()
    upd.blueprint_bcf = tmp_path / "db.bcf"
    upd.blueprint_bcf.write_text("bcf")
    upd.input_file = tmp_path / "input.bcf"
    upd.input_file.write_text("bcf")
    called = {"count": 0}

    def _ensure_indexed(_path):
        called["count"] += 1

    upd.ensure_indexed = _ensure_indexed  # type: ignore[assignment]
    upd._validate_input_files()
    assert called["count"] == 2
    assert any("Input validation successful" in msg for msg in upd.logger.debugs)


def test_updater_add_logs_file_not_found(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)

    class _Logger:
        def __init__(self):
            self.errors = []

        def info(self, *_args, **_kwargs):
            pass

        def error(self, msg):
            self.errors.append(msg)

    upd.logger = _Logger()
    upd.db_info = {"input_files": []}
    upd.input_md5 = "abc"
    upd.input_file = tmp_path / "input.bcf"

    monkeypatch.setattr(upd_mod, "check_duplicate_md5", lambda **_: False)
    monkeypatch.setattr(upd, "_validate_input_files", lambda: (_ for _ in ()).throw(FileNotFoundError("missing")))

    with pytest.raises(FileNotFoundError):
        upd.add()
    assert any("missing" in msg for msg in upd.logger.errors)


def test_updater_add_logs_value_error(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)

    class _Logger:
        def __init__(self):
            self.errors = []

        def info(self, *_args, **_kwargs):
            pass

        def error(self, msg):
            self.errors.append(msg)

    upd.logger = _Logger()
    upd.db_info = {"input_files": []}
    upd.input_md5 = "abc"
    upd.input_file = tmp_path / "input.bcf"

    monkeypatch.setattr(upd_mod, "check_duplicate_md5", lambda **_: False)
    monkeypatch.setattr(upd, "_validate_input_files", lambda: (_ for _ in ()).throw(ValueError("bad")))

    with pytest.raises(ValueError):
        upd.add()
    assert any("bad" in msg for msg in upd.logger.errors)


def test_updater_merge_variants_logs_summary(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)

    class _Logger:
        def __init__(self):
            self.infos = []
            self.debugs = []

        def info(self, msg):
            self.infos.append(msg)

        def debug(self, msg):
            self.debugs.append(msg)

    upd.logger = _Logger()
    upd.cached_output = _DummyCacheOutput(tmp_path / "cache", valid=True)
    upd.cached_output.create_structure()
    upd.blueprint_bcf = tmp_path / "db.bcf"
    upd.blueprint_bcf.write_bytes(b"bcf")
    upd.input_file = tmp_path / "input.bcf"
    upd.input_md5 = "md5"
    upd.info_file = tmp_path / "sources.info"
    upd.db_info = {"input_files": []}
    upd.bcftools_path = Path("/usr/bin/bcftools")
    upd.nx_workflow = _DummyWorkflow()
    upd.debug = False

    stats = iter(
        [
            {"number of records": "10"},
            {"number of records": "13"},
        ]
    )
    monkeypatch.setattr(upd_mod, "get_bcf_stats", lambda *_args, **_kwargs: next(stats))

    upd._merge_variants()
    assert any("Database update summary" in msg for msg in upd.logger.infos)
def test_updater_add_duplicate_skips(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = None
    upd.db_info = {"input_files": []}
    upd.input_md5 = "abc"
    upd.input_file = tmp_path / "input.bcf"

    called = {"merge": False}

    monkeypatch.setattr(upd_mod, "check_duplicate_md5", lambda **_: True)
    monkeypatch.setattr(upd, "_validate_input_files", lambda: (_ for _ in ()).throw(AssertionError))
    monkeypatch.setattr(upd, "_merge_variants", lambda: called.__setitem__("merge", True))

    upd.add()
    assert called["merge"] is False


def test_updater_merge_variants_updates_info(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = None
    upd.cached_output = _DummyCacheOutput(tmp_path / "cache", valid=True)
    upd.cached_output.create_structure()
    upd.blueprint_bcf = tmp_path / "db.bcf"
    upd.blueprint_bcf.write_bytes(b"bcf")
    upd.input_file = tmp_path / "input.bcf"
    upd.input_md5 = "md5"
    upd.info_file = tmp_path / "sources.info"
    upd.db_info = {"input_files": []}
    upd.bcftools_path = Path("/usr/bin/bcftools")
    upd.nx_workflow = _DummyWorkflow()
    upd.debug = False

    stats = iter(
        [
            {"number of records": "10", "foo": "NA"},
            {"number of records": "13", "foo": "NA"},
        ]
    )
    monkeypatch.setattr(upd_mod, "get_bcf_stats", lambda *_args, **_kwargs: next(stats))

    upd._merge_variants()
    assert upd.info_file.exists()
    data = json.loads(upd.info_file.read_text())
    assert data["input_files"][0]["md5"] == "md5"
    assert upd.nx_workflow.ran is True
    assert upd.nx_workflow.cleaned is True


def test_updater_merge_variants_failure(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = None
    upd.blueprint_bcf = tmp_path / "db.bcf"
    upd.blueprint_bcf.write_text("bcf")
    upd.input_file = tmp_path / "input.bcf"
    upd.input_md5 = "md5"
    upd.info_file = tmp_path / "sources.info"
    upd.db_info = {"input_files": []}
    upd.bcftools_path = Path("/usr/bin/bcftools")

    class _WF:
        def run(self, **_kwargs):
            raise upd_mod.subprocess.CalledProcessError(1, "bcftools")

        def cleanup_work_dir(self):
            pass

    upd.nx_workflow = _WF()
    monkeypatch.setattr(upd_mod, "get_bcf_stats", lambda *_args, **_kwargs: {"number of records": "1"})

    with pytest.raises(RuntimeError):
        upd._merge_variants()


def test_initializer_init_creates_minimal_params(monkeypatch, tmp_path):
    input_file = tmp_path / "input.bcf"
    input_file.write_text("bcf")
    (tmp_path / "input.bcf.csi").write_text("idx")
    output_dir = tmp_path / "cache"

    created = {}

    def _create_workflow(**kwargs):
        created.update(kwargs)
        return _DummyWorkflow()

    monkeypatch.setattr(db_base, "create_workflow", _create_workflow)

    init = init_mod.DatabaseInitializer(
        input_file=input_file,
        bcftools_path=Path("/usr/bin/bcftools"),
        params_file=None,
        output_dir=output_dir,
        verbosity=1,
        force=False,
        debug=False,
        normalize=False,
    )

    data = yaml.safe_load(init.config_yaml.read_text())
    assert data["bcftools_cmd"] == "/usr/bin/bcftools"
    assert data["genome_build"] == "UNKNOWN"
    assert created["input_file"] == input_file


def test_initializer_init_uses_params_file(monkeypatch, tmp_path):
    input_file = tmp_path / "input.bcf"
    input_file.write_text("bcf")
    (tmp_path / "input.bcf.csi").write_text("idx")
    output_dir = tmp_path / "cache"
    params_file = tmp_path / "params.yaml"
    params_file.write_text(
        "annotation_tool_cmd: bcftools\n"
        "bcftools_cmd: bcftools\n"
        "temp_dir: /tmp\n"
        "threads: 4\n"
        "genome_build: GRCh38\n"
    )

    created = {}

    def _create_workflow(**kwargs):
        created.update(kwargs)
        return _DummyWorkflow()

    monkeypatch.setattr(db_base, "create_workflow", _create_workflow)

    init = init_mod.DatabaseInitializer(
        input_file=input_file,
        bcftools_path=Path("/usr/bin/bcftools"),
        params_file=params_file,
        output_dir=output_dir,
        verbosity=1,
        force=False,
        debug=False,
        normalize=False,
    )

    assert init.config_yaml.read_text() == params_file.read_text()
    assert created["params_file"] == params_file


def test_updater_init_with_params_file(monkeypatch, tmp_path):
    db_path = tmp_path / "cache"
    (db_path / "blueprint").mkdir(parents=True)
    (db_path / "cache").mkdir()
    (db_path / "workflow").mkdir()

    input_file = tmp_path / "input.bcf"
    input_file.write_text("bcf")
    (tmp_path / "input.bcf.csi").write_text("idx")
    params_file = tmp_path / "params.yaml"
    params_file.write_text(
        "annotation_tool_cmd: bcftools\n"
        "bcftools_cmd: bcftools\n"
        "temp_dir: /tmp\n"
        "threads: 2\n"
        "genome_build: GRCh38\n"
    )

    monkeypatch.setattr(upd_mod, "compute_md5", lambda *_: "abc")
    monkeypatch.setattr(db_base, "create_workflow", lambda **_kwargs: _DummyWorkflow())

    upd = upd_mod.DatabaseUpdater(
        db_path=db_path,
        input_file=input_file,
        bcftools_path=Path("/usr/bin/bcftools"),
        params_file=params_file,
    )

    expected = db_path / "blueprint" / "add_abc.yaml"
    assert upd.params_file == expected
    assert expected.read_text() == params_file.read_text()


def test_updater_init_uses_init_yaml(monkeypatch, tmp_path):
    db_path = tmp_path / "cache"
    (db_path / "blueprint").mkdir(parents=True)
    (db_path / "cache").mkdir()
    workflow_dir = db_path / "workflow"
    workflow_dir.mkdir()
    init_yaml = workflow_dir / "init.yaml"
    init_yaml.write_text(
        "annotation_tool_cmd: bcftools\n"
        "bcftools_cmd: bcftools\n"
        "temp_dir: /tmp\n"
        "threads: 1\n"
        "genome_build: GRCh38\n"
    )

    input_file = tmp_path / "input.bcf"
    input_file.write_text("bcf")
    (tmp_path / "input.bcf.csi").write_text("idx")

    monkeypatch.setattr(upd_mod, "compute_md5", lambda *_: "abc")
    monkeypatch.setattr(db_base, "create_workflow", lambda **_kwargs: _DummyWorkflow())

    upd = upd_mod.DatabaseUpdater(
        db_path=db_path,
        input_file=input_file,
        bcftools_path=Path("/usr/bin/bcftools"),
        params_file=None,
    )

    assert upd.params_file == init_yaml


def test_updater_init_creates_extend_yaml(monkeypatch, tmp_path):
    db_path = tmp_path / "cache"
    (db_path / "blueprint").mkdir(parents=True)
    (db_path / "cache").mkdir()
    (db_path / "workflow").mkdir()

    input_file = tmp_path / "input.bcf"
    input_file.write_text("bcf")
    (tmp_path / "input.bcf.csi").write_text("idx")

    monkeypatch.setattr(upd_mod, "compute_md5", lambda *_: "abc")
    monkeypatch.setattr(db_base, "create_workflow", lambda **_kwargs: _DummyWorkflow())

    upd = upd_mod.DatabaseUpdater(
        db_path=db_path,
        input_file=input_file,
        bcftools_path=Path("/usr/bin/bcftools"),
        params_file=None,
    )

    assert upd.params_file == db_path / "workflow" / "extend.yaml"
    data = yaml.safe_load(upd.params_file.read_text())
    assert data["bcftools_cmd"] == "/usr/bin/bcftools"
    assert data["genome_build"] == "UNKNOWN"


def test_updater_validate_input_files_ok(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = None
    upd.blueprint_bcf = tmp_path / "db.bcf"
    upd.input_file = tmp_path / "input.bcf"
    upd.blueprint_bcf.write_text("bcf")
    upd.input_file.write_text("bcf")
    (tmp_path / "db.bcf.csi").write_text("idx")
    (tmp_path / "input.bcf.csi").write_text("idx")

    called = {"count": 0}

    def _ensure(path):
        called["count"] += 1

    upd.ensure_indexed = _ensure  # type: ignore[method-assign]
    upd._validate_input_files()
    assert called["count"] == 2


def test_updater_add_calls_merge(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = None
    upd.db_info = {"input_files": []}
    upd.input_md5 = "abc"
    upd.input_file = tmp_path / "input.bcf"

    monkeypatch.setattr(upd_mod, "check_duplicate_md5", lambda **_: False)

    called = {"validated": False, "merged": False}

    upd._validate_input_files = lambda: called.__setitem__("validated", True)  # type: ignore[assignment]
    upd._merge_variants = lambda: called.__setitem__("merged", True)  # type: ignore[assignment]

    upd.add()
    assert called["validated"] is True
    assert called["merged"] is True


def test_initializer_initialize_calls_create(monkeypatch, tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.logger = None
    init.input_file = tmp_path / "input.bcf"
    init.input_file.write_text("bcf")
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.blueprint_bcf.parent.mkdir(parents=True, exist_ok=True)

    called = {"created": False}
    init._create_database = lambda: called.__setitem__("created", True)  # type: ignore[assignment]

    init.initialize()
    assert called["created"] is True


def test_initializer_initialize_missing_input(tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.logger = None
    init.input_file = tmp_path / "missing.bcf"
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.blueprint_bcf.parent.mkdir(parents=True, exist_ok=True)
    with pytest.raises(FileNotFoundError):
        init.initialize()


def test_initializer_initialize_existing_blueprint(tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.logger = None
    init.input_file = tmp_path / "input.bcf"
    init.input_file.write_text("bcf")
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.blueprint_bcf.parent.mkdir(parents=True, exist_ok=True)
    init.blueprint_bcf.write_text("bcf")
    with pytest.raises(FileExistsError):
        init.initialize()


def test_initializer_validate_inputs_blueprint_exists(tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.logger = None
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.blueprint_bcf.parent.mkdir(parents=True, exist_ok=True)
    init.blueprint_bcf.write_text("bcf")
    init.input_file = tmp_path / "input.bcf"
    init.input_file.write_text("bcf")
    init.ensure_indexed = lambda *_args, **_kwargs: None  # type: ignore[assignment]
    with pytest.raises(FileExistsError):
        init._validate_inputs()


def test_initializer_validate_inputs_missing_index(tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.logger = None
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.input_file = tmp_path / "input.bcf"
    init.input_file.write_text("bcf")

    def _ensure(_path):
        raise RuntimeError("no index")

    init.ensure_indexed = _ensure  # type: ignore[assignment]
    with pytest.raises(RuntimeError):
        init._validate_inputs()


def test_initializer_create_database_missing_info_file(monkeypatch, tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.logger = None
    init.cached_output = _DummyCacheOutput(tmp_path / "cache", valid=True)
    init.cached_output.create_structure()
    init.cache_name = "demo"
    init.input_file = tmp_path / "input.bcf"
    init.input_file.write_bytes(b"bcf")
    init.info_file = tmp_path / "blueprint" / "sources.info"
    init.info_file.parent.mkdir(parents=True, exist_ok=True)
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.nx_workflow = _DummyWorkflow()
    init.debug = True

    monkeypatch.setattr(init_mod, "compute_md5", lambda *_: "abc123")
    monkeypatch.setattr(init, "_log_contigs", lambda: None)

    init._create_database()
    assert init.info_file.exists()


def test_initializer_create_database_workflow_error(monkeypatch, tmp_path):
    init = init_mod.DatabaseInitializer.__new__(init_mod.DatabaseInitializer)
    init.logger = None
    init.cached_output = _DummyCacheOutput(tmp_path / "cache", valid=True)
    init.cached_output.create_structure()
    init.cache_name = "demo"
    init.input_file = tmp_path / "input.bcf"
    init.input_file.write_bytes(b"bcf")
    init.info_file = tmp_path / "blueprint" / "sources.info"
    init.info_file.parent.mkdir(parents=True, exist_ok=True)
    init.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    init.debug = True

    class _WF:
        def run(self, **_kwargs):
            raise RuntimeError("workflow failed")

        def cleanup_work_dir(self):
            pass

    init.nx_workflow = _WF()

    monkeypatch.setattr(init_mod, "compute_md5", lambda *_: "abc123")
    monkeypatch.setattr(init, "_log_contigs", lambda: None)

    with pytest.raises(RuntimeError):
        init._create_database()


def test_updater_validate_input_files_missing_input(tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = None
    upd.blueprint_bcf = tmp_path / "db.bcf"
    upd.blueprint_bcf.write_text("bcf")
    upd.input_file = tmp_path / "missing.bcf"
    with pytest.raises(FileNotFoundError):
        upd._validate_input_files()


def test_updater_merge_variants_non_int_stats(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = None
    upd.cached_output = _DummyCacheOutput(tmp_path / "cache", valid=True)
    upd.cached_output.create_structure()
    upd.blueprint_bcf = tmp_path / "db.bcf"
    upd.blueprint_bcf.write_text("bcf")
    upd.input_file = tmp_path / "input.bcf"
    upd.input_md5 = "md5"
    upd.info_file = tmp_path / "sources.info"
    upd.db_info = {"input_files": []}
    upd.bcftools_path = Path("/usr/bin/bcftools")
    upd.nx_workflow = type("WF", (), {"run": lambda *_a, **_k: None, "cleanup_work_dir": lambda *_: None})()
    upd.debug = True

    stats = iter([{"number of records": "NA"}, {"number of records": "NA"}])
    monkeypatch.setattr(upd_mod, "get_bcf_stats", lambda *_args, **_kwargs: next(stats))

    upd._merge_variants()
    assert upd.info_file.exists()


def test_updater_merge_variants_logs_drop(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = None
    upd.cached_output = _DummyCacheOutput(tmp_path / "cache", valid=True)
    upd.cached_output.create_structure()
    upd.blueprint_bcf = tmp_path / "db.bcf"
    upd.blueprint_bcf.write_text("bcf")
    upd.input_file = tmp_path / "input.bcf"
    upd.input_md5 = "md5"
    upd.info_file = tmp_path / "sources.info"
    upd.db_info = {"input_files": []}
    upd.bcftools_path = Path("/usr/bin/bcftools")
    upd.nx_workflow = type("WF", (), {"run": lambda *_a, **_k: None, "cleanup_work_dir": lambda *_: None})()
    upd.debug = False

    stats = iter([{"number of records": "1"}, {"number of records": "3"}])
    monkeypatch.setattr(upd_mod, "get_bcf_stats", lambda *_args, **_kwargs: next(stats))

    upd._merge_variants()
    assert upd.info_file.exists()


def test_updater_merge_variants_stats_missing_keys(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = None
    upd.cached_output = _DummyCacheOutput(tmp_path / "cache", valid=True)
    upd.cached_output.create_structure()
    upd.blueprint_bcf = tmp_path / "db.bcf"
    upd.blueprint_bcf.write_text("bcf")
    upd.input_file = tmp_path / "input.bcf"
    upd.input_md5 = "md5"
    upd.info_file = tmp_path / "sources.info"
    upd.db_info = {"input_files": []}
    upd.bcftools_path = Path("/usr/bin/bcftools")
    upd.nx_workflow = type("WF", (), {"run": lambda *_a, **_k: None, "cleanup_work_dir": lambda *_: None})()
    upd.debug = True

    stats = iter([{"only_pre": "1"}, {"only_post": "2"}])
    monkeypatch.setattr(upd_mod, "get_bcf_stats", lambda *_args, **_kwargs: next(stats))

    upd._merge_variants()
    assert upd.info_file.exists()


def test_updater_add_duplicate_logs(monkeypatch, tmp_path):
    class _Logger:
        def __init__(self):
            self.warnings = []

        def warning(self, msg):
            self.warnings.append(msg)

        def info(self, msg):
            pass

    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = _Logger()
    upd.db_info = {"input_files": []}
    upd.input_md5 = "abc"
    upd.input_file = tmp_path / "input.bcf"

    monkeypatch.setattr(upd_mod, "check_duplicate_md5", lambda **_: True)
    upd.add()
    assert any("duplicate" in msg.lower() for msg in upd.logger.warnings)


def test_updater_add_validate_raises_propagates(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = None
    upd.db_info = {"input_files": []}
    upd.input_md5 = "abc"
    upd.input_file = tmp_path / "input.bcf"

    monkeypatch.setattr(upd_mod, "check_duplicate_md5", lambda **_: False)
    upd._validate_input_files = lambda: (_ for _ in ()).throw(FileNotFoundError("missing"))  # type: ignore[assignment]

    with pytest.raises(FileNotFoundError):
        upd.add()


def test_updater_merge_variants_calledprocesserror(monkeypatch, tmp_path):
    upd = upd_mod.DatabaseUpdater.__new__(upd_mod.DatabaseUpdater)
    upd.logger = None
    upd.blueprint_bcf = tmp_path / "db.bcf"
    upd.blueprint_bcf.write_text("bcf")
    upd.input_file = tmp_path / "input.bcf"
    upd.input_md5 = "md5"
    upd.info_file = tmp_path / "sources.info"
    upd.db_info = {"input_files": []}
    upd.bcftools_path = Path("/usr/bin/bcftools")

    class _WF:
        def run(self, **_kwargs):
            raise upd_mod.subprocess.CalledProcessError(1, "bcftools")

        def cleanup_work_dir(self):
            pass

    upd.nx_workflow = _WF()
    monkeypatch.setattr(upd_mod, "get_bcf_stats", lambda *_args, **_kwargs: {"number of records": "1"})

    with pytest.raises(RuntimeError):
        upd._merge_variants()
