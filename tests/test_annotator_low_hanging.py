# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

from pathlib import Path
import subprocess
import sys
import types

import pytest

from vcfcache.database import base as db_base
from vcfcache.database.annotator import DatabaseAnnotator, VCFAnnotator


def test_wavg():
    assert VCFAnnotator.wavg(0.1, 0.2, 10, 5) == pytest.approx((0.1 * 10 + 0.2 * 5) / 15)
    assert VCFAnnotator.wavg(None, None, 0, 0) is None
    assert VCFAnnotator.wavg(None, 0.2, 0, 5) == 0.2
    assert VCFAnnotator.wavg(0.1, None, 5, 0) == 0.1


def test_log_contig_overlap_reports(monkeypatch):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.cache_file = Path("cache.bcf")
    annotator.input_vcf = Path("input.bcf")

    class _Logger:
        def __init__(self):
            self.messages = []

        def info(self, msg):
            self.messages.append(msg)

    annotator.logger = _Logger()

    def _list_contigs(bcf):
        if bcf == annotator.cache_file:
            return ["1", "2", "3"]
        return ["2", "3", "4"]

    annotator._list_contigs = _list_contigs  # type: ignore[method-assign]

    annotator._log_contig_overlap()
    assert any("overlap=2" in msg for msg in annotator.logger.messages)
    assert any("Overlapping contigs" in msg for msg in annotator.logger.messages)


def test_log_contig_overlap_no_overlap(monkeypatch):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.cache_file = Path("cache.bcf")
    annotator.input_vcf = Path("input.bcf")
    annotator.logger = None

    def _list_contigs(bcf):
        if bcf == annotator.cache_file:
            return ["1"]
        return ["2"]

    annotator._list_contigs = _list_contigs  # type: ignore[method-assign]

    with pytest.raises(RuntimeError):
        annotator._log_contig_overlap()


def test_process_variant_annotations():
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    info = {
        "clinvar_clnsig": ["Pathogenic", "Likely"],
        "gnomadg_ac": 10,
        "gnomade_ac": 5,
        "gnomadg_af": 0.1,
        "gnomade_af": 0.2,
    }
    annotator._process_variant_annotations(info)

    assert info["clinvar_clnsig"] == "Pathogenic, Likely"
    assert info["gnomad_af"] == pytest.approx((0.1 * 0.1 + 0.2 * 5) / (0.1 + 5))
    assert "gnomadg_ac" not in info
    assert "gnomade_ac" not in info
    assert "gnomadg_af" not in info
    assert "gnomade_af" not in info


def test_list_contigs_parses_output(monkeypatch, tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.bcftools_path = Path("/usr/bin/bcftools")

    annotator.ensure_indexed = lambda *_a, **_k: None  # type: ignore[assignment]

    class _Res:
        stdout = "1\t100\n2\t200\n\n"

    monkeypatch.setattr("vcfcache.database.annotator.subprocess.run", lambda *a, **k: _Res())
    contigs = annotator._list_contigs(tmp_path / "x.bcf")
    assert contigs == ["1", "2"]


def test_process_region_basic(monkeypatch):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.logger = None
    annotator.BASES = {"A", "C", "G", "T"}
    annotator.INFO_FIELDS = ["GT", "DP", "AF", "clinvar_clnsig"]
    annotator.TRANSCRIPT_KEYS = ["SYMBOL", "PICK"]

    class _InfoEntry:
        description = "Format: SYMBOL|PICK"

    class _Header:
        info = {"CSQ": _InfoEntry()}

    class _VariantFile:
        def __init__(self, _path):
            self.header = _Header()

        def fetch(self, region=None):
            class _Rec:
                chrom = "chr1"
                pos = 10
                ref = "A"
                alts = ["T"]
                info = {"CSQ": ["GENE1|1"], "clinvar_clnsig": ["Pathogenic"]}
                samples = [{"AD": [10, 5], "DP": 15, "GT": (0, 1)}]
            return [_Rec()]

    monkeypatch.setattr("vcfcache.database.annotator.pysam.VariantFile", _VariantFile)
    df = annotator._process_region(("fake.bcf", "chr1:1-100"))
    assert len(df) == 1
    row = df.iloc[0]
    assert row["CHROM"] == "chr1"
    assert row["POS"] == 10
    assert row["REF"] == "A"
    assert row["ALT"] == "T"
    assert row["SYMBOL"] == "GENE1"
    assert bool(row["PICK"]) is True


def test_process_region_import_error(monkeypatch):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.logger = None

    def _importer(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "pandas":
            raise ImportError("no pandas")
        return __import__(name, globals, locals, fromlist, level)

    monkeypatch.setattr("builtins.__import__", _importer)
    with pytest.raises(ImportError):
        annotator._process_region(("fake.bcf", "chr1:1-10"))


def test_convert_to_parquet_success(monkeypatch, tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.logger = None
    annotator.output_dir = tmp_path
    annotator.input_vcf = tmp_path / "sample.bcf"

    class _Header:
        contigs = {"1": None, "2": None}

    class _VF:
        def __init__(self, _path):
            self.header = _Header()

    class _DF:
        def __init__(self, empty=False):
            self.empty = empty

        def __len__(self):
            return 1

        def __len__(self):
            return 1

    class _PD:
        @staticmethod
        def DataFrame(_rows):
            return _DF(empty=False)

        @staticmethod
        def concat(frames, ignore_index=True):
            return _DF(empty=False)

    class _PA:
        class Table:
            @staticmethod
            def from_pandas(_df):
                return "table"

    class _PQ:
        @staticmethod
        def write_table(table, output_file, **_kwargs):
            Path(output_file).write_text("parquet")

    def _process_region(args):
        return _DF(empty=False)

    class _Pool:
        def __init__(self, _threads):
            pass

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def map(self, func, args_list):
            return [func(a) for a in args_list]

    monkeypatch.setattr("vcfcache.database.annotator.pysam.VariantFile", _VF)
    pd_mod = types.ModuleType("pandas")
    pd_mod.DataFrame = _PD.DataFrame
    pd_mod.concat = _PD.concat
    pa_mod = types.ModuleType("pyarrow")
    pa_mod.Table = _PA.Table
    pq_mod = types.ModuleType("pyarrow.parquet")
    pq_mod.write_table = _PQ.write_table
    monkeypatch.setitem(sys.modules, "pandas", pd_mod)
    monkeypatch.setitem(sys.modules, "pyarrow", pa_mod)
    monkeypatch.setitem(sys.modules, "pyarrow.parquet", pq_mod)
    monkeypatch.setattr("vcfcache.database.annotator.Pool", _Pool)
    monkeypatch.setattr(annotator, "_process_region", _process_region)

    out = annotator._convert_to_parquet(tmp_path / "annotated.bcf", threads=2)
    assert out.exists()


def test_convert_to_parquet_no_variants(monkeypatch, tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.logger = None
    annotator.output_dir = tmp_path
    annotator.input_vcf = tmp_path / "sample.bcf"

    class _Header:
        contigs = {"1": None}

    class _VF:
        def __init__(self, _path):
            self.header = _Header()

    class _DF:
        def __init__(self, empty=True):
            self.empty = empty

    class _PD:
        @staticmethod
        def concat(frames, ignore_index=True):
            return _DF(empty=True)

    class _PA:
        class Table:
            @staticmethod
            def from_pandas(_df):
                return "table"

    class _PQ:
        @staticmethod
        def write_table(table, output_file, **_kwargs):
            Path(output_file).write_text("parquet")

    class _Pool:
        def __init__(self, _threads):
            pass

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def map(self, func, args_list):
            return [func(a) for a in args_list]

    monkeypatch.setattr("vcfcache.database.annotator.pysam.VariantFile", _VF)
    pd_mod = types.ModuleType("pandas")
    pd_mod.concat = _PD.concat
    pa_mod = types.ModuleType("pyarrow")
    pa_mod.Table = _PA.Table
    pq_mod = types.ModuleType("pyarrow.parquet")
    pq_mod.write_table = _PQ.write_table
    monkeypatch.setitem(sys.modules, "pandas", pd_mod)
    monkeypatch.setitem(sys.modules, "pyarrow", pa_mod)
    monkeypatch.setitem(sys.modules, "pyarrow.parquet", pq_mod)
    monkeypatch.setattr("vcfcache.database.annotator.Pool", _Pool)
    monkeypatch.setattr(annotator, "_process_region", lambda *_: _DF(empty=True))

    with pytest.raises(ValueError):
        annotator._convert_to_parquet(tmp_path / "annotated.bcf", threads=1)


def test_validate_and_extract_sample_name_vcf_gz_missing_index(tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.input_vcf = tmp_path / "sample.vcf.gz"
    annotator.input_vcf.write_text("vcf")
    with pytest.raises(FileNotFoundError):
        annotator._validate_and_extract_sample_name()


def test_validate_and_extract_sample_name_no_extension(tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.input_vcf = tmp_path / "sample"
    annotator.input_vcf.write_text("vcf")
    with pytest.raises(ValueError):
        annotator._validate_and_extract_sample_name()


def test_validate_and_extract_sample_name_vcf_tbi(tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.input_vcf = tmp_path / "sample.vcf"
    annotator.input_vcf.write_text("vcf")
    (tmp_path / "sample.vcf.tbi").write_text("idx")
    name, ext = annotator._validate_and_extract_sample_name()
    assert name == "sample"
    assert ext == ".vcf"


def test_validate_and_extract_sample_name_vcf_gz_tbi(tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.input_vcf = tmp_path / "sample.vcf.gz"
    annotator.input_vcf.write_text("vcf")
    (tmp_path / "sample.vcf.gz.tbi").write_text("idx")
    name, ext = annotator._validate_and_extract_sample_name()
    assert name == "sample"
    assert ext == ".vcf.gz"


def test_validate_and_extract_sample_name_bcf_csi(tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.input_vcf = tmp_path / "sample.bcf"
    annotator.input_vcf.write_text("bcf")
    (tmp_path / "sample.bcf.csi").write_text("idx")
    name, ext = annotator._validate_and_extract_sample_name()
    assert name == "sample"
    assert ext == ".bcf"


def test_validate_and_extract_sample_name_bad_extension(tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.input_vcf = tmp_path / "sample.txt"
    annotator.input_vcf.write_text("x")
    with pytest.raises(ValueError):
        annotator._validate_and_extract_sample_name()


def test_setup_output_existing_invalid(tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    root = tmp_path / "out"
    root.mkdir()
    annotator.output_annotations = type(
        "Out", (), {"root_dir": root, "validate_structure": lambda *_: False}
    )()
    annotator.logger = None
    with pytest.raises(FileNotFoundError):
        annotator._setup_output(force=False)


def test_validate_inputs_missing_input(tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.logger = None
    annotator.input_vcf = tmp_path / "missing.bcf"
    annotator.ensure_indexed = lambda *_args, **_kwargs: None  # type: ignore[assignment]
    with pytest.raises(FileNotFoundError):
        annotator._validate_inputs()


def test_setup_output_force_removes_existing(tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    root = tmp_path / "out"
    root.mkdir()
    marker = root / "old.txt"
    marker.write_text("old")

    class _Out:
        root_dir = root

        @staticmethod
        def validate_structure():
            return True

        @staticmethod
        def create_structure():
            (root / "workflow").mkdir(parents=True, exist_ok=True)

    annotator.output_annotations = _Out()
    annotator.logger = None

    annotator._setup_output(force=True)
    assert (root / "workflow").exists()


def test_convert_to_parquet_logs(monkeypatch, tmp_path):
    class _Logger:
        def __init__(self):
            self.infos = []
            self.debugs = []

        def info(self, msg):
            self.infos.append(msg)

        def debug(self, msg):
            self.debugs.append(msg)

    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.logger = _Logger()
    annotator.output_dir = tmp_path
    annotator.input_vcf = tmp_path / "sample.bcf"

    class _Header:
        contigs = {"1": None}

    class _VF:
        def __init__(self, _path):
            self.header = _Header()

    class _DF:
        def __init__(self, empty=False):
            self.empty = empty

        def __len__(self):
            return 1

    class _PD:
        @staticmethod
        def concat(frames, ignore_index=True):
            return _DF(empty=False)

    class _PA:
        class Table:
            @staticmethod
            def from_pandas(_df):
                return "table"

    class _PQ:
        @staticmethod
        def write_table(table, output_file, **_kwargs):
            Path(output_file).write_text("parquet")

    class _Pool:
        def __init__(self, _threads):
            pass

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def map(self, func, args_list):
            return [func(a) for a in args_list]

    monkeypatch.setattr("vcfcache.database.annotator.pysam.VariantFile", _VF)
    pd_mod = types.ModuleType("pandas")
    pd_mod.concat = _PD.concat
    pa_mod = types.ModuleType("pyarrow")
    pa_mod.Table = _PA.Table
    pq_mod = types.ModuleType("pyarrow.parquet")
    pq_mod.write_table = _PQ.write_table
    monkeypatch.setitem(sys.modules, "pandas", pd_mod)
    monkeypatch.setitem(sys.modules, "pyarrow", pa_mod)
    monkeypatch.setitem(sys.modules, "pyarrow.parquet", pq_mod)
    monkeypatch.setattr("vcfcache.database.annotator.Pool", _Pool)
    monkeypatch.setattr(annotator, "_process_region", lambda *_: _DF(empty=False))

    out = annotator._convert_to_parquet(tmp_path / "annotated.bcf", threads=1)
    assert out.exists()
    assert any("Converting BCF to Parquet" in msg for msg in annotator.logger.infos)


def test_annotate_calls_cleanup(monkeypatch, tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.logger = None
    annotator.debug = False
    annotator.no_stats = False
    annotator.output_vcf = tmp_path / "out.bcf"
    annotator.cache_file = tmp_path / "cache.bcf"
    annotator.cache_file.write_text("bcf")
    annotator.output_dir = tmp_path
    annotator.output_annotations = type("Out", (), {"root_dir": tmp_path})()
    annotator.nx_workflow = type(
        "WF", (), {"run": lambda *_a, **_k: None, "cleanup_work_dir": lambda *_: None}
    )()
    monkeypatch.setattr(annotator, "_write_compare_stats", lambda *_a, **_k: None)
    annotator._convert_to_parquet = lambda *_args, **_kwargs: tmp_path / "out.parquet"  # type: ignore[assignment]
    annotator.annotate(uncached=True, convert_parquet=True)


def test_annotate_no_stats_skips_outputs(monkeypatch, tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.logger = None
    annotator.debug = False
    annotator.no_stats = True
    annotator.output_vcf = tmp_path / "out.bcf"
    annotator.cache_file = tmp_path / "cache.bcf"
    annotator.cache_file.write_text("bcf")
    annotator.output_dir = tmp_path / "stats"
    annotator.output_dir.mkdir()
    annotator.output_annotations = type("Out", (), {"root_dir": annotator.output_dir})()
    annotator.nx_workflow = type(
        "WF", (), {"run": lambda *_a, **_k: None, "cleanup_work_dir": lambda *_: None}
    )()

    called = {"completion": False, "compare": False}

    def _flag(**_kwargs):
        called["completion"] = True

    monkeypatch.setattr("vcfcache.utils.completion.write_completion_flag", _flag)
    annotator._write_compare_stats = lambda *_a, **_k: called.__setitem__("compare", True)  # type: ignore[assignment]

    annotator.annotate(uncached=True)

    assert not called["completion"]
    assert not called["compare"]
    assert not annotator.output_dir.exists()


def test_annotate_raises_on_workflow_error(tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.logger = None
    annotator.debug = True
    annotator.no_stats = False
    annotator.output_vcf = tmp_path / "out.bcf"
    annotator.cache_file = tmp_path / "cache.bcf"
    annotator.cache_file.write_text("bcf")

    class _WF:
        def run(self, **_kwargs):
            raise RuntimeError("boom")

        def cleanup_work_dir(self):
            pass

    annotator.nx_workflow = _WF()
    with pytest.raises(RuntimeError):
        annotator.annotate()


def test_validate_inputs_success(tmp_path):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.logger = None
    annotator.input_vcf = tmp_path / "sample.bcf"
    annotator.input_vcf.write_text("bcf")
    called = {"indexed": False}

    def _ensure_indexed(_path):
        called["indexed"] = True

    annotator.ensure_indexed = _ensure_indexed  # type: ignore[assignment]
    annotator._validate_inputs()
    assert called["indexed"] is True


def test_process_region_handles_exception(monkeypatch):
    annotator = VCFAnnotator.__new__(VCFAnnotator)
    annotator.logger = None

    class _VF:
        def __init__(self, _path):
            pass

        def fetch(self, region=None):
            raise RuntimeError("bad region")

    monkeypatch.setattr("vcfcache.database.annotator.pysam.VariantFile", _VF)
    with pytest.raises(RuntimeError):
        annotator._process_region(("fake.bcf", "chr1:1-10"))


def test_preprocess_annotation_config(tmp_path):
    annotator = DatabaseAnnotator.__new__(DatabaseAnnotator)
    annotator.output_dir = tmp_path

    class _Logger:
        def __init__(self):
            self.debugs = []

        def debug(self, msg):
            self.debugs.append(msg)

    annotator.logger = _Logger()

    config = tmp_path / "annotation.yaml"
    config.write_text(
        "cmd: $INPUT_BCF ${OUTPUT_BCF} \\\\${INPUT_BCF} \\\\${OUTPUT_BCF} "
        "$AUXILIARY_DIR ${AUXILIARY_DIR}\n"
    )

    output = annotator._preprocess_annotation_config(config)
    text = output.read_text()
    assert "\\$INPUT_BCF" in text
    assert "\\${OUTPUT_BCF}" in text
    assert "\\$AUXILIARY_DIR" in text
    assert "\\\\${" not in text


def test_setup_annotation_cache_invalid_existing(tmp_path):
    annotator = DatabaseAnnotator.__new__(DatabaseAnnotator)

    class _DummyAnno:
        def __init__(self, root: Path):
            self.annotation_dir = root
            self.created = False

        def create_structure(self):
            self.annotation_dir.mkdir(parents=True, exist_ok=True)
            self.created = True

    class _DummyCache:
        @staticmethod
        def validate_structure():
            return False

    root = tmp_path / "cache" / "anno"
    root.mkdir(parents=True)
    annotator.cached_annotations = _DummyAnno(root)
    annotator.cached_output = _DummyCache()
    annotator.logger = None

    with pytest.raises(FileNotFoundError):
        annotator._setup_annotation_cache(force=False)


def test_database_annotator_annotate_error_cleans_output(monkeypatch, tmp_path):
    annotator = DatabaseAnnotator.__new__(DatabaseAnnotator)
    annotator.output_dir = tmp_path / "out"
    annotator.output_dir.mkdir(parents=True)
    annotator.info_file = tmp_path / "sources.info"
    annotator.info_file.write_text("info")
    annotator.info_snapshot_file = annotator.output_dir / "blueprint_snapshot.info"
    annotator.blueprint_bcf = tmp_path / "blueprint" / "vcfcache.bcf"
    annotator.blueprint_bcf.parent.mkdir(parents=True, exist_ok=True)
    annotator.blueprint_bcf.write_text("bcf")

    class _Logger:
        def error(self, *_args, **_kwargs):
            pass

        def warning(self, *_args, **_kwargs):
            pass

        def info(self, *_args, **_kwargs):
            pass

    class _WF:
        def run(self, **_kwargs):
            raise subprocess.CalledProcessError(1, "cmd", stderr="boom")

        def cleanup_work_dir(self):
            pass

    annotator.logger = _Logger()
    annotator.nx_workflow = _WF()
    annotator.debug = False

    monkeypatch.setattr("vcfcache.database.annotator.sys.exit", lambda code=0: (_ for _ in ()).throw(SystemExit(code)))

    with pytest.raises(SystemExit):
        annotator.annotate()

    assert not annotator.output_dir.exists()


def test_database_annotator_init_minimal(monkeypatch, tmp_path):
    db_root = tmp_path / "cache_root"
    (db_root / "cache").mkdir(parents=True)
    (db_root / "blueprint").mkdir()
    (db_root / "workflow").mkdir()

    params_file = tmp_path / "params.yaml"
    params_file.write_text(
        "annotation_tool_cmd: bcftools\n"
        "bcftools_cmd: bcftools\n"
        "temp_dir: /tmp\n"
        "threads: 1\n"
        "genome_build: GRCh38\n"
    )
    anno_config = tmp_path / "annotation.yaml"
    anno_config.write_text(
        "annotation_cmd: $INPUT_BCF\n"
        "must_contain_info_tag: CSQ\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )

    monkeypatch.setattr(db_base, "create_workflow", lambda **_kwargs: type("WF", (), {"run": lambda *_a, **_k: None, "cleanup_work_dir": lambda *_: None})())

    annotator = DatabaseAnnotator(
        annotation_name="anno1",
        db_path=db_root,
        anno_config_file=anno_config,
        bcftools_path=Path("/usr/bin/bcftools"),
        params_file=params_file,
    )

    assert annotator.output_dir.exists()
    assert annotator.params_file.read_text() == params_file.read_text()
    assert annotator.anno_config_file.exists()


def test_vcf_annotator_init_minimal(monkeypatch, tmp_path):
    cache_root = tmp_path / "cache_root"
    (cache_root / "blueprint").mkdir(parents=True)
    (cache_root / "cache").mkdir()
    workflow_dir = cache_root / "workflow"
    workflow_dir.mkdir()
    (workflow_dir / "init.yaml").write_text("threads: 1\n")

    anno_dir = cache_root / "cache" / "anno1"
    anno_dir.mkdir(parents=True)
    (anno_dir / "annotation.yaml").write_text(
        "annotation_cmd: echo ok\n"
        "must_contain_info_tag: CSQ\n"
        "required_tool_version: 1.0\n"
        "genome_build: GRCh38\n"
    )
    (anno_dir / "params.snapshot.yaml").write_text(
        "annotation_tool_cmd: bcftools\n"
        "bcftools_cmd: bcftools\n"
        "temp_dir: /tmp\n"
        "threads: 1\n"
        "genome_build: GRCh38\n"
    )
    (anno_dir / "vcfcache_annotated.bcf").write_text("bcf")

    input_vcf = tmp_path / "sample.bcf"
    input_vcf.write_text("bcf")
    (tmp_path / "sample.bcf.csi").write_text("idx")

    monkeypatch.setattr(db_base, "create_workflow", lambda **_kwargs: type("WF", (), {"run": lambda *_a, **_k: None, "cleanup_work_dir": lambda *_: None})())
    annotator = VCFAnnotator(
        input_vcf=input_vcf,
        annotation_db=anno_dir,
        output_file=tmp_path / "out.bcf",
        stats_dir=tmp_path / "out_stats",
        bcftools_path=Path("/usr/bin/bcftools"),
    )

    assert annotator.output_dir.exists()
    assert annotator.annotation_name == "anno1"
