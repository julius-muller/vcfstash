"""Tests for vcfcache compare command."""

import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from vcfcache.compare import (
    compare_runs,
    parse_workflow_log,
    find_output_bcf,
    format_time,
    read_compare_stats,
)


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary output directory structure."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return output_dir


@pytest.fixture
def completion_flag_uncached(temp_output_dir):
    """Create a completion flag for uncached run."""
    workflow_dir = temp_output_dir / "workflow"
    workflow_dir.mkdir()
    (workflow_dir / "params.snapshot.yaml").write_text("threads: 4\n")
    flag_file = temp_output_dir / ".vcfcache_complete"
    flag_file.write_text(
        "command: annotate\n"
        "mode: uncached\n"
        "version: 0.4.2\n"
        "commit: abc123\n"
        "timestamp: 2026-01-07T10:00:00\n"
    )
    (temp_output_dir / "compare_stats.yaml").write_text(
        "command: annotate\n"
        "mode: uncached\n"
        "output_file: /tmp/out_uncached.bcf\n"
        "input_name: sample.bcf\n"
        "cache_name: demo_cache\n"
        "annotation_yaml_md5: abc123\n"
        "genome_build_params: GRCh38\n"
        "genome_build_annotation: GRCh38\n"
        "vcfcache_version: 0.4.2\n"
        "variant_counts:\n"
        "  total_output: 10\n"
        "  annotated_output: 10\n"
        "  dropped_variants: 0\n"
        "variant_md5:\n"
        "  top10: deadbeef\n"
        "  bottom10: deadbeef\n"
    )
    return temp_output_dir


@pytest.fixture
def completion_flag_cached(tmp_path):
    """Create a completion flag for cached run."""
    output_dir = tmp_path / "cached_output"
    output_dir.mkdir()
    workflow_dir = output_dir / "workflow"
    workflow_dir.mkdir()
    (workflow_dir / "params.snapshot.yaml").write_text("threads: 4\n")
    flag_file = output_dir / ".vcfcache_complete"
    flag_file.write_text(
        "command: annotate\n"
        "mode: cached\n"
        "version: 0.4.2\n"
        "commit: abc123\n"
        "timestamp: 2026-01-07T10:05:00\n"
    )
    (output_dir / "compare_stats.yaml").write_text(
        "command: annotate\n"
        "mode: cached\n"
        "output_file: /tmp/out_cached.bcf\n"
        "input_name: sample.bcf\n"
        "cache_name: demo_cache\n"
        "annotation_yaml_md5: abc123\n"
        "genome_build_params: GRCh38\n"
        "genome_build_annotation: GRCh38\n"
        "vcfcache_version: 0.4.2\n"
        "variant_counts:\n"
        "  total_output: 10\n"
        "  annotated_output: 9\n"
        "  dropped_variants: 1\n"
        "variant_md5:\n"
        "  top10: deadbeef\n"
        "  bottom10: deadbeef\n"
    )
    return output_dir


@pytest.fixture
def workflow_log_with_timing(temp_output_dir):
    """Create a workflow.log file with timing information."""
    workflow_log = temp_output_dir / "workflow.log"
    workflow_log.write_text(
        "[2026-01-06 18:49:03] INFO Starting annotation...\n"
        "[2026-01-06 18:49:36] INFO Command completed in 32.733s: bcftools norm\n"
        "[2026-01-06 18:51:25] INFO Command completed in 109.151s: bcftools annotate\n"
        "[2026-01-06 20:05:36] INFO Workflow completed successfully in 4592.2s\n"
    )
    return temp_output_dir


@pytest.fixture
def output_bcf(temp_output_dir):
    """Create a dummy output BCF file."""
    bcf_file = temp_output_dir / "annotated_sample.bcf"
    bcf_file.write_bytes(b"dummy BCF content for testing")
    return bcf_file


def test_parse_workflow_log(workflow_log_with_timing):
    """Test parsing workflow.log for timing and steps."""
    total_time, steps = parse_workflow_log(workflow_log_with_timing)
    assert total_time == 4592.2
    assert len(steps) == 2
    assert steps[0]["command"] == "bcftools norm"
    assert steps[0]["duration"] == 32.733
    assert steps[1]["command"] == "bcftools annotate"
    assert steps[1]["duration"] == 109.151


def test_parse_workflow_log_missing(tmp_path):
    """Test parsing workflow.log when file doesn't exist."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    total_time, steps = parse_workflow_log(output_dir)
    assert total_time is None
    assert steps == []


def test_format_time():
    """Test time formatting."""
    assert format_time(45.6) == "45.6s"
    assert format_time(125.5) == "2m 5.5s"
    assert format_time(3665.0) == "1h 1m 5.0s"
    assert format_time(7322.5) == "2h 2m 2.5s"


def test_read_compare_stats(tmp_path):
    """Test reading compare_stats.yaml."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    stats_file = output_dir / "compare_stats.yaml"
    stats_file.write_text("mode: cached\nvcfcache_version: 0.4.2\n")

    stats = read_compare_stats(output_dir)
    assert stats["mode"] == "cached"
    assert stats["vcfcache_version"] == "0.4.2"


def test_read_compare_stats_missing(tmp_path):
    """Test reading compare_stats.yaml when file doesn't exist."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    stats = read_compare_stats(output_dir)
    assert stats == {}


def test_find_output_bcf(output_bcf):
    """Test finding output BCF file."""
    output_dir = output_bcf.parent
    found_bcf = find_output_bcf(output_dir)

    assert found_bcf == output_bcf


def test_find_output_bcf_missing(tmp_path):
    """Test finding output BCF when it doesn't exist."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    found_bcf = find_output_bcf(output_dir)
    assert found_bcf is None


def test_compare_runs_missing_directory(tmp_path):
    """Test compare_runs with missing directory."""
    dir1 = tmp_path / "nonexistent"
    dir2 = tmp_path / "also_nonexistent"

    with pytest.raises(FileNotFoundError, match="Directory not found"):
        compare_runs(dir1, dir2)


def test_compare_runs_missing_stats(tmp_path):
    """Test compare_runs with missing stats."""
    dir1 = tmp_path / "dir1"
    dir2 = tmp_path / "dir2"
    dir1.mkdir()
    dir2.mkdir()

    with pytest.raises(ValueError, match="compare_stats.yaml"):
        compare_runs(dir1, dir2)


def test_compare_runs_missing_timing(completion_flag_uncached, completion_flag_cached):
    """Test compare_runs when timing information is missing."""
    compare_runs(completion_flag_uncached, completion_flag_cached)


def test_compare_runs_success(
    completion_flag_uncached,
    completion_flag_cached,
    capsys,
):
    """Test successful comparison of two runs."""
    # Add timing to uncached run
    workflow_log1 = completion_flag_uncached / "workflow.log"
    workflow_log1.write_text(
        "[2026-01-06 20:06:44] INFO Command completed in 100.00s: bcftools norm\n"
        "[2026-01-07 02:58:55] INFO Workflow completed successfully in 150.00s\n"
    )

    # Add timing to cached run
    workflow_log2 = completion_flag_cached / "workflow.log"
    workflow_log2.write_text(
        "[2026-01-06 18:49:03] INFO Command completed in 10.00s: bcftools norm\n"
        "[2026-01-06 18:49:36] INFO Command completed in 20.00s: bcftools annotate\n"
        "[2026-01-06 20:05:36] INFO Workflow completed successfully in 50.00s\n"
    )

    # Run comparison
    compare_runs(completion_flag_uncached, completion_flag_cached)

    # Check output
    captured = capsys.readouterr()
    assert "VCFcache Run Comparison" in captured.out
    assert "Comparator A" in captured.out
    assert "Comparator B" in captured.out
    assert "Top10 MD5" in captured.out


def test_compare_runs_different_outputs(
    completion_flag_uncached,
    completion_flag_cached,
    capsys,
):
    """Test comparison when output files differ."""
    # Add timing to both runs
    workflow_log1 = completion_flag_uncached / "workflow.log"
    workflow_log1.write_text(
        "[2026-01-07 02:58:55] INFO Workflow completed successfully in 150.00s\n"
    )

    workflow_log2 = completion_flag_cached / "workflow.log"
    workflow_log2.write_text(
        "[2026-01-06 20:05:36] INFO Workflow completed successfully in 50.00s\n"
    )

    # Run comparison
    compare_runs(completion_flag_uncached, completion_flag_cached)

    # Check output
    captured = capsys.readouterr()
    assert "VCFcache Run Comparison" in captured.out


def test_compare_runs_same_mode(tmp_path, capsys):
    """Test comparison when both runs are in same mode (no warning - useful for comparing different caches)."""
    # Create two runs both in cached mode
    dir1 = tmp_path / "cached1"
    dir2 = tmp_path / "cached2"
    dir1.mkdir()
    dir2.mkdir()
    (dir1 / "workflow").mkdir()
    (dir2 / "workflow").mkdir()
    (dir1 / "workflow" / "params.snapshot.yaml").write_text("threads: 2\n")
    (dir2 / "workflow" / "params.snapshot.yaml").write_text("threads: 2\n")

    # Create completion flags with same mode
    flag1 = dir1 / ".vcfcache_complete"
    flag1.write_text(
        "command: annotate\n"
        "mode: cached\n"
        "version: 0.4.2\n"
        "commit: abc123\n"
        "timestamp: 2026-01-07T10:00:00\n"
    )

    flag2 = dir2 / ".vcfcache_complete"
    flag2.write_text(
        "command: annotate\n"
        "mode: cached\n"
        "version: 0.4.2\n"
        "commit: abc123\n"
        "timestamp: 2026-01-07T10:05:00\n"
    )
    (dir1 / "compare_stats.yaml").write_text(
        "command: annotate\n"
        "mode: cached\n"
        "output_file: /tmp/out1.bcf\n"
        "input_name: sample.bcf\n"
        "cache_name: demo_cache\n"
        "annotation_yaml_md5: abc123\n"
        "genome_build_params: GRCh38\n"
        "genome_build_annotation: GRCh38\n"
        "vcfcache_version: 0.4.2\n"
        "variant_counts:\n"
        "  total_output: 10\n"
        "  annotated_output: 9\n"
        "  dropped_variants: 1\n"
        "variant_md5:\n"
        "  top10: deadbeef\n"
        "  bottom10: deadbeef\n"
    )
    (dir2 / "compare_stats.yaml").write_text(
        "command: annotate\n"
        "mode: cached\n"
        "output_file: /tmp/out2.bcf\n"
        "input_name: sample.bcf\n"
        "cache_name: demo_cache\n"
        "annotation_yaml_md5: abc123\n"
        "genome_build_params: GRCh38\n"
        "genome_build_annotation: GRCh38\n"
        "vcfcache_version: 0.4.2\n"
        "variant_counts:\n"
        "  total_output: 10\n"
        "  annotated_output: 9\n"
        "  dropped_variants: 1\n"
        "variant_md5:\n"
        "  top10: deadbeef\n"
        "  bottom10: deadbeef\n"
    )

    # Add timing to both
    (dir1 / "workflow.log").write_text(
        "[2026-01-06 20:05:36] INFO Workflow completed successfully in 100.00s\n"
    )
    (dir2 / "workflow.log").write_text(
        "[2026-01-06 20:05:36] INFO Workflow completed successfully in 95.00s\n"
    )

    # Run comparison (should work without warning - useful for comparing different caches)
    compare_runs(dir1, dir2)

    # Check output does NOT contain warning about same mode
    captured = capsys.readouterr()
    assert "VCFcache Run Comparison" in captured.out
    assert "may not be meaningful" not in captured.out


def test_compare_cli_integration(tmp_path):
    """Test that compare command runs end-to-end with CLI flags."""
    # Create two valid output directories with completion flags
    dir1 = tmp_path / "run1"
    dir2 = tmp_path / "run2"
    dir1.mkdir()
    dir2.mkdir()
    (dir1 / "workflow").mkdir()
    (dir2 / "workflow").mkdir()
    (dir1 / "workflow" / "params.snapshot.yaml").write_text("threads: 2\n")
    (dir2 / "workflow" / "params.snapshot.yaml").write_text("threads: 2\n")

    # Create minimal completion flags and stats
    for d, mode in [(dir1, "uncached"), (dir2, "cached")]:
        flag = d / ".vcfcache_complete"
        flag.write_text(
            f"command: annotate\n"
            f"mode: {mode}\n"
            f"version: 0.4.2\n"
            f"commit: test\n"
            f"timestamp: 2026-01-07T10:00:00\n"
        )
        (d / "compare_stats.yaml").write_text(
            "command: annotate\n"
            f"mode: {mode}\n"
            "output_file: /tmp/out.bcf\n"
            "input_name: sample.bcf\n"
            "cache_name: demo_cache\n"
            "annotation_yaml_md5: abc123\n"
            "genome_build_params: GRCh38\n"
            "genome_build_annotation: GRCh38\n"
            "vcfcache_version: 0.4.2\n"
            "variant_counts:\n"
            "  total_output: 10\n"
            "  annotated_output: 10\n"
            "  dropped_variants: 0\n"
            "variant_md5:\n"
            "  top10: deadbeef\n"
            "  bottom10: deadbeef\n"
        )
        # Add timing with proper format
        (d / "workflow.log").write_text(
            "[2026-01-06 20:05:36] INFO Workflow completed successfully in 100.00s\n"
        )
        # Add output file
        (d / "output.bcf").write_bytes(b"test")

    # Test with various flags (should not crash with AttributeError)
    import subprocess as sp

    # Test with --verbose
    result = sp.run(
        ["vcfcache", "compare", str(dir1), str(dir2), "--verbose"],
        capture_output=True,
        text=True,
    )
    # Should not crash with AttributeError about 'verbose'
    assert "AttributeError" not in result.stderr
    assert result.returncode == 0

    # Test with --quiet
    result = sp.run(
        ["vcfcache", "compare", str(dir1), str(dir2), "--quiet"],
        capture_output=True,
        text=True,
    )
    assert "AttributeError" not in result.stderr
    assert result.returncode == 0

    # Test with --debug
    result = sp.run(
        ["vcfcache", "compare", str(dir1), str(dir2), "--debug"],
        capture_output=True,
        text=True,
    )
    assert "AttributeError" not in result.stderr
    assert result.returncode == 0
