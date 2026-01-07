"""Tests for vcfcache compare command."""

import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from vcfcache.compare import (
    compare_runs,
    parse_workflow_log,
    find_output_bcf,
    get_md5,
    format_time,
    count_variants,
    read_params_yaml,
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
    flag_file = temp_output_dir / ".vcfcache_complete"
    flag_file.write_text(
        "command: annotate\n"
        "mode: uncached\n"
        "version: 0.4.2\n"
        "commit: abc123\n"
        "timestamp: 2026-01-07T10:00:00\n"
    )
    return temp_output_dir


@pytest.fixture
def completion_flag_cached(tmp_path):
    """Create a completion flag for cached run."""
    output_dir = tmp_path / "cached_output"
    output_dir.mkdir()
    flag_file = output_dir / ".vcfcache_complete"
    flag_file.write_text(
        "command: annotate\n"
        "mode: cached\n"
        "version: 0.4.2\n"
        "commit: abc123\n"
        "timestamp: 2026-01-07T10:05:00\n"
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


def test_get_md5(tmp_path):
    """Test MD5 calculation."""
    test_file = tmp_path / "test.txt"
    test_content = b"test content"
    test_file.write_bytes(test_content)

    expected_md5 = hashlib.md5(test_content).hexdigest()
    actual_md5 = get_md5(test_file)

    assert actual_md5 == expected_md5


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


def test_count_variants_missing_file(tmp_path):
    """Test counting variants when file doesn't exist."""
    bcf_file = tmp_path / "nonexistent.bcf"
    count = count_variants(bcf_file)
    assert count is None


def test_read_params_yaml(tmp_path):
    """Test reading params.snapshot.yaml."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    workflow_dir = output_dir / "workflow"
    workflow_dir.mkdir()

    # Create params.snapshot.yaml
    params_file = workflow_dir / "params.snapshot.yaml"
    params_file.write_text(
        "genome_build: GRCh38\n"
        "threads: 8\n"
        "other_param: value\n"
    )

    params = read_params_yaml(output_dir)
    assert params["genome_build"] == "GRCh38"
    assert params["threads"] == 8


def test_read_params_yaml_missing(tmp_path):
    """Test reading params.yaml when file doesn't exist."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    params = read_params_yaml(output_dir)
    assert params == {}


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


def test_compare_runs_missing_completion_flag(tmp_path):
    """Test compare_runs with missing completion flag."""
    dir1 = tmp_path / "dir1"
    dir2 = tmp_path / "dir2"
    dir1.mkdir()
    dir2.mkdir()

    # The error is raised as ValueError with message about no completion flag
    with pytest.raises(ValueError, match="No completion flag found"):
        compare_runs(dir1, dir2)


def test_compare_runs_missing_timing(completion_flag_uncached, completion_flag_cached):
    """Test compare_runs when timing information is missing."""
    # Completion flags exist, but no timing files
    with pytest.raises(ValueError, match="Could not extract timing"):
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

    # Add output BCF files with identical content
    bcf1 = completion_flag_uncached / "output1.bcf"
    bcf2 = completion_flag_cached / "output2.bcf"
    bcf_content = b"test BCF content"
    bcf1.write_bytes(bcf_content)
    bcf2.write_bytes(bcf_content)

    # Run comparison
    compare_runs(completion_flag_uncached, completion_flag_cached)

    # Check output
    captured = capsys.readouterr()
    assert "VCFcache Annotation Comparison" in captured.out
    assert "Speed-up:" in captured.out
    assert "Time Saved:" in captured.out
    assert "Output files are IDENTICAL" in captured.out or "✓" in captured.out


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

    # Add output BCF files with DIFFERENT content
    bcf1 = completion_flag_uncached / "output1.bcf"
    bcf2 = completion_flag_cached / "output2.bcf"
    bcf1.write_bytes(b"content A")
    bcf2.write_bytes(b"content B")

    # Run comparison
    compare_runs(completion_flag_uncached, completion_flag_cached)

    # Check output
    captured = capsys.readouterr()
    assert "WARNING: MD5 checksums DIFFER" in captured.out or "✗" in captured.out


def test_compare_runs_same_mode(tmp_path, capsys):
    """Test comparison when both runs are in same mode (no warning - useful for comparing different caches)."""
    # Create two runs both in cached mode
    dir1 = tmp_path / "cached1"
    dir2 = tmp_path / "cached2"
    dir1.mkdir()
    dir2.mkdir()

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

    # Add timing to both
    (dir1 / "workflow.log").write_text(
        "[2026-01-06 20:05:36] INFO Workflow completed successfully in 100.00s\n"
    )
    (dir2 / "workflow.log").write_text(
        "[2026-01-06 20:05:36] INFO Workflow completed successfully in 95.00s\n"
    )

    # Add identical output files
    bcf_content = b"test content"
    (dir1 / "output.bcf").write_bytes(bcf_content)
    (dir2 / "output.bcf").write_bytes(bcf_content)

    # Run comparison (should work without warning - useful for comparing different caches)
    compare_runs(dir1, dir2)

    # Check output does NOT contain warning about same mode
    captured = capsys.readouterr()
    assert "VCFcache Annotation Comparison" in captured.out
    assert "may not be meaningful" not in captured.out


def test_compare_cli_integration(tmp_path):
    """Test that compare command runs end-to-end with CLI flags."""
    # Create two valid output directories with completion flags
    dir1 = tmp_path / "run1"
    dir2 = tmp_path / "run2"
    dir1.mkdir()
    dir2.mkdir()

    # Create minimal completion flags
    for d, mode in [(dir1, "uncached"), (dir2, "cached")]:
        flag = d / ".vcfcache_complete"
        flag.write_text(
            f"command: annotate\n"
            f"mode: {mode}\n"
            f"version: 0.4.2\n"
            f"commit: test\n"
            f"timestamp: 2026-01-07T10:00:00\n"
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
