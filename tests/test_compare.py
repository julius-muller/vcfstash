"""Tests for vcfcache compare command."""

import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from vcfcache.compare import compare_runs, extract_timing, find_output_bcf, get_md5


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
        "Starting annotation...\n"
        "Processing variants...\n"
        "Command completed in 123.45s\n"
        "Done.\n"
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


def test_extract_timing_from_workflow_log(workflow_log_with_timing):
    """Test extracting timing from workflow.log."""
    timing = extract_timing(workflow_log_with_timing)
    assert timing == 123.45


def test_extract_timing_from_timing_file(tmp_path):
    """Test extracting timing from work/timing.txt when workflow.log exists but has no timing."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    # Create workflow.log without timing info
    workflow_log = output_dir / "workflow.log"
    workflow_log.write_text("Some log without timing\n")

    # Create timing in work/timing.txt
    work_dir = output_dir / "work"
    work_dir.mkdir()
    timing_file = work_dir / "timing.txt"
    timing_file.write_text("98.76")

    timing = extract_timing(output_dir)
    assert timing == 98.76


def test_extract_timing_missing(tmp_path):
    """Test extracting timing when no timing files exist."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    timing = extract_timing(output_dir)
    assert timing is None


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
    workflow_log1.write_text("Command completed in 150.00s\n")

    # Add timing to cached run
    workflow_log2 = completion_flag_cached / "workflow.log"
    workflow_log2.write_text("Command completed in 50.00s\n")

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
    assert "Time saved:" in captured.out
    assert "Output files are identical" in captured.out or "✓" in captured.out


def test_compare_runs_different_outputs(
    completion_flag_uncached,
    completion_flag_cached,
    capsys,
):
    """Test comparison when output files differ."""
    # Add timing to both runs
    workflow_log1 = completion_flag_uncached / "workflow.log"
    workflow_log1.write_text("Command completed in 150.00s\n")

    workflow_log2 = completion_flag_cached / "workflow.log"
    workflow_log2.write_text("Command completed in 50.00s\n")

    # Add output BCF files with DIFFERENT content
    bcf1 = completion_flag_uncached / "output1.bcf"
    bcf2 = completion_flag_cached / "output2.bcf"
    bcf1.write_bytes(b"content A")
    bcf2.write_bytes(b"content B")

    # Run comparison
    compare_runs(completion_flag_uncached, completion_flag_cached)

    # Check output
    captured = capsys.readouterr()
    assert "WARNING: Output files differ" in captured.out or "✗" in captured.out


def test_compare_runs_same_mode(tmp_path, capsys):
    """Test comparison when both runs are in same mode."""
    # Create two uncached runs
    dir1 = tmp_path / "uncached1"
    dir2 = tmp_path / "uncached2"
    dir1.mkdir()
    dir2.mkdir()

    # Create completion flags with same mode
    flag1 = dir1 / ".vcfcache_complete"
    flag1.write_text(
        "command: annotate\n"
        "mode: uncached\n"
        "version: 0.4.2\n"
        "commit: abc123\n"
        "timestamp: 2026-01-07T10:00:00\n"
    )

    flag2 = dir2 / ".vcfcache_complete"
    flag2.write_text(
        "command: annotate\n"
        "mode: uncached\n"
        "version: 0.4.2\n"
        "commit: abc123\n"
        "timestamp: 2026-01-07T10:05:00\n"
    )

    # Add timing to both
    (dir1 / "workflow.log").write_text("Command completed in 100.00s\n")
    (dir2 / "workflow.log").write_text("Command completed in 95.00s\n")

    # Add identical output files
    bcf_content = b"test content"
    (dir1 / "output.bcf").write_bytes(bcf_content)
    (dir2 / "output.bcf").write_bytes(bcf_content)

    # Run comparison
    compare_runs(dir1, dir2)

    # Check output contains warning about same mode
    captured = capsys.readouterr()
    assert "Both runs are in" in captured.out


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
        # Add timing
        (d / "workflow.log").write_text("Command completed in 100.00s\n")
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
