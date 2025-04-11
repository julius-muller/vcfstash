# Add to tests/test_annotation_results.py

import os
from pathlib import Path
import subprocess

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TEST_VCF = os.path.join(TEST_DATA_DIR, "nodata", "crayz_db.bcf")
TEST_CONFIG = os.path.join(os.path.dirname(__file__), "config", "nextflow_test.config")
TEST_ANNO_CONFIG = os.path.join(os.path.dirname(__file__), "config", "annotation.config")
VCFSTASH_CMD = os.path.join(os.path.dirname(os.path.dirname(__file__)), "vcfstash.py")
EXPECTED_OUTPUT_DIR = os.path.join(TEST_DATA_DIR, "expected_output")

def test_cached_vs_uncached_annotation():
    """Test that cached and uncached annotation results match except for headers."""
    cached_bcf = Path(EXPECTED_OUTPUT_DIR) / "annotate_result" / "cached" / "sample4_vst.bcf"
    uncached_bcf = Path(EXPECTED_OUTPUT_DIR) / "annotate_result" / "uncached" / "sample4_vst.bcf"

    assert cached_bcf.exists(), "Cached reference BCF not found"
    assert uncached_bcf.exists(), "Uncached reference BCF not found"

    # Read BCF contents excluding headers
    def get_variants(bcf_path):
        bcf_text = subprocess.run(
            ["bcftools", "view", str(bcf_path)],
            capture_output=True,
            text=True,
            check=True
        ).stdout
        return [line for line in bcf_text.splitlines()
                if line and not line.startswith('#')]

    cached_variants = get_variants(cached_bcf)
    uncached_variants = get_variants(uncached_bcf)

    assert len(cached_variants) > 0, "No variants found in cached BCF"
    assert len(uncached_variants) > 0, "No variants found in uncached BCF"

    # Find and display the first difference
    for i, (cached, uncached) in enumerate(zip(cached_variants, uncached_variants)):
        if cached != uncached:
            raise AssertionError(
                f"First difference at line {i + 1}:\n"
                f"Cached:   {cached}\n"
                f"Uncached: {uncached}"
            )

    # If lengths differ, show which file has extra lines
    if len(cached_variants) != len(uncached_variants):
        raise AssertionError(
            f"Number of variants differs: "
            f"cached={len(cached_variants)}, uncached={len(uncached_variants)}"
        )
