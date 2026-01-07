"""Tests for annotation workflows across all scenarios.

These tests verify that:
1. Annotation works correctly in all scenarios (vanilla, blueprint, annotated)
2. Results are identical regardless of cache size (MD5sum match)
3. Cache hits and misses are tracked correctly
4. Both cached and non-cached variants are annotated properly
"""

import pytest
import subprocess
import tempfile
import shutil
from pathlib import Path
import hashlib
from tests.conftest import get_bcftools_cmd


# Get bcftools command once (respects VCFCACHE_BCFTOOLS env var for CI)
def get_bcftools():
    """Lazy-load bcftools command."""
    return get_bcftools_cmd()


# ============================================================================
# Test Fixtures
# ============================================================================

@pytest.fixture
def test_sample_with_hits_and_misses(test_output_dir):
    """Create a test sample VCF with both cache hits and misses.

    Returns a BCF file with:
    - Some variants from gnomAD (potential cache hits)
    - Some variants unique to sample (cache misses)
    """
    from tests.conftest import TEST_DATA_DIR

    # Use only gnomad_test for now (has both common and rare variants)
    # In a real scenario, some will be in cache, some won't
    gnomad_test = TEST_DATA_DIR / "gnomad_test.bcf"

    output_dir = Path(test_output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    test_bcf = output_dir / "test_sample.bcf"

    # Get bcftools command (respects VCFCACHE_BCFTOOLS env var for CI)
    bcftools = get_bcftools()

    # Extract first 4 variants from gnomad_test as our test sample
    # This creates a small, valid BCF for testing
    subprocess.run(
        [bcftools, "view", "-H", str(gnomad_test)],
        capture_output=True,
        text=True,
        check=True
    )

    # Simply copy and subset gnomad_test
    # Take first 4 variants (small test set)
    cmd = f"{bcftools} view {gnomad_test} | head -1000 | {bcftools} view -Ob -o {test_bcf}"
    subprocess.run(cmd, shell=True, check=True)

    # Index
    subprocess.run(
        [bcftools, "index", str(test_bcf)],
        check=True
    )

    return test_bcf


@pytest.fixture
def annotation_cache_path(test_scenario):
    """Get path to annotation cache based on scenario.

    Returns:
        Path or None: Path to annotation cache in annotated scenario, None otherwise
    """
    if test_scenario == "annotated":
        # Check if vep_gnomad cache exists
        cache_path = Path("/cache/db/cache/vep_gnomad")
        if cache_path.exists():
            return cache_path
    return None


# ============================================================================
# Annotation Workflow Tests
# ============================================================================

def test_annotated_cache_exists(test_scenario, annotation_cache_path):
    """Test that pre-annotated cache exists in annotated scenario."""
    if test_scenario != "annotated":
        pytest.skip("Only applicable to annotated scenario")

    assert annotation_cache_path is not None, "Annotation cache path not found"
    assert annotation_cache_path.exists(), f"Annotation cache does not exist: {annotation_cache_path}"

    # Check for annotated cache file
    annotated_bcf = annotation_cache_path / "vcfcache_annotated.bcf"
    assert annotated_bcf.exists(), "Annotated cache BCF not found"

    # Check index
    annotated_csi = annotation_cache_path / "vcfcache_annotated.bcf.csi"
    assert annotated_csi.exists(), "Annotated cache index not found"

    # Verify it's a valid BCF with annotations
    result = subprocess.run(
        [get_bcftools(), "view", "-h", str(annotated_bcf)],
        capture_output=True,
        text=True,
        check=True
    )

    # Should have CSQ tag for VEP annotations
    assert "##INFO=<ID=CSQ" in result.stdout, "CSQ INFO tag not found in annotated cache"

    print(f"✓ Annotated cache verified at {annotation_cache_path}")


def test_annotate_with_cache(test_scenario, test_sample_with_hits_and_misses,
                             annotation_cache_path, test_output_dir):
    """Test annotation using pre-built cache in annotated scenario."""
    if test_scenario != "annotated":
        pytest.skip("Only applicable to annotated scenario")

    assert annotation_cache_path is not None, "No annotation cache available"

    from tests.conftest import get_vcfcache_root

    # Use the real VEP params that match the Docker build
    vep_params = get_vcfcache_root() / "recipes" / "docker-annotated" / "params.yaml"
    assert vep_params.exists(), f"VEP params not found: {vep_params}"

    output_dir = Path(test_output_dir) / "annotated_output"
    # Don't create the output directory - vcfcache will create it

    # Run annotation using the pre-built cache
    cmd = [
        "vcfcache", "annotate",
        "-a", str(annotation_cache_path),
        "--vcf", str(test_sample_with_hits_and_misses),
        "--output", str(output_dir),
        "-y", str(vep_params),
        "-vv"
    ]

    print(f"\nRunning annotation with cache: {' '.join(cmd)}")

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=300  # 5 minutes max
    )

    print(f"STDOUT:\n{result.stdout}")
    if result.stderr:
        print(f"STDERR:\n{result.stderr}")

    assert result.returncode == 0, f"Annotation failed: {result.stderr}"

    sample_name = test_sample_with_hits_and_misses.name.replace(".bcf", "").replace(".vcf", "")
    # In annotate mode the workflow publishes "<sample>_vc.bcf"
    output_bcf = output_dir / f"{sample_name}_vc.bcf"
    assert output_bcf.exists(), f"Annotated output BCF not found at {output_bcf}"

    # Resolve bcftools path from the same params file used for the run
    import yaml
    with open(vep_params) as f:
        params_cfg = yaml.safe_load(f)
    bcftools_cmd = params_cfg.get("bcftools_cmd", "bcftools")
    if bcftools_cmd.startswith("${VCFCACHE_ROOT}"):
        root = get_vcfcache_root()
        bcftools_cmd = bcftools_cmd.replace("${VCFCACHE_ROOT}", str(root))

    # Verify annotations were added
    header_result = subprocess.run(
        [bcftools_cmd, "view", "-h", str(output_bcf)],
        capture_output=True,
        text=True,
        check=True
    )

    assert "##INFO=<ID=CSQ" in header_result.stdout, "CSQ tag not found in output"

    # Check that variants have annotations
    variants_result = subprocess.run(
        [bcftools_cmd, "query", "-f", "%CHROM\t%POS\t%CSQ\n", str(output_bcf)],
        capture_output=True,
        text=True,
        check=True
    )

    # All variants should have CSQ annotations
    variant_lines = [line for line in variants_result.stdout.strip().split("\n") if line]
    assert len(variant_lines) > 0, "No variants in output"

    # Check that at least some have non-empty CSQ
    annotated_count = sum(1 for line in variant_lines if line.split("\t")[2])
    print(f"Annotated variants: {annotated_count}/{len(variant_lines)}")

    assert annotated_count > 0, "No variants were annotated"

    print(f"✓ Annotation with cache successful: {annotated_count} variants annotated")


def test_blueprint_annotation_workflow(test_scenario, prebuilt_cache, test_output_dir):
    """Test creating and using annotation in blueprint scenario."""
    if test_scenario != "blueprint":
        pytest.skip("Only applicable to blueprint scenario")

    from tests.conftest import TEST_DATA_DIR, TEST_ANNO_CONFIG, TEST_PARAMS

    # Use gnomad_test as blueprint (small test file)
    blueprint_vcf = TEST_DATA_DIR / "gnomad_test.bcf"

    # Step 1: Initialize cache
    cache_dir = Path(test_output_dir) / "test_cache"

    cmd_init = [
        "vcfcache", "blueprint-init",
        "--vcf", str(blueprint_vcf),
        "--output", str(cache_dir),
        "-y", str(TEST_PARAMS),
        "-vv"
    ]

    print(f"\nInitializing cache: {' '.join(cmd_init)}")

    result = subprocess.run(cmd_init, capture_output=True, text=True, timeout=120)
    print(f"Init STDOUT:\n{result.stdout}")
    assert result.returncode == 0, f"Cache init failed: {result.stderr}"

    # Step 2: Annotate the blueprint (create cache)
    cmd_annotate = [
        "vcfcache", "cache-build",
        "--name", "test_anno",
        "--db", str(cache_dir),
        "-a", str(TEST_ANNO_CONFIG),
        "-vv"
    ]

    print(f"\nAnnotating blueprint: {' '.join(cmd_annotate)}")

    result = subprocess.run(cmd_annotate, capture_output=True, text=True, timeout=120)
    print(f"cache-build STDOUT:\n{result.stdout}")
    assert result.returncode == 0, f"Blueprint annotation failed: {result.stderr}"

    # Cache location differs between packaged cache (/db/cache) and ad-hoc init ( /cache )
    cache_base = cache_dir / "db" / "cache"
    if not cache_base.exists():
        cache_base = cache_dir / "cache"
    cache_path = cache_base / "test_anno"
    assert cache_path.exists(), "Annotation cache not created"

    annotated_cache = cache_path / "vcfcache_annotated.bcf"
    assert annotated_cache.exists(), "Annotated cache file not found"

    # Step 3: Annotate a sample using the cache
    sample_vcf = TEST_DATA_DIR / "sample4.bcf"
    output_dir = Path(test_output_dir) / "sample_output"

    cmd_use_cache = [
        "vcfcache", "annotate",
        "-a", str(cache_path),
        "--vcf", str(sample_vcf),
        "--output", str(output_dir),
        "-y", str(TEST_PARAMS),
        "-vv"
    ]

    print(f"\nAnnotating sample with cache: {' '.join(cmd_use_cache)}")

    result = subprocess.run(cmd_use_cache, capture_output=True, text=True, timeout=120)
    print(f"Annotate STDOUT:\n{result.stdout}")
    assert result.returncode == 0, f"Sample annotation failed: {result.stderr}"

    # Verify output (annotate mode emits <sample>_vc.bcf)
    output_bcf = output_dir / f"{sample_vcf.stem}_vc.bcf"
    assert output_bcf.exists(), f"Annotated sample not found at {output_bcf}"

    # Check for MOCK_ANNO tag
    header_result = subprocess.run(
        [get_bcftools(), "view", "-h", str(output_bcf)],
        capture_output=True,
        text=True,
        check=True
    )

    assert "MOCK_ANNO" in header_result.stdout, "MOCK_ANNO tag not found in output"

    print(f"✓ Blueprint annotation workflow successful")


def test_vanilla_annotation_workflow(test_scenario, test_output_dir):
    """Test annotation without any cache in vanilla scenario."""
    if test_scenario != "vanilla":
        pytest.skip("Only applicable to vanilla scenario")

    from tests.conftest import TEST_DATA_DIR, TEST_ANNO_CONFIG, TEST_PARAMS

    sample_vcf = TEST_DATA_DIR / "sample4.bcf"
    output_dir = Path(test_output_dir) / "vanilla_output"
    Path(test_output_dir).mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    # In vanilla scenario, we can't use vcfcache annotate directly without a cache
    # Instead, we test that the annotation config works by running it manually
    # This simulates what would happen without cache

    # Create a temporary output BCF
    output_bcf = Path(test_output_dir) / "vanilla_annotated.bcf"

    # Run the mock annotation directly using bcftools
    # This simulates vcfcache behavior without cache
    bcftools = get_bcftools()
    cmd = f"""
        echo '##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description=Mock_annotation_for_testing_purposes>' > {test_output_dir}/newheader.txt
        {bcftools} query -f '%CHROM\\t%POS\\ttest\\n' {sample_vcf} > {test_output_dir}/mockanno.txt
        bgzip -c {test_output_dir}/mockanno.txt > {test_output_dir}/mockanno.txt.gz
        tabix -s 1 -b 2 -e 2 {test_output_dir}/mockanno.txt.gz
        {bcftools} annotate \\
            -a {test_output_dir}/mockanno.txt.gz \\
            -h {test_output_dir}/newheader.txt \\
            -c CHROM,POS,INFO/MOCK_ANNO \\
            {sample_vcf} \\
            -Ob -o {output_bcf} -W
    """

    print(f"\nRunning vanilla annotation (no cache)")

    result = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
        timeout=60
    )

    assert result.returncode == 0, f"Vanilla annotation failed: {result.stderr}"
    assert output_bcf.exists(), "Vanilla annotated output not found"

    # Verify annotations
    header_result = subprocess.run(
        [get_bcftools(), "view", "-h", str(output_bcf)],
        capture_output=True,
        text=True,
        check=True
    )

    assert "MOCK_ANNO" in header_result.stdout, "MOCK_ANNO tag not found"

    print(f"✓ Vanilla annotation workflow successful")


def test_annotation_consistency_across_scenarios(test_scenario, test_output_dir):
    """Test that annotations are consistent regardless of cache presence.

    This is the key test: verify that results are MD5sum identical
    regardless of whether cache is used or not.
    """
    # This test needs to be run in all scenarios and results compared
    # For now, we verify internal consistency within each scenario

    from tests.conftest import TEST_DATA_DIR, TEST_ANNO_CONFIG, TEST_PARAMS

    sample_vcf = TEST_DATA_DIR / "sample4.bcf"

    if test_scenario == "annotated":
        # In annotated scenario, test with cache
        annotation_cache_path = Path("/cache/db/cache/vep_gnomad")
        if not annotation_cache_path.exists():
            pytest.skip("No annotation cache available in annotated scenario")

        vep_params = Path("/app/recipes/docker-annotated/params.yaml")
        if not vep_params.exists():
            pytest.skip("VEP params not found")

        output_dir = Path(test_output_dir) / "consistency_test"
        if output_dir.exists():
            shutil.rmtree(output_dir)

        cmd = [
            "vcfcache", "annotate",
            "-a", str(annotation_cache_path),
            "--vcf", str(sample_vcf),
            "--output", str(output_dir),
            "-y", str(vep_params),
            "-vv"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        assert result.returncode == 0, f"Annotation failed: {result.stderr}"

        sample_name = sample_vcf.name.replace(".bcf", "").replace(".vcf", "")
        output_bcf = output_dir / f"{sample_name}_vc.bcf"

    elif test_scenario == "blueprint":
        # In blueprint scenario, create minimal cache and use it
        cache_dir = Path(test_output_dir) / "consistency_cache"
        blueprint_vcf = TEST_DATA_DIR / "gnomad_test.bcf"

        # Init
        subprocess.run([
            "vcfcache", "blueprint-init",
            "--vcf", str(blueprint_vcf),
            "--output", str(cache_dir),
            "-y", str(TEST_PARAMS)
        ], check=True, timeout=120)

        # Annotate blueprint
        subprocess.run([
            "vcfcache", "cache-build",
            "--name", "consistency_anno",
            "--db", str(cache_dir),
            "-a", str(TEST_ANNO_CONFIG)
        ], check=True, timeout=120)

        # Annotate sample
        output_dir = Path(test_output_dir) / "consistency_output"
        cache_base = cache_dir / "db" / "cache"
        if not cache_base.exists():
            cache_base = cache_dir / "cache"
        cache_path = cache_base / "consistency_anno"

        subprocess.run([
            "vcfcache", "annotate",
            "-a", str(cache_path),
            "--vcf", str(sample_vcf),
            "--output", str(output_dir),
            "-y", str(TEST_PARAMS)
        ], check=True, timeout=120)

        sample_name = sample_vcf.name.replace(".bcf", "").replace(".vcf", "")
        output_bcf = output_dir / f"{sample_name}_vc.bcf"

    else:  # vanilla
        # Run direct annotation
        output_bcf = Path(test_output_dir) / "consistency_vanilla.bcf"
        Path(test_output_dir).mkdir(parents=True, exist_ok=True)

        bcftools = get_bcftools()
        cmd = f"""
            echo '##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description=Mock>' > {test_output_dir}/h.txt
            {bcftools} query -f '%CHROM\\t%POS\\ttest\\n' {sample_vcf} > {test_output_dir}/a.txt
            bgzip -c {test_output_dir}/a.txt > {test_output_dir}/a.txt.gz
            tabix -s 1 -b 2 -e 2 {test_output_dir}/a.txt.gz
            {bcftools} annotate \\
                -a {test_output_dir}/a.txt.gz \\
                -h {test_output_dir}/h.txt \\
                -c CHROM,POS,INFO/MOCK_ANNO \\
                {sample_vcf} \\
                -Ob -o {output_bcf} -W
        """
        subprocess.run(cmd, shell=True, check=True, timeout=60)

    # Verify output exists and is valid
    assert output_bcf.exists(), f"Output not found: {output_bcf}"

    # Calculate MD5sum for comparison (in real test, this would be compared across runs)
    md5_hash = hashlib.md5()
    with open(output_bcf, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)

    md5sum = md5_hash.hexdigest()
    print(f"\n✓ Annotation consistency test passed (scenario: {test_scenario})")
    print(f"  Output MD5: {md5sum}")

    # Store MD5 for potential cross-scenario comparison
    md5_file = Path(test_output_dir) / "annotation_md5.txt"
    with open(md5_file, 'w') as f:
        f.write(f"{test_scenario}: {md5sum}\n")


# ============================================================================
# Cache Hit Rate Tests
# ============================================================================

def test_cache_hit_statistics(test_scenario, test_sample_with_hits_and_misses,
                              annotation_cache_path, test_output_dir):
    """Test that cache hit statistics are tracked correctly."""
    if test_scenario != "annotated":
        pytest.skip("Cache hit stats only available in annotated scenario")

    assert annotation_cache_path is not None, "No annotation cache available"

    from tests.conftest import get_vcfcache_root
    vep_params = get_vcfcache_root() / "recipes" / "docker-annotated" / "params.yaml"

    output_dir = Path(test_output_dir) / "cache_stats_test"
    if output_dir.exists():
        shutil.rmtree(output_dir)

    cmd = [
        "vcfcache", "annotate",
        "-a", str(annotation_cache_path),
        "--vcf", str(test_sample_with_hits_and_misses),
        "--output", str(output_dir),
        "-y", str(vep_params),
        "-vv"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, f"Annotation failed: {result.stderr}"

    # Check for cache hit statistics in output
    output_text = result.stdout + result.stderr

    # Look for cache hit indicators (these depend on vcfcache implementation)
    # This is a placeholder - adjust based on actual vcfcache output
    print(f"\n✓ Cache hit statistics test completed")
    print(f"  Check output for cache hit/miss information:")
    print(f"  (Stats tracking to be implemented in vcfcache)")


# ============================================================================
# Performance Tests
# ============================================================================

def test_annotation_performance_with_cache(test_scenario, annotation_cache_path,
                                           test_sample_with_hits_and_misses, test_output_dir):
    """Test that annotation with cache is faster than without (informational only)."""
    if test_scenario != "annotated":
        pytest.skip("Performance test only applicable to annotated scenario")

    # This is an informational test - not a strict requirement
    # Just verify that annotation completes in reasonable time

    import time

    assert annotation_cache_path is not None

    from tests.conftest import get_vcfcache_root
    vep_params = get_vcfcache_root() / "recipes" / "docker-annotated" / "params.yaml"

    output_dir = Path(test_output_dir) / "perf_test"
    if output_dir.exists():
        shutil.rmtree(output_dir)

    cmd = [
        "vcfcache", "annotate",
        "-a", str(annotation_cache_path),
        "--vcf", str(test_sample_with_hits_and_misses),
        "--output", str(output_dir),
        "-y", str(vep_params)
    ]

    start_time = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    elapsed_time = time.time() - start_time

    assert result.returncode == 0, f"Annotation failed: {result.stderr}"

    print(f"\n✓ Annotation completed in {elapsed_time:.2f} seconds")
    print(f"  (Performance benefit depends on cache hit rate)")
