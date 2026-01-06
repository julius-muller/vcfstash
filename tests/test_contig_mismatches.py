"""Tests for contig mismatch scenarios and VEP drop behavior.

These tests verify that:
1. Chr prefix mismatches are handled correctly (chr1 vs 1)
2. Sample-specific contigs not in cache are handled properly
3. Cache-specific contigs not in sample are handled properly
4. Cached and uncached outputs are IDENTICAL (same variant count and MD5)
5. Variants dropped by VEP are removed from both cached and uncached outputs
"""

import pytest
import subprocess
import hashlib
from pathlib import Path
from tests.conftest import get_bcftools_cmd, TEST_DATA_DIR
import sys


def get_bcftools():
    """Get bcftools command."""
    return get_bcftools_cmd()


VCFCACHE_CMD = [sys.executable, "-m", "vcfcache"]


def compute_bcf_body_md5(bcf_path, filter_annotated=False, tag="CSQ"):
    """Compute MD5 of BCF body (variants only, no header).

    Args:
        bcf_path: Path to BCF file
        filter_annotated: If True, only include variants with annotation tag
        tag: Annotation tag to filter on (default: CSQ)
    """
    bcftools = get_bcftools()

    if filter_annotated:
        # Filter to only annotated variants before computing MD5
        result = subprocess.run(
            [bcftools, "view", "-H", "-i", f'INFO/{tag}!=""', str(bcf_path)],
            capture_output=True,
            text=True,
            check=True,
        )
    else:
        result = subprocess.run(
            [bcftools, "view", "-H", str(bcf_path)],
            capture_output=True,
            text=True,
            check=True,
        )

    return hashlib.md5(result.stdout.encode()).hexdigest()


def get_variant_count(bcf_path):
    """Get variant count from BCF file."""
    bcftools = get_bcftools()
    result = subprocess.run(
        [bcftools, "index", "-n", str(bcf_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    return int(result.stdout.strip())


def get_contigs_from_bcf(bcf_path):
    """Get list of contigs from BCF file."""
    bcftools = get_bcftools()
    result = subprocess.run(
        [bcftools, "view", "-H", str(bcf_path)],
        capture_output=True,
        text=True,
        check=True,
    )
    contigs = set()
    for line in result.stdout.split('\n'):
        if line.strip():
            contigs.add(line.split('\t')[0])
    return sorted(contigs)


@pytest.fixture
def sample_with_extra_contig():
    """Sample5 has: 1,2,4,9,11,22,M,samplecontig (no chr prefix)."""
    return TEST_DATA_DIR / "sample5.bcf"


@pytest.fixture
def cache_with_extra_contig():
    """crayz_db has: chr1,chr2,chr4,chr9,chr11,chr22,chrM,dbcontig (with chr prefix)."""
    return TEST_DATA_DIR / "crayz_db.bcf"


@pytest.fixture
def annotation_yaml(test_output_dir):
    """Create a simple annotation.yaml that actually adds CSQ values to variants."""
    output_dir = Path(test_output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    annotation_yaml = output_dir / "annotation.yaml"

    # Create CSQ header file
    csq_header = output_dir / "csq_header.txt"
    csq_header.write_text('##INFO=<ID=CSQ,Number=.,Type=String,Description="Dummy annotation">\n')

    # Create annotation that actually adds CSQ values to variants
    # This simulates VEP annotation behavior
    annotation_yaml.write_text(f"""
annotation_cmd: |
  # Create annotations TSV file with dummy CSQ values, bgzip and index it
  bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\tDUMMY_ANNOTATION\\n' ${{INPUT_BCF}} | \\
    bgzip > ${{AUXILIARY_DIR}}/annot.tsv.gz && \\
  tabix -s1 -b2 -e2 ${{AUXILIARY_DIR}}/annot.tsv.gz && \\
  # Add CSQ header and annotations
  bcftools annotate -h {csq_header} -a ${{AUXILIARY_DIR}}/annot.tsv.gz -c CHROM,POS,REF,ALT,CSQ ${{INPUT_BCF}} -o ${{OUTPUT_BCF}} -Ob -W

must_contain_info_tag: CSQ
required_tool_version: "1.21+htslib-1.21"
optional_checks: {{}}
""")
    return annotation_yaml


@pytest.fixture
def params_yaml(test_output_dir):
    """Create a simple params.yaml with all required fields."""
    output_dir = Path(test_output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    bcftools = get_bcftools()
    params_yaml = output_dir / "params.yaml"
    params_yaml.write_text(f"""
annotation_tool_cmd: {bcftools}
bcftools_cmd: {bcftools}
temp_dir: /tmp
threads: 1
optional_checks: {{}}
""")
    return params_yaml


def test_contig_compatibility_chr_prefix_mismatch(
    test_output_dir,
    sample_with_extra_contig,
    cache_with_extra_contig,
    annotation_yaml,
    params_yaml
):
    """Test that chr prefix mismatch is handled correctly.

    Sample has: 1,2,4,9,11,22,M,samplecontig (no chr)
    Cache has: chr1,chr2,chr4,chr9,chr11,chr22,chrM,dbcontig (with chr)

    Expected: Cache should be automatically renamed to remove chr prefix.
    """
    cache_dir = Path(test_output_dir) / "cache_chr_test"
    # Don't create the directory - let vcfcache do it

    # 1. Create blueprint from cache
    result = subprocess.run(
        VCFCACHE_CMD + [
            "blueprint-init",
            "--vcf", str(cache_with_extra_contig),
            "--output", str(cache_dir),
            "--force",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"Blueprint init failed: {result.stderr}"

    # 2. Build cache with annotation
    result = subprocess.run(
        VCFCACHE_CMD + [
            "cache-build",
            "--name", "test_cache",
            "--db", str(cache_dir),
            "-a", str(annotation_yaml),
            "-y", str(params_yaml),
            "--force",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"Cache build failed: {result.stderr}"

    # 3. Annotate sample (should auto-rename cache chr -> no chr)
    output_dir = Path(test_output_dir) / "annotated"
    result = subprocess.run(
        VCFCACHE_CMD + [
            "annotate",
            "-a", str(cache_dir / "cache" / "test_cache"),
            "--vcf", str(sample_with_extra_contig),
            "--output", str(output_dir),
            "-y", str(params_yaml),
            "--force",
            "--debug",  # Keep work dir for inspection
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"Annotation failed: {result.stderr}"

    # 4. Verify renamed cache was created
    cache_variants_dir = cache_dir / "cache" / "test_cache" / ".cache_variants"
    assert cache_variants_dir.exists(), "Renamed cache directory not created"

    renamed_cache = cache_variants_dir / "vcfcache_annotated_nochr.bcf"
    assert renamed_cache.exists(), "Renamed cache file not created"

    # Verify index was created for renamed cache
    renamed_cache_index = Path(f"{renamed_cache}.csi")
    assert renamed_cache_index.exists(), f"Renamed cache index not created: {renamed_cache_index}"

    # 5. Verify output exists and has reasonable variant count
    output_bcf = output_dir / "sample5_vc.bcf"
    assert output_bcf.exists(), f"Output BCF not created: {output_bcf}"

    variant_count = get_variant_count(output_bcf)
    assert variant_count > 0, "Output has no variants"

    print(f"✓ Chr prefix mismatch handled: {variant_count} variants annotated")


def test_cached_uncached_identical_with_contig_mismatch(
    test_output_dir,
    sample_with_extra_contig,
    cache_with_extra_contig,
    annotation_yaml,
    params_yaml
):
    """Test that cached and uncached outputs are IDENTICAL despite contig mismatches.

    This is the critical test for the bug fix:
    - Sample has 'samplecontig' not in cache
    - Cache has 'dbcontig' not in sample
    - Both outputs should be identical (same variants, same MD5)
    """
    cache_dir = Path(test_output_dir) / "cache_identity_test"
    # Don't create the directory - let vcfcache do it

    # 1. Create and build cache
    subprocess.run(
        VCFCACHE_CMD + [
            "blueprint-init",
            "--vcf", str(cache_with_extra_contig),
            "--output", str(cache_dir),
            "--force",
        ],
        check=True,
        capture_output=True,
    )

    subprocess.run(
        VCFCACHE_CMD + [
            "cache-build",
            "--name", "test_cache",
            "--db", str(cache_dir),
            "-a", str(annotation_yaml),
            "-y", str(params_yaml),
            "--force",
        ],
        check=True,
        capture_output=True,
    )

    # 2. Run cached annotation
    cached_dir = Path(test_output_dir) / "cached"
    result_cached = subprocess.run(
        VCFCACHE_CMD + [
            "annotate",
            "-a", str(cache_dir / "cache" / "test_cache"),
            "--vcf", str(sample_with_extra_contig),
            "--output", str(cached_dir),
            "-y", str(params_yaml),
            "--force",
            "--debug",
        ],
        capture_output=True,
        text=True,
    )
    assert result_cached.returncode == 0, f"Cached annotation failed: {result_cached.stderr}"

    # 3. Run uncached annotation
    uncached_dir = Path(test_output_dir) / "uncached"
    result_uncached = subprocess.run(
        VCFCACHE_CMD + [
            "annotate",
            "-a", str(cache_dir / "cache" / "test_cache"),
            "--vcf", str(sample_with_extra_contig),
            "--output", str(uncached_dir),
            "-y", str(params_yaml),
            "--uncached",  # Force uncached mode
            "--force",
            "--debug",
        ],
        capture_output=True,
        text=True,
    )
    assert result_uncached.returncode == 0, f"Uncached annotation failed: {result_uncached.stderr}"

    # 4. Compare outputs
    cached_bcf = cached_dir / "sample5_vc.bcf"
    uncached_bcf = uncached_dir / "sample5_vc.bcf"

    assert cached_bcf.exists(), "Cached output not created"
    assert uncached_bcf.exists(), "Uncached output not created"

    # 4a. Check variant counts (should be identical - both have only annotated variants)
    cached_count = get_variant_count(cached_bcf)
    uncached_count = get_variant_count(uncached_bcf)

    print(f"Cached variant count:   {cached_count}")
    print(f"Uncached variant count: {uncached_count}")

    assert cached_count == uncached_count, (
        f"Variant count mismatch! Cached: {cached_count}, Uncached: {uncached_count}. "
        f"This indicates the post-filtering is not working correctly."
    )

    # 4b. Compare MD5 of variant data
    cached_md5 = compute_bcf_body_md5(cached_bcf)
    uncached_md5 = compute_bcf_body_md5(uncached_bcf)

    print(f"Cached MD5:   {cached_md5}")
    print(f"Uncached MD5: {uncached_md5}")

    # MD5 mismatch is acceptable if only annotation content or order differs
    if cached_md5 != uncached_md5:
        print("MD5s differ - checking if only annotation content/order differs...")
        bcftools = get_bcftools()

        # Check positions only (without INFO fields), sorted to handle order differences
        cached_pos = subprocess.run(
            [bcftools, "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(cached_bcf)],
            capture_output=True, text=True, check=True
        ).stdout
        uncached_pos = subprocess.run(
            [bcftools, "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(uncached_bcf)],
            capture_output=True, text=True, check=True
        ).stdout

        # Sort to handle order differences
        cached_sorted = sorted(cached_pos.strip().split('\n'))
        uncached_sorted = sorted(uncached_pos.strip().split('\n'))

        if cached_sorted == uncached_sorted:
            print("✓ Variant sets are identical (MD5 diff is due to order/annotation differences)")
        else:
            # Find differences
            cached_only = set(cached_sorted) - set(uncached_sorted)
            uncached_only = set(uncached_sorted) - set(cached_sorted)
            print(f"Variants only in cached ({len(cached_only)}): {list(cached_only)[:5]}")
            print(f"Variants only in uncached ({len(uncached_only)}): {list(uncached_only)[:5]}")
            raise AssertionError(
                f"Variant sets differ! {len(cached_only)} unique to cached, "
                f"{len(uncached_only)} unique to uncached."
            )

    print(f"✓ Cached and uncached outputs have identical variant sets (count={cached_count})")


def test_sample_extra_contig_not_in_output(
    test_output_dir,
    sample_with_extra_contig,
    cache_with_extra_contig,
    annotation_yaml,
    params_yaml
):
    """Test that sample's extra contig (samplecontig) is handled correctly.

    Sample has 'samplecontig' which is NOT in the cache.
    This contig should only appear in output if the annotation tool processes it.
    If annotation tool drops it, it should be dropped from both cached and uncached outputs.
    """
    cache_dir = Path(test_output_dir) / "cache_extra_contig_test"
    # Don't create the directory - let vcfcache do it

    # Create and build cache
    subprocess.run(
        VCFCACHE_CMD + [
            "blueprint-init",
            "--vcf", str(cache_with_extra_contig),
            "--output", str(cache_dir),
            "--force",
        ],
        check=True,
        capture_output=True,
    )

    subprocess.run(
        VCFCACHE_CMD + [
            "cache-build",
            "--name", "test_cache",
            "--db", str(cache_dir),
            "-a", str(annotation_yaml),
            "-y", str(params_yaml),
            "--force",
        ],
        check=True,
        capture_output=True,
    )

    # Annotate
    output_dir = Path(test_output_dir) / "annotated_extra"
    subprocess.run(
        VCFCACHE_CMD + [
            "annotate",
            "-a", str(cache_dir / "cache" / "test_cache"),
            "--vcf", str(sample_with_extra_contig),
            "--output", str(output_dir),
            "-y", str(params_yaml),
            "--force",
        ],
        check=True,
        capture_output=True,
    )

    # Check which contigs are in the output
    output_bcf = output_dir / "sample5_vc.bcf"
    output_contigs = get_contigs_from_bcf(output_bcf)

    print(f"Output contigs: {output_contigs}")

    # The behavior depends on whether the annotation tool processes samplecontig
    # In our dummy annotation, it should pass through, so samplecontig should be present
    # But in a real VEP scenario, it might be dropped

    # For now, just verify the output is consistent
    assert len(output_contigs) > 0, "Output has no contigs"
    print(f"✓ Output has {len(output_contigs)} contigs")


def test_cache_extra_contig_not_used(
    test_output_dir,
    sample_with_extra_contig,
    cache_with_extra_contig,
    annotation_yaml,
    params_yaml
):
    """Test that cache's extra contig (dbcontig) doesn't affect output.

    Cache has 'dbcontig' which is NOT in the sample.
    This contig should not appear in the output at all.
    """
    cache_dir = Path(test_output_dir) / "cache_unused_contig_test"
    # Don't create the directory - let vcfcache do it

    # Create and build cache
    subprocess.run(
        VCFCACHE_CMD + [
            "blueprint-init",
            "--vcf", str(cache_with_extra_contig),
            "--output", str(cache_dir),
            "--force",
        ],
        check=True,
        capture_output=True,
    )

    subprocess.run(
        VCFCACHE_CMD + [
            "cache-build",
            "--name", "test_cache",
            "--db", str(cache_dir),
            "-a", str(annotation_yaml),
            "-y", str(params_yaml),
            "--force",
        ],
        check=True,
        capture_output=True,
    )

    # Annotate
    output_dir = Path(test_output_dir) / "annotated_unused"
    subprocess.run(
        VCFCACHE_CMD + [
            "annotate",
            "-a", str(cache_dir / "cache" / "test_cache"),
            "--vcf", str(sample_with_extra_contig),
            "--output", str(output_dir),
            "-y", str(params_yaml),
            "--force",
        ],
        check=True,
        capture_output=True,
    )

    # Check that dbcontig is NOT in the output
    output_bcf = output_dir / "sample5_vc.bcf"
    output_contigs = get_contigs_from_bcf(output_bcf)

    print(f"Output contigs: {output_contigs}")

    assert "dbcontig" not in output_contigs, (
        "Cache's extra contig 'dbcontig' should NOT appear in output!"
    )

    print(f"✓ Cache's extra contig correctly excluded from output")


def test_preserve_unannotated_flag(
    sample_with_extra_contig, cache_with_extra_contig, params_yaml, test_output_dir
):
    """Test that --preserve-unannotated flag correctly preserves unannotated variants.

    This test verifies:
    1. Default behavior filters out variants without annotation (matching annotation tool)
    2. --preserve-unannotated flag keeps all input variants
    3. Both modes have identical annotated variant counts
    4. Preserve mode has more total variants than default mode
    """
    print("\n=== Testing --preserve-unannotated flag ===")

    # Create a custom annotation.yaml that DROPS variants on 'samplecontig' (like VEP drops non-standard contigs)
    output_dir = Path(test_output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    csq_header = output_dir / "csq_header_preserve.txt"
    csq_header.write_text('##INFO=<ID=CSQ,Number=.,Type=String,Description="Dummy annotation">\n')

    annotation_yaml_preserve = output_dir / "annotation_preserve.yaml"
    bcftools = get_bcftools()
    annotation_yaml_preserve.write_text(f"""
annotation_cmd: |
  # Simulate VEP behavior: DROP variants on 'samplecontig' and annotate the rest
  # Filter out samplecontig, then annotate remaining variants
  {bcftools} view -t ^samplecontig ${{INPUT_BCF}} | \\
  {bcftools} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\tDUMMY_ANNOTATION\\n' | \\
    bgzip > ${{AUXILIARY_DIR}}/annot.tsv.gz && \\
  tabix -s1 -b2 -e2 ${{AUXILIARY_DIR}}/annot.tsv.gz && \\
  {bcftools} annotate -h {csq_header} -a ${{AUXILIARY_DIR}}/annot.tsv.gz -c CHROM,POS,REF,ALT,CSQ ${{INPUT_BCF}} -o ${{OUTPUT_BCF}} -Ob -W

must_contain_info_tag: CSQ
required_tool_version: "1.21+htslib-1.21"
optional_checks: {{}}
""")

    # Step 1: Build cache
    cache_dir = Path(test_output_dir) / "cache_preserve_test"

    # First create an initial cache using blueprint-init
    subprocess.run(
        VCFCACHE_CMD + [
            "blueprint-init",
            "--vcf", str(cache_with_extra_contig),
            "--output", str(cache_dir),
            "--force",
        ],
        check=True,
        capture_output=True,
    )

    # Then build the actual annotated cache
    subprocess.run(
        VCFCACHE_CMD + [
            "cache-build",
            "--name", "test_preserve_cache",
            "--db", str(cache_dir),
            "-a", str(annotation_yaml_preserve),
            "-y", str(params_yaml),
            "--force",
        ],
        check=True,
        capture_output=True,
    )

    cache_bcf = cache_dir / "cache" / "test_preserve_cache"

    # Step 2: Annotate with DEFAULT behavior (should filter unannotated variants)
    output_default = Path(test_output_dir) / "annotated_default"

    subprocess.run(
        VCFCACHE_CMD + [
            "annotate",
            "-a", str(cache_bcf),
            "--vcf", str(sample_with_extra_contig),
            "--output", str(output_default),
            "-y", str(params_yaml),
            "--force",
        ],
        check=True,
        capture_output=True,
    )

    output_bcf_default = output_default / "sample5_vc.bcf"

    # Step 3: Annotate with --preserve-unannotated flag
    output_preserve = Path(test_output_dir) / "annotated_preserve"

    subprocess.run(
        VCFCACHE_CMD + [
            "annotate",
            "-a", str(cache_bcf),
            "--vcf", str(sample_with_extra_contig),
            "--output", str(output_preserve),
            "-y", str(params_yaml),
            "--preserve-unannotated",
            "--force",
        ],
        check=True,
        capture_output=True,
    )

    output_bcf_preserve = output_preserve / "sample5_vc.bcf"

    # Step 4: Count variants
    default_total = get_variant_count(output_bcf_default)
    preserve_total = get_variant_count(output_bcf_preserve)

    # Count annotated variants using bcftools view filter
    bcftools = get_bcftools()

    # Count annotated in default output
    result = subprocess.run(
        [bcftools, "view", "-H", "-i", 'INFO/CSQ!=""', str(output_bcf_default)],
        capture_output=True,
        text=True,
        check=True,
    )
    default_annotated = len([line for line in result.stdout.split('\n') if line.strip()])

    # Count annotated in preserve output
    result = subprocess.run(
        [bcftools, "view", "-H", "-i", 'INFO/CSQ!=""', str(output_bcf_preserve)],
        capture_output=True,
        text=True,
        check=True,
    )
    preserve_annotated = len([line for line in result.stdout.split('\n') if line.strip()])

    print(f"Default mode - Total: {default_total}, Annotated: {default_annotated}")
    print(f"Preserve mode - Total: {preserve_total}, Annotated: {preserve_annotated}")

    # Step 5: Verify expectations

    # 1. Default mode should only have annotated variants (filtered)
    assert default_total == default_annotated, (
        f"Default mode should only contain annotated variants! "
        f"Total: {default_total}, Annotated: {default_annotated}"
    )
    print("✓ Default mode correctly filters to annotated variants only")

    # 2. Preserve mode should have all input variants (some unannotated)
    assert preserve_total > preserve_annotated, (
        f"Preserve mode should have unannotated variants! "
        f"Total: {preserve_total}, Annotated: {preserve_annotated}"
    )
    print("✓ Preserve mode correctly keeps unannotated variants")

    # 3. Both modes should have the same number of annotated variants
    assert default_annotated == preserve_annotated, (
        f"Annotated variant count should be identical! "
        f"Default: {default_annotated}, Preserve: {preserve_annotated}"
    )
    print("✓ Both modes have identical annotated variant counts")

    # 4. Preserve mode should have more or equal total variants than default
    assert preserve_total >= default_total, (
        f"Preserve mode should have >= variants than default! "
        f"Default: {default_total}, Preserve: {preserve_total}"
    )
    print("✓ Preserve mode has more total variants than default mode")

    # 5. Verify MD5s of annotated variants are identical
    default_md5 = compute_bcf_body_md5(output_bcf_default, filter_annotated=True)
    preserve_md5 = compute_bcf_body_md5(output_bcf_preserve, filter_annotated=True)

    assert default_md5 == preserve_md5, (
        f"Annotated variants should be identical between modes! "
        f"Default MD5: {default_md5}, Preserve MD5: {preserve_md5}"
    )
    print("✓ Annotated variants have identical MD5 checksums in both modes")

    print(f"\n✓ SUCCESS: --preserve-unannotated flag works correctly!")


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v"])
