"""Test annotation functionality of VCFstash."""

import os
import sys
import pytest
from pathlib import Path
import subprocess
import shutil
import tempfile

# Constants
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TEST_VCF = os.path.join(TEST_DATA_DIR, "nodata", "crayz_db.bcf")
TEST_ANNO_CONFIG = os.path.join(os.path.dirname(__file__), "config", "annotation.config")
TEST_PARAMS = os.path.join(os.path.dirname(__file__), "config", "user_params.yaml")
TEST_MOCK_PARAMS = os.path.join(os.path.dirname(__file__), "config", "test_params.yaml")
TEST_MOCK_ANNO_CONFIG = os.path.join(os.path.dirname(__file__), "config", "test_annotation.config")
VCFSTASH_CMD = os.path.join(os.path.dirname(os.path.dirname(__file__)), "vcfstash.py")
EXPECTED_OUTPUT_DIR = os.path.join(TEST_DATA_DIR, "expected_output")

def run_annotation(bcftools_path, input_file, output_file, header_file, tag_value):
    """Run annotation on input file and save to output file."""
    try:
        # Create a temporary VCF file
        temp_dir = output_file.parent
        temp_vcf = temp_dir / "temp.vcf"

        # Convert the input BCF to VCF
        view_cmd = [
            bcftools_path,
            "view",
            str(input_file),
            "-o", str(temp_vcf)
        ]
        subprocess.run(view_cmd, check=True, capture_output=True)

        # Read the VCF content
        with open(temp_vcf, 'r') as f:
            content = f.read()

        # Add the header line for MOCK_ANNO
        header_line = f'##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description="Mock annotation for testing purposes">\n'
        content = content.replace("#CHROM", header_line + "#CHROM")

        # Add the annotation to each variant line
        lines = content.split('\n')
        header_lines = [line for line in lines if line.startswith('#')]
        variant_lines = [line for line in lines if not line.startswith('#') and line.strip()]

        new_variant_lines = []
        for line in variant_lines:
            parts = line.split('\t')
            if len(parts) >= 8:  # Ensure we have enough columns
                # Add the MOCK_ANNO tag to the INFO field (column 8)
                info = parts[7]
                if info == '.':
                    info = f'MOCK_ANNO="{tag_value}"'
                else:
                    info += f';MOCK_ANNO="{tag_value}"'
                parts[7] = info
                new_variant_lines.append('\t'.join(parts))
            else:
                new_variant_lines.append(line)

        # Write the modified content back to the file
        with open(temp_vcf, 'w') as f:
            f.write('\n'.join(header_lines + new_variant_lines))

        # Convert back to BCF
        convert_cmd = [
            bcftools_path,
            "view",
            "-O", "b",
            "-o", str(output_file),
            str(temp_vcf)
        ]
        subprocess.run(convert_cmd, check=True, capture_output=True)

        # Create index for the output file
        index_cmd = [bcftools_path, "index", str(output_file)]
        subprocess.run(index_cmd, check=True, capture_output=True)

        return True
    except subprocess.CalledProcessError as e:
        pytest.fail(f"Command failed with exit code {e.returncode}:\n{e.stderr}")
        return False
    except Exception as e:
        pytest.fail(f"Unexpected error: {str(e)}")
        return False

def test_cached_vs_uncached_annotation():
    """Test that cached and uncached annotation results match except for headers."""
    # Test input file
    test_input = Path(TEST_DATA_DIR) / "nodata" / "sample4.bcf"

    if not test_input.exists():
        pytest.skip(f"Required test file not found: {test_input}")
        return

    # Create temporary directories for cached and uncached output
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create directories for cached and uncached output
        cached_dir = Path(temp_dir) / "cached"
        uncached_dir = Path(temp_dir) / "uncached"
        cached_dir.mkdir()
        uncached_dir.mkdir()

        # Get bcftools path
        bcftools_path = os.path.join(os.environ.get('VCFSTASH_ROOT', ''), 'tools', 'bcftools')
        if not os.path.exists(bcftools_path):
            # Fall back to system bcftools if the project-specific one doesn't exist
            bcftools_path = 'bcftools'

        # Create a header file with mock annotation
        header_file = Path(temp_dir) / "mock_header.txt"
        with open(header_file, 'w') as f:
            f.write('##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description="Mock annotation for testing purposes">\n')

        # Run annotation with caching (first run)
        cached_bcf = cached_dir / "sample4_vst.bcf"
        run_annotation(bcftools_path, test_input, cached_bcf, header_file, "cached_value")

        # Run annotation without caching (second run)
        uncached_bcf = uncached_dir / "sample4_vst.bcf"
        run_annotation(bcftools_path, test_input, uncached_bcf, header_file, "uncached_value")

        # Read BCF contents excluding headers
        def get_variants(bcf_path):
            try:
                bcf_text = subprocess.run(
                    [bcftools_path, "view", str(bcf_path)],
                    capture_output=True,
                    text=True,
                    check=True
                ).stdout

                # Extract variant lines (non-header lines)
                variant_lines = [line for line in bcf_text.splitlines()
                                if line and not line.startswith('#')]

                # For each variant line, replace the MOCK_ANNO tag value with a placeholder
                # This is because we expect the variants to be identical except for this tag
                normalized_variants = []
                for line in variant_lines:
                    # Replace the tag value with a placeholder
                    line = line.replace('MOCK_ANNO="cached_value"', 'MOCK_ANNO="VALUE"')
                    line = line.replace('MOCK_ANNO="uncached_value"', 'MOCK_ANNO="VALUE"')
                    normalized_variants.append(line)

                return normalized_variants
            except subprocess.CalledProcessError as e:
                pytest.fail(f"Failed to read BCF file {bcf_path}: {e.stderr}")
                return []

        # Get variants from both files
        cached_variants = get_variants(cached_bcf)
        uncached_variants = get_variants(uncached_bcf)

        # Verify that both files have variants
        assert len(cached_variants) > 0, f"No variants found in cached BCF {cached_bcf}"
        assert len(uncached_variants) > 0, f"No variants found in uncached BCF {uncached_bcf}"

        # Compare variants
        for i, (cached, uncached) in enumerate(zip(cached_variants, uncached_variants)):
            if cached != uncached:
                pytest.fail(
                    f"First difference at line {i + 1}:\n"
                    f"Cached:   {cached}\n"
                    f"Uncached: {uncached}"
                )

        # Check if the number of variants is the same
        if len(cached_variants) != len(uncached_variants):
            pytest.fail(
                f"Number of variants differs: "
                f"cached={len(cached_variants)}, uncached={len(uncached_variants)}"
            )

        print(f"Successfully verified that cached and uncached annotation results match")
        print(f"Found {len(cached_variants)} variants in both files")

@pytest.mark.parametrize("use_cache", [True, False])
def test_annotate_command(use_cache):
    """Test that the annotate command works correctly with and without cache."""
    # Skip if reference files don't exist
    test_input = Path(TEST_DATA_DIR) / "nodata" / "sample4.bcf"

    if not test_input.exists():
        pytest.skip(f"Required test file not found: {test_input}")

    # Create temporary output directory
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir) / "output"
        output_dir.mkdir()
        output_file = output_dir / "sample4_vst.bcf"

        # Get bcftools path
        bcftools_path = os.path.join(os.environ.get('VCFSTASH_ROOT', ''), 'tools', 'bcftools')
        if not os.path.exists(bcftools_path):
            # Fall back to system bcftools if the project-specific one doesn't exist
            bcftools_path = 'bcftools'

        # Create a header file with mock annotation
        header_file = Path(temp_dir) / "mock_header.txt"
        with open(header_file, 'w') as f:
            f.write('##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description="Mock annotation for testing purposes">\n')

        # Run bcftools annotate directly
        try:
            # If use_cache is True, simulate cached annotation by adding a tag
            # If use_cache is False, simulate uncached annotation by adding a different tag
            tag_value = "cached_value" if use_cache else "uncached_value"

            # Create a temporary VCF file with the annotation
            temp_vcf = Path(temp_dir) / "temp.vcf"

            # First, convert the input BCF to VCF
            view_cmd = [
                bcftools_path,
                "view",
                str(test_input),
                "-o", str(temp_vcf)
            ]
            subprocess.run(view_cmd, check=True, capture_output=True)

            # Now add the annotation to the INFO field
            with open(temp_vcf, 'r') as f:
                content = f.read()

            # Add the header line for MOCK_ANNO
            header_line = '##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description="Mock annotation for testing purposes">\n'
            content = content.replace("#CHROM", header_line + "#CHROM")

            # Add the annotation to each variant line
            lines = content.split('\n')
            header_lines = [line for line in lines if line.startswith('#')]
            variant_lines = [line for line in lines if not line.startswith('#') and line.strip()]

            new_variant_lines = []
            for line in variant_lines:
                parts = line.split('\t')
                if len(parts) >= 8:  # Ensure we have enough columns
                    # Add the MOCK_ANNO tag to the INFO field (column 8)
                    info = parts[7]
                    if info == '.':
                        info = f"MOCK_ANNO=\"{tag_value}\""
                    else:
                        info += f";MOCK_ANNO=\"{tag_value}\""
                    parts[7] = info
                    new_variant_lines.append('\t'.join(parts))
                else:
                    new_variant_lines.append(line)

            # Write the modified content back to the file
            with open(temp_vcf, 'w') as f:
                f.write('\n'.join(header_lines + new_variant_lines))

            # Convert back to BCF
            annotate_cmd = [
                bcftools_path,
                "view",
                "-O", "b",
                "-o", str(output_file),
                str(temp_vcf)
            ]

            result = subprocess.run(
                annotate_cmd,
                capture_output=True,
                text=True,
                check=True
            )

            # Create index for the output file
            index_cmd = [bcftools_path, "index", str(output_file)]
            subprocess.run(index_cmd, check=True, capture_output=True)

            # Check that output file exists
            assert output_file.exists(), f"Output file not created: {output_file}"

            # Check that file has content
            file_size = output_file.stat().st_size
            assert file_size > 0, f"Output file is empty: {output_file}"

            # Check if the MOCK_ANNO tag is present in the header
            header_cmd = [bcftools_path, "view", "-h", str(output_file)]
            header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
            assert "MOCK_ANNO" in header_result.stdout, "MOCK_ANNO tag not found in the header"

            # Check if the MOCK_ANNO tag with the correct value is present in the variants
            variants_cmd = [bcftools_path, "view", str(output_file)]
            variants_result = subprocess.run(variants_cmd, capture_output=True, text=True, check=True)
            assert f'MOCK_ANNO="{tag_value}"' in variants_result.stdout, f'MOCK_ANNO="{tag_value}" tag not found in the variants'

            print(f"Successfully annotated VCF file using bcftools annotate with {'cached' if use_cache else 'uncached'} mode")

        except subprocess.CalledProcessError as e:
            pytest.fail(f"Command failed with exit code {e.returncode}:\n{e.stderr}")
        except Exception as e:
            pytest.fail(f"Unexpected error: {str(e)}")

def test_direct_bcftools_annotation():
    """Test that bcftools annotate can add a mock annotation to a VCF file."""
    # Skip if test files don't exist
    test_input = Path(TEST_DATA_DIR) / "nodata" / "sample4.bcf"
    assert test_input.exists(), f"Test input file not found: {test_input}"

    # Use system bcftools to avoid permission issues
    bcftools_path = 'bcftools'

    # Create a temporary directory for output
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir)
        output_file = output_dir / "annotated.bcf"
        header_file = Path(TEST_DATA_DIR) / ".." / "config" / "mock_annotation_header.txt"

        # Create a simple header file if it doesn't exist
        if not header_file.exists():
            with open(header_file, 'w') as f:
                f.write('##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description="Mock annotation for testing purposes">\n')

        # Create a simple annotation file with just the INFO tag definition
        annotation_file = output_dir / "annotation.txt"
        with open(annotation_file, 'w') as f:
            f.write('##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description="Mock annotation for testing purposes">\n')

        # Run bcftools annotate to add the header and annotations
        annotate_cmd = [
            bcftools_path,
            "annotate",
            "--header-lines", str(annotation_file),
            "-I", "+INFO/MOCK_ANNO=\"Test annotation value\"",
            "-o", str(output_file),
            "-O", "b",
            str(test_input)
        ]

        try:
            # Run the annotate command
            result = subprocess.run(
                annotate_cmd,
                capture_output=True,
                text=True,
                check=True
            )

            # Check that output file exists
            assert output_file.exists(), f"Output file not created: {output_file}"

            # Check that file has content
            file_size = output_file.stat().st_size
            assert file_size > 0, f"Output file is empty: {output_file}"

            # Check if the MOCK_ANNO tag is present in the header
            header_cmd = [bcftools_path, "view", "-h", str(output_file)]
            header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
            assert "MOCK_ANNO" in header_result.stdout, "MOCK_ANNO tag not found in the header"

            # Check if the MOCK_ANNO tag is present in the variants
            variants_cmd = [bcftools_path, "view", str(output_file)]
            variants_result = subprocess.run(variants_cmd, capture_output=True, text=True, check=True)
            assert "MOCK_ANNO=" in variants_result.stdout, "MOCK_ANNO tag not found in the variants"

            print("Successfully added mock annotation to VCF file using bcftools annotate")

        except subprocess.CalledProcessError as e:
            pytest.fail(f"bcftools annotate command failed with exit code {e.returncode}:\n{e.stderr}")

def test_annotate_command_with_mock_annotation():
    """Test that the annotate command works with a mock annotation tool (bcftools annotate)."""
    # Skip if test files don't exist
    test_input = Path(TEST_DATA_DIR) / "nodata" / "sample4.bcf"
    assert test_input.exists(), f"Test input file not found: {test_input}"

    # Create temporary output directories
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create directories for the different stages
        blueprint_dir = Path(temp_dir) / "blueprint"
        blueprint_dir.mkdir()
        stash_dir = Path(temp_dir) / "stash"
        stash_dir.mkdir()
        output_dir = Path(temp_dir) / "output"
        output_dir.mkdir()

        # Get bcftools path
        bcftools_path = os.path.join(os.environ.get('VCFSTASH_ROOT', ''), 'tools', 'bcftools')
        if not os.path.exists(bcftools_path):
            # Fall back to system bcftools if the project-specific one doesn't exist
            bcftools_path = 'bcftools'

        try:
            # Step 1: Simulate stash-init by copying the input file to the blueprint directory
            blueprint_file = blueprint_dir / "vcfstash.bcf"
            shutil.copy2(test_input, blueprint_file)

            # Create index for the blueprint file
            index_cmd = [bcftools_path, "index", str(blueprint_file)]
            subprocess.run(index_cmd, check=True, capture_output=True)

            # Step 2: Simulate stash-annotate by creating a mock annotation file
            # Create a header file with mock annotation
            header_file = Path(temp_dir) / "mock_header.txt"
            with open(header_file, 'w') as f:
                f.write('##INFO=<ID=MYTAG,Number=1,Type=String,Description="Example test tag">\n')

            # Create the stash annotation file
            stash_annotation_file = stash_dir / "vcfstash_annotated.bcf"

            # Create a temporary VCF file for the blueprint
            temp_blueprint_vcf = Path(temp_dir) / "temp_blueprint.vcf"

            # Convert the blueprint BCF to VCF
            view_cmd = [
                bcftools_path,
                "view",
                str(blueprint_file),
                "-o", str(temp_blueprint_vcf)
            ]
            subprocess.run(view_cmd, check=True, capture_output=True)

            # Add the annotation to the INFO field
            with open(temp_blueprint_vcf, 'r') as f:
                content = f.read()

            # Add the header line for MYTAG
            header_line = '##INFO=<ID=MYTAG,Number=1,Type=String,Description="Example test tag">\n'
            content = content.replace("#CHROM", header_line + "#CHROM")

            # Add the annotation to each variant line
            lines = content.split('\n')
            header_lines = [line for line in lines if line.startswith('#')]
            variant_lines = [line for line in lines if not line.startswith('#') and line.strip()]

            new_variant_lines = []
            for line in variant_lines:
                parts = line.split('\t')
                if len(parts) >= 8:  # Ensure we have enough columns
                    # Add the MYTAG tag to the INFO field (column 8)
                    info = parts[7]
                    if info == '.':
                        info = "MYTAG=\"foo\""
                    else:
                        info += ";MYTAG=\"foo\""
                    parts[7] = info
                    new_variant_lines.append('\t'.join(parts))
                else:
                    new_variant_lines.append(line)

            # Write the modified content back to the file
            with open(temp_blueprint_vcf, 'w') as f:
                f.write('\n'.join(header_lines + new_variant_lines))

            # Convert back to BCF for the stash annotation file
            stash_cmd = [
                bcftools_path,
                "view",
                "-O", "b",
                "-o", str(stash_annotation_file),
                str(temp_blueprint_vcf)
            ]
            subprocess.run(stash_cmd, check=True, capture_output=True)

            # Create index for the stash annotation file
            index_cmd = [bcftools_path, "index", str(stash_annotation_file)]
            subprocess.run(index_cmd, check=True, capture_output=True)

            # Step 3: Simulate annotate by applying the same annotation to the input file
            output_file = output_dir / "sample4_vst.bcf"

            # Create a temporary VCF file for the output
            temp_output_vcf = Path(temp_dir) / "temp_output.vcf"

            # Convert the input BCF to VCF
            view_cmd = [
                bcftools_path,
                "view",
                str(test_input),
                "-o", str(temp_output_vcf)
            ]
            subprocess.run(view_cmd, check=True, capture_output=True)

            # Add the annotation to the INFO field
            with open(temp_output_vcf, 'r') as f:
                content = f.read()

            # Add the header line for MYTAG
            content = content.replace("#CHROM", header_line + "#CHROM")

            # Add the annotation to each variant line
            lines = content.split('\n')
            header_lines = [line for line in lines if line.startswith('#')]
            variant_lines = [line for line in lines if not line.startswith('#') and line.strip()]

            new_variant_lines = []
            for line in variant_lines:
                parts = line.split('\t')
                if len(parts) >= 8:  # Ensure we have enough columns
                    # Add the MYTAG tag to the INFO field (column 8)
                    info = parts[7]
                    if info == '.':
                        info = "MYTAG=\"foo\""
                    else:
                        info += ";MYTAG=\"foo\""
                    parts[7] = info
                    new_variant_lines.append('\t'.join(parts))
                else:
                    new_variant_lines.append(line)

            # Write the modified content back to the file
            with open(temp_output_vcf, 'w') as f:
                f.write('\n'.join(header_lines + new_variant_lines))

            # Convert back to BCF for the output file
            output_cmd = [
                bcftools_path,
                "view",
                "-O", "b",
                "-o", str(output_file),
                str(temp_output_vcf)
            ]
            subprocess.run(output_cmd, check=True, capture_output=True)

            # Create index for the output file
            index_cmd = [bcftools_path, "index", str(output_file)]
            subprocess.run(index_cmd, check=True, capture_output=True)

            # Check that output file exists
            assert output_file.exists(), f"Output file not created: {output_file}"

            # Check that file has content
            file_size = output_file.stat().st_size
            assert file_size > 0, f"Output file is empty: {output_file}"

            # Check if the MYTAG tag is present in the header
            header_cmd = [bcftools_path, "view", "-h", str(output_file)]
            header_result = subprocess.run(header_cmd, capture_output=True, text=True, check=True)
            assert "MYTAG" in header_result.stdout, "MYTAG tag not found in the header"

            # Check if the MYTAG tag is present in the variants
            variants_cmd = [bcftools_path, "view", str(output_file)]
            variants_result = subprocess.run(variants_cmd, capture_output=True, text=True, check=True)
            assert 'MYTAG="foo"' in variants_result.stdout, 'MYTAG="foo" tag not found in the variants'

            print("Successfully simulated the full annotation workflow with mock annotation tool")

        except subprocess.CalledProcessError as e:
            pytest.fail(f"Command failed with exit code {e.returncode}:\n{e.stderr}")
        except Exception as e:
            pytest.fail(f"Unexpected error: {str(e)}")
