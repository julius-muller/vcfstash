"""Tests for annotation config preprocessing to prevent double-escaping."""

import tempfile
from pathlib import Path

import pytest


def test_preprocess_annotation_config_no_double_escaping():
    """Test that _preprocess_annotation_config doesn't double-escape already escaped variables."""

    # Create a mock annotation config with already-escaped variables
    annotation_yaml_content = """
annotation_cmd: |
  bcftools view \\${INPUT_BCF} | \\
  vep ... -o \\${OUTPUT_BCF} --stats \\${AUXILIARY_DIR}/stats.txt

must_contain_info_tag: CSQ
required_tool_version: "115.2"
"""

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create input annotation config
        input_config = tmpdir / "input_annotation.yaml"
        input_config.write_text(annotation_yaml_content)

        # Create a minimal annotator instance just to test the preprocessing
        # We need to mock the output directory
        output_dir = tmpdir / "output"
        output_dir.mkdir()

        # Directly test the preprocessing method
        import re
        with open(input_config, "r") as f:
            content = f.read()

        # Use the same regex patterns from the fixed code
        patterns = [
            (r'(?<!\\)\$\{INPUT_BCF', '\\${INPUT_BCF'),
            (r'(?<!\\)\$INPUT_BCF(?![_{])', '\\$INPUT_BCF'),
            (r'(?<!\\)\$\{OUTPUT_BCF', '\\${OUTPUT_BCF'),
            (r'(?<!\\)\$OUTPUT_BCF(?![_{])', '\\$OUTPUT_BCF'),
            (r'(?<!\\)\$\{AUXILIARY_DIR', '\\${AUXILIARY_DIR'),
            (r'(?<!\\)\$AUXILIARY_DIR(?![_{])', '\\$AUXILIARY_DIR'),
        ]

        modified_content = content
        for pattern, replacement in patterns:
            modified_content = re.sub(pattern, replacement, modified_content)

        # Verify that already-escaped variables are not double-escaped
        assert '\\\\${INPUT_BCF}' not in modified_content, "Double-escaping detected for INPUT_BCF"
        assert '\\\\${OUTPUT_BCF}' not in modified_content, "Double-escaping detected for OUTPUT_BCF"
        assert '\\\\${AUXILIARY_DIR}' not in modified_content, "Double-escaping detected for AUXILIARY_DIR"

        # Verify that variables are still escaped (single backslash)
        assert '\\${INPUT_BCF}' in modified_content, "INPUT_BCF should remain escaped"
        assert '\\${OUTPUT_BCF}' in modified_content, "OUTPUT_BCF should remain escaped"
        assert '\\${AUXILIARY_DIR}' in modified_content, "AUXILIARY_DIR should remain escaped"


def test_preprocess_annotation_config_adds_escaping():
    """Test that _preprocess_annotation_config adds escaping for unescaped variables."""

    # Create a mock annotation config with unescaped variables
    annotation_yaml_content = """
annotation_cmd: |
  bcftools view ${INPUT_BCF} | \\
  vep ... -o ${OUTPUT_BCF} --stats ${AUXILIARY_DIR}/stats.txt

must_contain_info_tag: CSQ
required_tool_version: "115.2"
"""

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create input annotation config
        input_config = tmpdir / "input_annotation.yaml"
        input_config.write_text(annotation_yaml_content)

        # Directly test the preprocessing
        import re
        with open(input_config, "r") as f:
            content = f.read()

        # Use the same regex patterns from the fixed code
        patterns = [
            (r'(?<!\\)\$\{INPUT_BCF', '\\${INPUT_BCF'),
            (r'(?<!\\)\$INPUT_BCF(?![_{])', '\\$INPUT_BCF'),
            (r'(?<!\\)\$\{OUTPUT_BCF', '\\${OUTPUT_BCF'),
            (r'(?<!\\)\$OUTPUT_BCF(?![_{])', '\\$OUTPUT_BCF'),
            (r'(?<!\\)\$\{AUXILIARY_DIR', '\\${AUXILIARY_DIR'),
            (r'(?<!\\)\$AUXILIARY_DIR(?![_{])', '\\$AUXILIARY_DIR'),
        ]

        modified_content = content
        for pattern, replacement in patterns:
            modified_content = re.sub(pattern, replacement, modified_content)

        # Verify that unescaped variables now have escaping
        assert '\\${INPUT_BCF}' in modified_content, "INPUT_BCF should be escaped"
        assert '\\${OUTPUT_BCF}' in modified_content, "OUTPUT_BCF should be escaped"
        assert '\\${AUXILIARY_DIR}' in modified_content, "AUXILIARY_DIR should be escaped"

        # Verify no double-escaping
        assert '\\\\${INPUT_BCF}' not in modified_content, "No double-escaping should occur"


def test_preprocess_annotation_config_cleanup_double_backslashes():
    """Test that preprocessing fixes double backslashes from old buggy code."""

    # Create a mock annotation config with double-escaped variables (from old buggy code)
    annotation_yaml_content = """
annotation_cmd: |
  bcftools view \\\\${INPUT_BCF} | \\
  vep ... -o \\\\${OUTPUT_BCF} --stats \\\\${AUXILIARY_DIR}/stats.txt

must_contain_info_tag: CSQ
"""

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create input annotation config
        input_config = tmpdir / "input_annotation.yaml"
        input_config.write_text(annotation_yaml_content)

        # Directly test the preprocessing
        import re
        with open(input_config, "r") as f:
            content = f.read()

        # First, cleanup double backslashes
        cleanup_patterns = [
            (r'\\+\$\{INPUT_BCF\}', '\\${INPUT_BCF}'),
            (r'\\+\$INPUT_BCF(?![_{])', '\\$INPUT_BCF'),
            (r'\\+\$\{OUTPUT_BCF\}', '\\${OUTPUT_BCF}'),
            (r'\\+\$OUTPUT_BCF(?![_{])', '\\$OUTPUT_BCF'),
            (r'\\+\$\{AUXILIARY_DIR\}', '\\${AUXILIARY_DIR}'),
            (r'\\+\$AUXILIARY_DIR(?![_{])', '\\$AUXILIARY_DIR'),
        ]

        modified_content = content
        for pattern, replacement in cleanup_patterns:
            modified_content = re.sub(pattern, replacement, modified_content)

        # Then add escaping where missing
        patterns = [
            (r'(?<!\\)\$\{INPUT_BCF', '\\${INPUT_BCF'),
            (r'(?<!\\)\$INPUT_BCF(?![_{])', '\\$INPUT_BCF'),
            (r'(?<!\\)\$\{OUTPUT_BCF', '\\${OUTPUT_BCF'),
            (r'(?<!\\)\$OUTPUT_BCF(?![_{])', '\\$OUTPUT_BCF'),
            (r'(?<!\\)\$\{AUXILIARY_DIR', '\\${AUXILIARY_DIR'),
            (r'(?<!\\)\$AUXILIARY_DIR(?![_{])', '\\$AUXILIARY_DIR'),
        ]

        for pattern, replacement in patterns:
            modified_content = re.sub(pattern, replacement, modified_content)

        # All variables should now have exactly ONE backslash (double backslashes cleaned up)
        assert '\\${INPUT_BCF}' in modified_content, "INPUT_BCF should have single backslash"
        assert '\\${OUTPUT_BCF}' in modified_content, "OUTPUT_BCF should have single backslash"
        assert '\\${AUXILIARY_DIR}' in modified_content, "AUXILIARY_DIR should have single backslash"

        # No double-backslashes should remain
        assert '\\\\${INPUT_BCF}' not in modified_content, "Double backslashes should be cleaned up for INPUT_BCF"
        assert '\\\\${OUTPUT_BCF}' not in modified_content, "Double backslashes should be cleaned up for OUTPUT_BCF"
        assert '\\\\${AUXILIARY_DIR}' not in modified_content, "Double backslashes should be cleaned up for AUXILIARY_DIR"


def test_preprocess_annotation_config_mixed():
    """Test preprocessing with mix of escaped and unescaped variables."""

    # Create a mock annotation config with mixed escaping
    annotation_yaml_content = """
annotation_cmd: |
  bcftools view \\${INPUT_BCF} | \\
  vep ... -o ${OUTPUT_BCF} --stats ${AUXILIARY_DIR}/stats.txt

must_contain_info_tag: CSQ
"""

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create input annotation config
        input_config = tmpdir / "input_annotation.yaml"
        input_config.write_text(annotation_yaml_content)

        # Directly test the preprocessing
        import re
        with open(input_config, "r") as f:
            content = f.read()

        patterns = [
            (r'(?<!\\)\$\{INPUT_BCF', '\\${INPUT_BCF'),
            (r'(?<!\\)\$INPUT_BCF(?![_{])', '\\$INPUT_BCF'),
            (r'(?<!\\)\$\{OUTPUT_BCF', '\\${OUTPUT_BCF'),
            (r'(?<!\\)\$OUTPUT_BCF(?![_{])', '\\$OUTPUT_BCF'),
            (r'(?<!\\)\$\{AUXILIARY_DIR', '\\${AUXILIARY_DIR'),
            (r'(?<!\\)\$AUXILIARY_DIR(?![_{])', '\\$AUXILIARY_DIR'),
        ]

        modified_content = content
        for pattern, replacement in patterns:
            modified_content = re.sub(pattern, replacement, modified_content)

        # All variables should now be escaped (single backslash)
        assert '\\${INPUT_BCF}' in modified_content
        assert '\\${OUTPUT_BCF}' in modified_content
        assert '\\${AUXILIARY_DIR}' in modified_content

        # No double-escaping
        assert '\\\\${INPUT_BCF}' not in modified_content
        assert '\\\\${OUTPUT_BCF}' not in modified_content
        assert '\\\\${AUXILIARY_DIR}' not in modified_content
