"""Tests for YAML schema validation."""

import pytest
from pathlib import Path
from vcfcache.utils.schemas import (
    ParamsYAMLSchema,
    AnnotationYAMLSchema,
    validate_params_yaml,
    validate_annotation_yaml,
)


class TestParamsYAMLSchema:
    """Test params.yaml schema validation."""

    def test_valid_params_yaml(self):
        """Valid params.yaml should pass validation."""
        data = {
            "annotation_tool_cmd": "vep",
            "bcftools_cmd": "bcftools",
            "temp_dir": "/tmp",
            "threads": 2,
            "genome_build": "GRCh38",
            "tool_version_command": "vep --version",
            "optional_checks": {},
            "custom_field": "value",
        }
        is_valid, error = ParamsYAMLSchema.validate(data)
        assert is_valid
        assert error is None

    def test_minimal_valid_params_yaml(self):
        """Minimal params.yaml with only required fields should pass."""
        data = {
            "annotation_tool_cmd": "vep",
            "bcftools_cmd": "bcftools",
            "temp_dir": "/tmp",
            "threads": 1,
            "genome_build": "GRCh38",
        }
        is_valid, error = ParamsYAMLSchema.validate(data)
        assert is_valid
        assert error is None

    def test_params_yaml_missing_required_field(self):
        """params.yaml missing required field should fail."""
        data = {
            "annotation_tool_cmd": "vep",
            "bcftools_cmd": "bcftools",
            # Missing temp_dir
        }
        is_valid, error = ParamsYAMLSchema.validate(data)
        assert not is_valid
        assert "Missing required fields" in error
        assert "temp_dir" in error

    def test_params_yaml_with_annotation_fields(self):
        """params.yaml with annotation.yaml fields should fail with helpful message."""
        data = {
            "annotation_cmd": "echo test",
            "must_contain_info_tag": "CSQ",
            "required_tool_version": "115.2",
            "bcftools_cmd": "bcftools",
        }
        is_valid, error = ParamsYAMLSchema.validate(data)
        assert not is_valid
        assert "annotation.yaml fields" in error
        assert "annotation_cmd" in error
        assert "Did you accidentally provide annotation.yaml instead" in error

    def test_params_yaml_invalid_type(self):
        """params.yaml with invalid field type should fail."""
        data = {
            "annotation_tool_cmd": 123,  # Should be string
            "bcftools_cmd": "bcftools",
            "temp_dir": "/tmp",
            "threads": 1,
            "genome_build": "GRCh38",
        }
        is_valid, error = ParamsYAMLSchema.validate(data)
        assert not is_valid
        assert "must be a string" in error

    def test_params_yaml_not_dict(self):
        """params.yaml that's not a dictionary should fail."""
        data = "not a dict"
        is_valid, error = ParamsYAMLSchema.validate(data)
        assert not is_valid
        assert "must contain a dictionary" in error

    def test_params_yaml_invalid_optional_checks(self):
        """params.yaml with non-dict optional_checks should fail."""
        data = {
            "annotation_tool_cmd": "vep",
            "bcftools_cmd": "bcftools",
            "temp_dir": "/tmp",
            "threads": 1,
            "genome_build": "GRCh38",
            "optional_checks": "not a dict",
        }
        is_valid, error = ParamsYAMLSchema.validate(data)
        assert not is_valid
        assert "optional_checks" in error
        assert "must be a dictionary" in error

    def test_params_yaml_threads_not_int(self):
        """params.yaml with non-integer threads should fail."""
        data = {
            "annotation_tool_cmd": "vep",
            "bcftools_cmd": "bcftools",
            "temp_dir": "/tmp",
            "threads": "not an int",
            "genome_build": "GRCh38",
        }
        is_valid, error = ParamsYAMLSchema.validate(data)
        assert not is_valid
        assert "threads" in error
        assert "must be an integer" in error

    def test_params_yaml_threads_less_than_one(self):
        """params.yaml with threads < 1 should fail."""
        data = {
            "annotation_tool_cmd": "vep",
            "bcftools_cmd": "bcftools",
            "temp_dir": "/tmp",
            "threads": 0,
            "genome_build": "GRCh38",
        }
        is_valid, error = ParamsYAMLSchema.validate(data)
        assert not is_valid
        assert "threads" in error
        assert ">= 1" in error


class TestAnnotationYAMLSchema:
    """Test annotation.yaml schema validation."""

    def test_valid_annotation_yaml(self):
        """Valid annotation.yaml should pass validation."""
        data = {
            "annotation_cmd": "echo test",
            "must_contain_info_tag": "CSQ",
            "required_tool_version": "115.2",
            "genome_build": "GRCh38",
            "optional_checks": {},
        }
        is_valid, error = AnnotationYAMLSchema.validate(data)
        assert is_valid
        assert error is None

    def test_minimal_valid_annotation_yaml(self):
        """Minimal annotation.yaml with only required fields should pass."""
        data = {
            "annotation_cmd": "echo test",
            "must_contain_info_tag": "CSQ",
            "required_tool_version": "115.2",
            "genome_build": "GRCh38",
        }
        is_valid, error = AnnotationYAMLSchema.validate(data)
        assert is_valid
        assert error is None

    def test_annotation_yaml_missing_required_field(self):
        """annotation.yaml missing required field should fail."""
        data = {
            "annotation_cmd": "echo test",
            "must_contain_info_tag": "CSQ",
            # Missing required_tool_version
        }
        is_valid, error = AnnotationYAMLSchema.validate(data)
        assert not is_valid
        assert "Missing required fields" in error
        assert "required_tool_version" in error

    def test_annotation_yaml_with_params_fields(self):
        """annotation.yaml with params.yaml fields should fail with helpful message."""
        data = {
            "annotation_tool_cmd": "vep",
            "bcftools_cmd": "bcftools",
            "temp_dir": "/tmp",
            "annotation_cmd": "echo test",
        }
        is_valid, error = AnnotationYAMLSchema.validate(data)
        assert not is_valid
        assert "params.yaml fields" in error
        assert "annotation_tool_cmd" in error or "bcftools_cmd" in error
        assert "Did you accidentally provide params.yaml instead" in error

    def test_annotation_yaml_invalid_type(self):
        """annotation.yaml with invalid field type should fail."""
        data = {
            "annotation_cmd": 123,  # Should be string
            "must_contain_info_tag": "CSQ",
            "required_tool_version": "115.2",
            "genome_build": "GRCh38",
        }
        is_valid, error = AnnotationYAMLSchema.validate(data)
        assert not is_valid
        assert "must be a string" in error

    def test_annotation_yaml_not_dict(self):
        """annotation.yaml that's not a dictionary should fail."""
        data = ["not", "a", "dict"]
        is_valid, error = AnnotationYAMLSchema.validate(data)
        assert not is_valid
        assert "must contain a dictionary" in error


class TestFileValidation:
    """Test file-based validation functions."""

    def test_validate_params_yaml_file(self, tmp_path):
        """Test validating a params.yaml file."""
        params_file = tmp_path / "params.yaml"
        params_file.write_text("""
annotation_tool_cmd: vep
bcftools_cmd: bcftools
temp_dir: /tmp
threads: 1
genome_build: GRCh38
""")
        is_valid, error, data = validate_params_yaml(params_file)
        assert is_valid
        assert error is None
        assert data["annotation_tool_cmd"] == "vep"

    def test_validate_annotation_yaml_file(self, tmp_path):
        """Test validating an annotation.yaml file."""
        anno_file = tmp_path / "annotation.yaml"
        anno_file.write_text("""
annotation_cmd: echo test
must_contain_info_tag: CSQ
required_tool_version: "115.2"
genome_build: GRCh38
""")
        is_valid, error, data = validate_annotation_yaml(anno_file)
        assert is_valid
        assert error is None
        assert data["must_contain_info_tag"] == "CSQ"

    def test_validate_nonexistent_file(self, tmp_path):
        """Test validating a non-existent file."""
        fake_file = tmp_path / "nonexistent.yaml"
        is_valid, error, data = validate_params_yaml(fake_file)
        assert not is_valid
        assert "Failed to load" in error
        assert data is None

    def test_validate_invalid_yaml(self, tmp_path):
        """Test validating malformed YAML."""
        bad_file = tmp_path / "bad.yaml"
        bad_file.write_text("not: valid: yaml: content:")
        is_valid, error, data = validate_params_yaml(bad_file)
        assert not is_valid
        assert "Failed to load" in error
        assert data is None


class TestSchemaDocumentation:
    """Test schema documentation methods."""

    def test_params_schema_documentation(self):
        """Test that params schema has documentation."""
        doc = ParamsYAMLSchema.get_documentation()
        assert "params.yaml Schema" in doc
        assert "annotation_tool_cmd" in doc
        assert "bcftools_cmd" in doc
        assert "temp_dir" in doc

    def test_annotation_schema_documentation(self):
        """Test that annotation schema has documentation."""
        doc = AnnotationYAMLSchema.get_documentation()
        assert "annotation.yaml Schema" in doc
        assert "annotation_cmd" in doc
        assert "must_contain_info_tag" in doc
        assert "required_tool_version" in doc


class TestCrossValidation:
    """Test that schemas correctly reject each other's files."""

    def test_reject_annotation_as_params(self, tmp_path):
        """Test that annotation.yaml is rejected when params.yaml is expected."""
        anno_file = tmp_path / "annotation.yaml"
        anno_file.write_text("""
annotation_cmd: echo test
must_contain_info_tag: CSQ
required_tool_version: "115.2"
genome_build: GRCh38
""")
        is_valid, error, data = validate_params_yaml(anno_file)
        assert not is_valid
        assert "annotation.yaml fields" in error or "Missing required fields" in error

    def test_reject_params_as_annotation(self, tmp_path):
        """Test that params.yaml is rejected when annotation.yaml is expected."""
        params_file = tmp_path / "params.yaml"
        params_file.write_text("""
annotation_tool_cmd: vep
bcftools_cmd: bcftools
temp_dir: /tmp
threads: 1
genome_build: GRCh38
""")
        is_valid, error, data = validate_annotation_yaml(params_file)
        assert not is_valid
        assert "params.yaml fields" in error or "Missing required fields" in error
