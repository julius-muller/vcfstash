# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""Schema definitions and validation for vcfcache YAML configuration files.

This module provides formal schema definitions for:
- params.yaml: Environment and tool configuration
- annotation.yaml: Annotation workflow definition

Each schema defines required and optional fields with validation methods.
"""

from pathlib import Path
from typing import Any, Dict, Optional, Set


class ParamsYAMLSchema:
    """Schema definition for params.yaml configuration files.

    params.yaml files contain environment variables, tool paths, and resources
    used during vcfcache operations. These parameters can be referenced in
    annotation.yaml files using ${params.key_name} syntax.

    Required Fields:
        annotation_tool_cmd (str): Command to invoke the annotation tool
        bcftools_cmd (str): Command to invoke bcftools
        temp_dir (str): Directory for temporary files
        threads (int): Number of threads for bcftools operations

    Optional Fields:
        tool_version_command (str): Command to check annotation tool version
        optional_checks (dict): Additional validation parameters
        [any_custom_key]: Custom parameters for use in annotation.yaml templates

    Example:
        ```yaml
        # Required
        annotation_tool_cmd: "vep"
        bcftools_cmd: "bcftools"
        temp_dir: "/tmp"
        threads: 1
        genome_build: "GRCh38"

        # Optional
        tool_version_command: "vep --version"
        optional_checks: {}

        # Custom parameters
        vep_cache_dir: "/opt/vep/cache"
        vep_plugins_dir: "/opt/vep/plugins"
        ```
    """

    # Required fields
    REQUIRED_FIELDS: Set[str] = {
        "annotation_tool_cmd",
        "bcftools_cmd",
        "temp_dir",
        "threads",
        "genome_build",
    }

    # Optional fields (known/documented)
    OPTIONAL_FIELDS: Set[str] = {
        "tool_version_command",
        "optional_checks",
    }

    # Fields that should NOT appear (these belong in annotation.yaml)
    FORBIDDEN_FIELDS: Set[str] = {
        "annotation_cmd",
        "must_contain_info_tag",
        "required_tool_version",
    }

    @classmethod
    def validate(cls, data: Dict[str, Any], file_path: Optional[Path] = None) -> tuple[bool, Optional[str]]:
        """Validate a dictionary against the params.yaml schema.

        Args:
            data: Dictionary loaded from YAML file
            file_path: Optional path to the file being validated (for error messages)

        Returns:
            (is_valid, error_message): Tuple of validation result and error message if invalid
        """
        if not isinstance(data, dict):
            return False, f"params.yaml must contain a dictionary, got {type(data).__name__}"

        # Check for forbidden fields FIRST (indicates annotation.yaml was provided instead)
        # This provides a clearer error message than "missing required fields"
        forbidden_found = cls.FORBIDDEN_FIELDS.intersection(data.keys())
        if forbidden_found:
            file_hint = f" in {file_path}" if file_path else ""
            return False, (
                f"Invalid params.yaml{file_hint}: Contains annotation.yaml fields: {sorted(forbidden_found)}. "
                f"Did you accidentally provide annotation.yaml instead of params.yaml? "
                f"Use -a for annotation.yaml files (part of cache-build), and -y for params.yaml files."
            )

        # Check for required fields
        missing = cls.REQUIRED_FIELDS - data.keys()
        if missing:
            file_hint = f" in {file_path}" if file_path else ""
            # If missing fields AND has annotation.yaml-like structure, give helpful hint
            has_anno_markers = bool({'annotation_cmd', 'must_contain_info_tag', 'required_tool_version'}.intersection(data.keys()))
            hint = " (This looks like an annotation.yaml file - use -a, not -y)" if has_anno_markers else ""
            return False, (
                f"Invalid params.yaml{file_hint}: Missing required fields: {sorted(missing)}. "
                f"Required fields are: {sorted(cls.REQUIRED_FIELDS)}{hint}"
            )

        # Validate types of required fields
        for field_name in cls.REQUIRED_FIELDS:
            if field_name == "threads":
                # threads must be an integer
                if not isinstance(data[field_name], int):
                    return False, (
                        f"Invalid params.yaml: Field 'threads' must be an integer, "
                        f"got {type(data[field_name]).__name__}"
                    )
                if data[field_name] < 1:
                    return False, "Invalid params.yaml: Field 'threads' must be >= 1"
            else:
                # Other required fields must be strings
                if not isinstance(data[field_name], str):
                    return False, (
                        f"Invalid params.yaml: Field '{field_name}' must be a string, "
                        f"got {type(data[field_name]).__name__}"
                    )
                if field_name == "genome_build" and not data[field_name].strip():
                    return False, "Invalid params.yaml: Field 'genome_build' must be a non-empty string"

        # Validate optional_checks if present
        if "optional_checks" in data and not isinstance(data["optional_checks"], dict):
            return False, (
                f"Invalid params.yaml: Field 'optional_checks' must be a dictionary, "
                f"got {type(data['optional_checks']).__name__}"
            )

        return True, None

    @classmethod
    def get_documentation(cls) -> str:
        """Get human-readable documentation of the schema."""
        return f"""
params.yaml Schema
==================

Required Fields:
{chr(10).join(f'  - {field}' for field in sorted(cls.REQUIRED_FIELDS))}

Optional Fields:
{chr(10).join(f'  - {field}' for field in sorted(cls.OPTIONAL_FIELDS))}
  - [any custom key]: For use in annotation.yaml templates

Custom fields can be referenced in annotation.yaml using ${{params.key_name}} syntax.

Example:
--------
annotation_tool_cmd: "vep"
bcftools_cmd: "bcftools"
temp_dir: "/tmp"
genome_build: "GRCh38"
tool_version_command: "vep --version"
optional_checks: {{}}
vep_cache_dir: "/opt/vep/cache"
"""


class AnnotationYAMLSchema:
    """Schema definition for annotation.yaml configuration files.

    annotation.yaml files define the annotation workflow - the shell commands
    to execute for annotating variants. These files are stored with the cache
    and define how annotation should be performed.

    Required Fields:
        annotation_cmd (str): Shell commands to annotate INPUT_BCF → OUTPUT_BCF
        must_contain_info_tag (str): INFO tag that must appear in annotated output
        required_tool_version (str): Version marker for the annotation tool

    Optional Fields:
        optional_checks (dict): Additional validation parameters that must match
                               between cache and runtime params.yaml

    Available Template Variables:
        ${INPUT_BCF}: Path to input BCF file
        ${OUTPUT_BCF}: Path to output BCF file
        ${AUXILIARY_DIR}: Directory for auxiliary files
        ${params.*}: Any value from params.yaml (e.g., ${params.bcftools_cmd})

    Example:
        ```yaml
        annotation_cmd: |
          ${params.annotation_tool_cmd} \\
            --cache ${params.vep_cache_dir} \\
            --input_file ${INPUT_BCF} \\
            --output_file ${OUTPUT_BCF}

        must_contain_info_tag: CSQ
        required_tool_version: "115.2"
        genome_build: "GRCh38"
        optional_checks: {}
        ```
    """

    # Required fields
    REQUIRED_FIELDS: Set[str] = {
        "annotation_cmd",
        "must_contain_info_tag",
        "required_tool_version",
        "genome_build",
    }

    # Optional fields
    OPTIONAL_FIELDS: Set[str] = {
        "optional_checks",
    }

    # Fields that should NOT appear (these belong in params.yaml)
    FORBIDDEN_FIELDS: Set[str] = {
        "annotation_tool_cmd",
        "bcftools_cmd",
        "temp_dir",
        "tool_version_command",
    }

    @classmethod
    def validate(cls, data: Dict[str, Any], file_path: Optional[Path] = None) -> tuple[bool, Optional[str]]:
        """Validate a dictionary against the annotation.yaml schema.

        Args:
            data: Dictionary loaded from YAML file
            file_path: Optional path to the file being validated (for error messages)

        Returns:
            (is_valid, error_message): Tuple of validation result and error message if invalid
        """
        if not isinstance(data, dict):
            return False, f"annotation.yaml must contain a dictionary, got {type(data).__name__}"

        # Check for forbidden fields (indicates params.yaml was provided instead)
        forbidden_found = cls.FORBIDDEN_FIELDS.intersection(data.keys())
        if forbidden_found:
            file_hint = f" in {file_path}" if file_path else ""
            return False, (
                f"Invalid annotation.yaml{file_hint}: Contains params.yaml fields: {forbidden_found}. "
                f"Did you accidentally provide params.yaml instead of annotation.yaml?"
            )

        # Check for required fields
        missing = cls.REQUIRED_FIELDS - data.keys()
        if missing:
            file_hint = f" in {file_path}" if file_path else ""
            return False, (
                f"Invalid annotation.yaml{file_hint}: Missing required fields: {missing}. "
                f"Required fields are: {cls.REQUIRED_FIELDS}"
            )

        # Validate types of required fields
        for field_name in cls.REQUIRED_FIELDS:
            if not isinstance(data[field_name], str):
                return False, (
                    f"Invalid annotation.yaml: Field '{field_name}' must be a string, "
                    f"got {type(data[field_name]).__name__}"
                )
            if field_name == "genome_build" and not data[field_name].strip():
                return False, "Invalid annotation.yaml: Field 'genome_build' must be a non-empty string"

        # Validate optional_checks if present
        if "optional_checks" in data and not isinstance(data["optional_checks"], dict):
            return False, (
                f"Invalid annotation.yaml: Field 'optional_checks' must be a dictionary, "
                f"got {type(data['optional_checks']).__name__}"
            )

        return True, None

    @classmethod
    def get_documentation(cls) -> str:
        """Get human-readable documentation of the schema."""
        return f"""
annotation.yaml Schema
======================

Required Fields:
{chr(10).join(f'  - {field}' for field in sorted(cls.REQUIRED_FIELDS))}

Optional Fields:
{chr(10).join(f'  - {field}' for field in sorted(cls.OPTIONAL_FIELDS))}

Available Template Variables:
  - ${{INPUT_BCF}}: Path to input BCF file
  - ${{OUTPUT_BCF}}: Path to output BCF file
  - ${{AUXILIARY_DIR}}: Directory for auxiliary files
  - ${{params.*}}: Any value from params.yaml

Example:
--------
annotation_cmd: |
  ${{params.annotation_tool_cmd}} \\
    --input ${{INPUT_BCF}} \\
    --output ${{OUTPUT_BCF}}

must_contain_info_tag: CSQ
required_tool_version: "115.2"
genome_build: "GRCh38"
optional_checks: {{}}
"""


def validate_params_yaml(file_path: Path) -> tuple[bool, Optional[str], Optional[Dict[str, Any]]]:
    """Load and validate a params.yaml file.

    Args:
        file_path: Path to the params.yaml file

    Returns:
        (is_valid, error_message, data): Validation result, error message if invalid, and loaded data
    """
    import yaml

    try:
        with file_path.open() as f:
            data = yaml.safe_load(f)
    except Exception as e:
        return False, f"Failed to load YAML file {file_path}: {e}", None

    is_valid, error = ParamsYAMLSchema.validate(data, file_path)
    return is_valid, error, data


def validate_annotation_yaml(file_path: Path) -> tuple[bool, Optional[str], Optional[Dict[str, Any]]]:
    """Load and validate an annotation.yaml file.

    Args:
        file_path: Path to the annotation.yaml file

    Returns:
        (is_valid, error_message, data): Validation result, error message if invalid, and loaded data
    """
    import yaml

    try:
        with file_path.open() as f:
            data = yaml.safe_load(f)
    except Exception as e:
        return False, f"Failed to load YAML file {file_path}: {e}", None

    is_valid, error = AnnotationYAMLSchema.validate(data, file_path)
    return is_valid, error, data
