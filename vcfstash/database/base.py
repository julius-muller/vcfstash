import importlib.resources
import json
import os
import re
import shutil
import subprocess
import tempfile
from logging import Logger
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pysam
import yaml

from vcfstash.database.outputs import StashOutput
from vcfstash.utils.logging import setup_logging
from vcfstash.utils.validation import validate_bcf_header


class NextflowWorkflow:
    """Base class for Nextflow workflow operations for a single input file.

    The NextflowWorkflow class is designed to handle the initialization and execution
    setup for Nextflow workflows. It covers the loading and validation of configuration
    files, YAML parameters, and environmental settings, which are essential for running
    workflow processes. This class supports additional features like custom annotation
    configuration and logging for monitoring workflow processes. Its primary purpose
    is to integrate setup functionality with Nextflow scripts and ensure validation
    for prerequisite parameters and file paths.

    Attributes:
        workflow_file (Path): The path to the Nextflow workflow file.
        workflow_dir (Path): The directory of the Nextflow workflow file.
        input_file (Path): The input file for the workflow.
        output_dir (Path): The directory where workflow output is stored.
        name (str): The unique identifier or name for the workflow instance.
        work_dir (Optional[Path]): The working directory for the workflow, automatically set up.
        logger (Logger): The logger instance for handling workflow logs.
        nf_config (Optional[Path]): The path to the main configuration file for the workflow.
        nf_config_content (Optional[dict]): Parsed content of the main configuration file.
        nfa_config (Optional[Path]): The path to the annotation configuration file.
        nfa_config_content (Optional[dict]): Parsed content of the annotation configuration file.
        params_file (Optional[Path]): The path to the YAML file containing workflow parameters.
        params_file_content (Optional[dict]): Parsed content of the YAML parameters file.
        workflow_dir_src (Path): The source directory for workflows within the package.
    """

    NXF_VERSION = "24.10.5"

    REQUIRED_PARAMS = (
        {  # these are the minimal yaml parameters required to run the workflow
            "bcftools_cmd": Path,
            "annotation_tool_cmd": str,
            "tool_version_command": str,
            "reference": Path,
            "chr_add": Path,
            "temp_dir": Path,
            "optional_checks": dict,
        }
    )

    """ Base class for Nextflow workflow operations for a single input file """

    def __init__(
        self,
        workflow: Path,
        input_file: Path,
        output_dir: Path,
        name: str,
        config_file: Optional[Path] = None,
        anno_config_file: Optional[Path] = None,
        params_file: Optional[Path] = None,
        verbosity: int = 0,
    ):

        self.workflow_file = Path(workflow).expanduser()
        self.workflow_dir = self.workflow_file.parent
        self.input_file = Path(input_file).expanduser()
        self.output_dir = Path(output_dir).expanduser()
        self.name = name
        self.work_dir: Optional[Path] = None
        # Set up logging
        log_file = Path(self.workflow_dir) / "workflow.log"
        self.logger = setup_logging(
            verbosity=verbosity,
            log_file=log_file if self.workflow_file.exists() else None,
        )

        if not self.workflow_file.exists():
            self.logger.error(f"Workflow file not found: {self.workflow_file}")
            raise FileNotFoundError(f"Workflow file not found: {self.workflow_file}")

        self.workflow_dir = self.workflow_file.parent
        # Convert Traversable to string to avoid type error
        workflow_resource = importlib.resources.files("vcfstash") / "workflow"
        self.workflow_dir_src = Path(str(workflow_resource))

        self.logger.info(f"Initializing Nextflow workflow in: {self.workflow_dir}")

        self.nf_config = self.nf_config_content = None
        if config_file:
            self.nf_config = Path(config_file).expanduser().resolve()
            if not self.nf_config.exists():
                self.logger.error(f"Nextflow config file not found: {self.nf_config}")
                raise FileNotFoundError(
                    f"Nextflow config file not found: {self.nf_config}"
                )
            self.nf_config_content = self.read_groovy_config(self.nf_config)
            self.validate_env_config()

        self.nfa_config = self.nfa_config_content = None
        if anno_config_file:
            self.nfa_config = Path(anno_config_file).expanduser().resolve()
            if not self.nfa_config.exists():
                self.logger.error(
                    f"Nextflow annotation config file not found: {self.nfa_config}"
                )
                raise FileNotFoundError(
                    f"Nextflow config file not found: {self.nfa_config}"
                )
            self.nfa_config_content = self.read_groovy_config(self.nfa_config)
            self.validate_annotation_config()

        self.params_file = self.params_file_content = None
        if params_file:
            self.params_file = Path(params_file).expanduser().resolve()
            if not self.params_file.exists():
                self.logger.error(f"Config yaml file not found: {self.params_file}")
                raise FileNotFoundError(
                    f"Config yaml file not found: {self.params_file}"
                )
            self.params_file_content = yaml.safe_load(self.params_file.read_text())
            self.validate_params_config()

        self.params_file = Path(params_file).expanduser() if params_file else None
        if self.params_file and not self.params_file.exists():
            self.logger.error(f"Parameters file not found: {self.params_file}")
            raise FileNotFoundError(f"Parameters file not found: {self.params_file}")

        self.work_dir = None
        self._setup_nextflow_skeleton()

    def validate_params_config(self):
        # Load and validate YAML parameters
        try:
            # Check for missing parameters
            missing_params = set(self.REQUIRED_PARAMS.keys()) - set(
                self.params_file_content.keys()
            )
            if missing_params:
                raise ValueError(
                    f"Missing required parameters in YAML: {', '.join(missing_params)}"
                )

            # Validate parameter types and paths
            for param, param_type in self.REQUIRED_PARAMS.items():
                value = self.params_file_content[param]
                if param == "optional_checks":
                    # Special handling for optional_checks
                    if not isinstance(value, dict):
                        raise TypeError(
                            f"Parameter {param} must be a dictionary/mapping"
                        )
                elif param_type == Path:
                    # Skip validation for paths with environment variables
                    if str(value).startswith("${"):
                        continue
                    path = Path(str(value)).expanduser()
                    if not path.exists():
                        raise FileNotFoundError(f"Path not found for {param}: {path}")
                elif not isinstance(value, param_type):
                    raise TypeError(
                        f"Parameter {param} must be of type {param_type.__name__}"
                    )

            # Validate optional checks against annotation config if available
            if hasattr(self, "nfa_config_content") and self.nfa_config_content:
                self._validate_optional_checks()

            self.logger.debug("YAML parameter validation successful")

        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing YAML file: {e}") from e
        except Exception as e:
            raise ValueError(f"Error validating YAML parameters: {e}") from e

    def _validate_optional_checks(self):
        """Validates that all keys in the optional_checks section of the YAML file
        have matching values in the annotation configuration's optional_checks section.
        """
        yaml_optional_checks = self.params_file_content.get("optional_checks", {})
        if not yaml_optional_checks:
            return  # No optional checks to validate

        # Check if the annotation config has an optional_checks section
        config_optional_checks = self.nfa_config_content.get("optional_checks", {})

        # Fall back to params section if optional_checks section doesn't exist (for backward compatibility)
        if not config_optional_checks:
            self.logger.debug(
                "No dedicated optional_checks section found, falling back to params section"
            )
            config_optional_checks = self.nfa_config_content.get("params", {})

        if not config_optional_checks:
            raise ValueError(
                "No section 'optional_checks' found in annotation config, please do not remove"
            )

        # Check each key-value pair from YAML optional_checks
        mismatches = []
        missing_keys = []

        for key, yaml_value in yaml_optional_checks.items():
            if key not in config_optional_checks:
                missing_keys.append(key)
                continue

            config_value = config_optional_checks[key]

            # Convert both to strings for comparison
            yaml_value_str = str(yaml_value)
            config_value_str = str(config_value)

            # Perform comparison (case-sensitive)
            if yaml_value_str != config_value_str:
                mismatches.append(
                    f"Value mismatch for '{key}': YAML has '{yaml_value_str}', config has '{config_value_str}'"
                )

        # Report any issues
        if missing_keys:
            self.logger.warning(
                f"Keys in optional_checks not found in config: {', '.join(missing_keys)}"
            )

        if mismatches:
            error_msg = "Optional checks validation failed:\n" + "\n".join(mismatches)
            self.logger.error(error_msg)
            raise ValueError(error_msg)

        self.logger.debug(
            f"Successfully validated {len(yaml_optional_checks)} optional checks"
        )

    def validate_env_config(self):
        """Validates the main environment configuration file (typically env_xxx.config) loaded into self.config_content.

        Checks for required parameters, valid paths, and proper configuration structure.

        Returns:
            bool: True if configuration is valid

        Raises:
            ValueError: If critical configuration issues are found
        """
        warnings = []

        # Check if config_content exists
        if not hasattr(self, "nf_config_content") or not self.nf_config_content:
            raise ValueError(
                f"No configuration for {self.nf_config} content loaded:\n{self.nf_config_content}"
            )

        # Extract params section
        config_params = self.nf_config_content.get("params", {})
        if config_params:
            err_msg = f"This config file {self.nf_config} should only contain a process section: {config_params}"
            self.logger.error(err_msg)
            raise ValueError(err_msg)

        # 1. Check required parameters
        required_params = [
            "bcftools_cmd",
            "chr_add",
            "annotation_tool_cmd",
            "tool_version_command",
        ]

        missing_params = [
            param for param in required_params if param not in config_params
        ]
        if missing_params:
            err_msg = f"Missing required parameters in environment nextflow config: {', '.join(missing_params)}"
            self.logger.error(err_msg)
            raise ValueError(err_msg)

        # 2. Check paths for existence
        file_paths = ["reference", "chr_add"]

        for path_param in file_paths:
            if path_param in config_params:
                path_str = config_params[path_param]
                if (
                    path_str
                    and not isinstance(path_str, bool)
                    and not path_str.startswith("${")
                ):
                    path = Path(path_str).expanduser().resolve()
                    if not path.exists():
                        warn_msg = f"Path defined in '{path_param}' does not exist: {str(path)}"
                        self.logger.warning(warn_msg)
                        warnings.append(warn_msg)

        # Print summary of validation
        if warnings:
            self.logger.warning(
                f"Configuration validation found {len(warnings)} warning(s)."
            )

        self.logger.info("Main configuration validation completed successfully.")

    def validate_annotation_config(self):
        """Validates the annotation configuration (annotation.config) loaded into self.anno_config_content.

        Checks for required parameters, annotation options, and proper configuration structure.

        Returns:
            bool: True if annotation configuration is valid

        Raises:
            ValueError: If critical configuration issues are found
        """
        # Check if anno_config_content exists
        if not hasattr(self, "nfa_config_content") or not self.nfa_config_content:
            raise ValueError(
                f"No annotation configuration content loaded:\n{self.nfa_config_content}"
            )

        # Extract params section
        anno_params = self.nfa_config_content.get("params", {})
        if not anno_params:
            err_msg = "No 'params' section found in annotation config file."
            self.logger.error(err_msg)
            raise ValueError(err_msg)

        # 1. Check required MD5 parameters
        required_params = [
            "must_contain_info_tag",
            "annotation_cmd",
            "required_tool_version",
        ]

        missing_params = [
            param for param in required_params if param not in anno_params
        ]
        if missing_params:
            err_msg = (
                f"Missing parameters in annotation config: {', '.join(missing_params)}"
            )
            self.logger.error(err_msg)
            raise ValueError(err_msg)

        self.logger.debug("Annotation configuration validation completed successfully.")

    def read_groovy_config(self, config_path: Path | str) -> dict:
        """Reads and parses a Groovy configuration file into a Python dictionary.

        Args:
            config_path (Path | str): Path to the Groovy configuration file.

        Returns:
            dict: Dictionary representation of the Groovy configuration file if it exists,
                  else empty dictionary.
        """
        try:
            config_file = Path(config_path)
            if not config_file.exists():
                self.logger.warning(f"Groovy config file not found: {config_path}")
                return {}

            content = config_file.read_text()
            return self._parse_groovy_content(content)
        except Exception as e:
            self.logger.error(f"Error reading Groovy config file {config_path}: {e}")
            return {}

    @staticmethod
    def _parse_groovy_content(content: str) -> dict:
        """Parses Groovy configuration content into a Python dictionary.

        Args:
            content (str): Content of the Groovy configuration file.

        Returns:
            dict: Dictionary representation of the Groovy configuration.
        """
        result: Dict[str, Any] = {}
        current_section = result
        section_stack = []

        # Track if we're in a list
        in_list = False
        current_key = None
        list_items = []

        # Split content into lines and process
        lines = content.split("\n")
        for line in lines:
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith("//"):
                continue

            # Check if we're in a multi-line list
            if in_list:
                if "]" in line:
                    # End of list
                    list_part = line.split("]")[0].strip()
                    if list_part:
                        # Add the last item if it's not empty
                        if list_part.endswith(","):
                            list_part = list_part[:-1].strip()
                        if list_part:
                            item = list_part.strip()
                            if (item.startswith('"') and item.endswith('"')) or (
                                item.startswith("'") and item.endswith("'")
                            ):
                                item = item[1:-1]
                            list_items.append(item)

                    if current_key is not None:
                        current_section[current_key] = list_items
                    in_list = False
                    list_items = []
                    current_key = None
                else:
                    # Continue adding to the list
                    if line:
                        item = line.strip()
                        if item.endswith(","):
                            item = item[:-1].strip()
                        if item:
                            if (item.startswith('"') and item.endswith('"')) or (
                                item.startswith("'") and item.endswith("'")
                            ):
                                item = item[1:-1]
                            list_items.append(item)
                continue

            # Handle section start
            section_match = re.match(r"(\w+)\s*\{", line)
            if section_match:
                section_name = section_match.group(1)
                if section_name not in current_section:
                    current_section[section_name] = {}
                section_stack.append(current_section)
                current_section = current_section[section_name]
                continue

            # Handle section end
            if line.startswith("}"):
                if section_stack:
                    current_section = section_stack.pop()
                continue

            # Handle key-value pairs
            kv_match = re.match(r"(\w+)\s*=\s*(.+)", line)
            if kv_match:
                key, value = kv_match.group(1), kv_match.group(2).strip()

                # Check for list start
                if value.startswith("[") and not value.endswith("]"):
                    # Start of a multi-line list
                    in_list = True
                    current_key = key
                    list_items = []

                    # Process first line of the list
                    first_item = value[1:].strip()
                    if first_item and not first_item.startswith("//"):
                        if first_item.endswith(","):
                            first_item = first_item[:-1].strip()
                        if first_item:
                            if (
                                first_item.startswith('"') and first_item.endswith('"')
                            ) or (
                                first_item.startswith("'") and first_item.endswith("'")
                            ):
                                first_item = first_item[1:-1]
                            list_items.append(first_item)
                    continue

                # Remove trailing commas
                if value.endswith(","):
                    value = value[:-1].strip()

                # Handle quoted strings
                if (value.startswith('"') and value.endswith('"')) or (
                    value.startswith("'") and value.endswith("'")
                ):
                    value = value[1:-1]
                # Handle numbers
                elif value.isdigit():
                    value = int(value)
                elif value.lower() == "true":
                    value = True
                elif value.lower() == "false":
                    value = False
                elif value.lower() == "null":
                    value = None
                # Handle single-line lists
                elif (
                    isinstance(value, str)
                    and value.startswith("[")
                    and value.endswith("]")
                ):
                    list_str = value[1:-1].strip()
                    if list_str:  # Non-empty list
                        items = []
                        # Simple split by comma for now
                        for item in re.split(r",\s*", list_str):
                            item = item.strip()
                            # Remove quotes from string items
                            if (item.startswith('"') and item.endswith('"')) or (
                                item.startswith("'") and item.endswith("'")
                            ):
                                item = item[1:-1]
                            items.append(item)
                        value = items
                    else:
                        value = []

                current_section[key] = value

        return result

    def _setup_nextflow_skeleton(self) -> None:
        """Sets up a skeleton .nextflow directory structure in workflow_dir and includes the jar file.

        Args:
            nxf_parent (Path): Path to the Nextflow parent dir
        """
        self.nxf_home = self.workflow_dir / ".nextflow"
        jar_dir = self.nxf_home / "framework" / self.NXF_VERSION
        jar_dir.mkdir(parents=True, exist_ok=True)

        # Create some other common subdirectories, however they should exist already as they are copied from workflow_dir_src at stash-init
        dirs = ["cache", "plugins", "plr"]
        for d in dirs:
            dir_path = Path(self.nxf_home / d)
            dir_path.mkdir(parents=True, exist_ok=True)

        jar_fl = jar_dir / f"nextflow-{self.NXF_VERSION}-one.jar"
        if not jar_fl.exists():
            self.logger.error(f"Nextflow JAR file not found: {jar_fl}")
            raise FileNotFoundError(f"Nextflow JAR file not found: {jar_fl}")

        # Create an empty history file
        with open(self.nxf_home / "history", "w"):
            pass

        self.logger.debug(f"Nextflow skeleton setup completed in {self.nxf_home}")

    def _get_temp_files(self) -> List[Path]:
        """Get list of temporary work directories"""
        if not self.work_dir:
            return []
        return [self.work_dir]

    def cleanup_work_dir(self) -> None:
        """Remove temporary work directory"""
        if not self.work_dir:
            return

        self.logger.debug("Cleaning up work directory")
        try:
            if self.work_dir.exists():
                self.logger.debug(f"Removing work directory: {self.work_dir}")
                shutil.rmtree(self.work_dir)
        except Exception as e:
            self.logger.warning(f"Failed to remove work directory {self.work_dir}: {e}")
        self.work_dir = None

    def warn_temp_files(self) -> None:
        """Print warning about existing temporary files"""
        temp_files = self._get_temp_files()
        if temp_files:
            self.logger.warning("Temporary files from failed run exist:")
            for path in temp_files:
                if path.exists():
                    self.logger.warning(f"- {path}")
            self.logger.warning("You may want to remove these files manually")

    def store_workflow_dag(self, run_dir: Path, cmd: List[str]) -> None:
        """Generate and store workflow DAG visualization"""
        try:
            self.logger.debug("Generating workflow DAG")
            try:
                del cmd[cmd.index("-with-trace")]
            except ValueError:
                pass

            dag_cmd = cmd + ["-preview", "-with-dag", str(run_dir / "flowchart.html")]
            subprocess.run(dag_cmd, check=True, capture_output=True)
            self.logger.info(f"Workflow DAG saved to: {run_dir}/flowchart.html")

        except subprocess.CalledProcessError as e:
            self.logger.warning(f"Failed to generate workflow DAG: {e}")

    def _expand_and_write_params(self, no_conf: bool) -> str:
        """Expand env vars in params YAML and write to a temp file."""
        if self.params_file is None:
            raise ValueError("params_file is required for _expand_and_write_params")
        with open(self.params_file) as f:
            content = f.read()
        expanded = os.path.expandvars(content)
        params_dict = yaml.safe_load(expanded)

        # Always set must_contain_info_tag if missing to default value
        if no_conf and "must_contain_info_tag" not in params_dict:
            params_dict["must_contain_info_tag"] = ""

        tmp = tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False)
        yaml.safe_dump(params_dict, tmp)
        tmp.close()
        return tmp.name

    def _create_work_dir(self, parent: Path, dirname: str = "work") -> None:
        """Create a temporary work directory for Nextflow"""
        self.work_dir = parent / dirname
        if self.work_dir is not None:
            if not self.work_dir.exists():
                self.work_dir.mkdir(parents=True, exist_ok=True)
                if self.logger:
                    self.logger.debug(f"Created work directory: {self.work_dir}")
            else:
                if self.logger:
                    self.logger.warning(
                        f"Work directory already exists: {self.work_dir}"
                    )

    def run(
        self,
        db_mode: str,
        nextflow_args: Optional[List[str]] = None,
        trace: bool = False,
        db_bcf: Optional[Path] = None,
        dag: bool = False,
        timeline: bool = False,
        report: bool = False,
        temp: Union[Path, str] = "/tmp",
    ) -> subprocess.CompletedProcess:
        """Run the Nextflow workflow."""
        # Create temporary work directory
        self._create_work_dir(self.output_dir)
        env = os.environ.copy()
        env["NXF_HOME"] = str(self.nxf_home)
        env["NXF_WORK"] = str(self.work_dir)
        env["NXF_TEMP"] = str(temp)
        env["NXF_VER"] = self.NXF_VERSION
        env["NXF_DISABLE_CHECK_LATEST"] = "1"
        env["NXF_OFFLINE"] = "true"
        if "VCFSTASH_ROOT" not in env:
            raise RuntimeError("VCFSTASH_ROOT environment variable is not set")
        # Set up the Nextflow executable
        nxf_exe = [
            "java",
            "-jar",
            str(
                Path(env["VCFSTASH_ROOT"])
                / f"vcfstash/workflow/.nextflow/framework/{self.NXF_VERSION}/nextflow-{self.NXF_VERSION}-one.jar"
            ),
        ]
        # this could also use the executable in PATH with nxf_exe = ['nextflow']

        # Clean Nextflow metadata
        subprocess.run(nxf_exe + ["clean", "-f"], check=False, cwd=self.output_dir)

        # Global options (applied before the command)
        global_opts = [
            "-log",
            str(
                self.output_dir / ".nextflow.log"
            ),  # todo: try '-bg' / evaluate '-c' vs '-C'
        ]

        if self.nfa_config:
            global_opts.append("-C")
            global_opts.append(str(self.nfa_config))

        # Run-specific options (applied after the "run" command)
        run_opts = [
            str(self.workflow_file),  # project name or repo URL
            "-offline",
            "--input",
            str(self.input_file),
            "--output",
            str(self.output_dir),
            "--db_mode",
            db_mode,
            "-w",
            str(self.work_dir),
            "-name",
            f"vcfstash_{self.name}",
            "-ansi-log",
            "true",  # enable colored output
        ]

        # Add params file if specified
        params_tfile = None
        if self.params_file:
            params_tfile = self._expand_and_write_params(
                no_conf=self.nfa_config is None
            )
            run_opts.extend(["-params-file", params_tfile])

        if db_bcf:
            run_opts.extend(["--db_bcf", str(db_bcf)])

        if trace:
            run_opts.extend(
                ["-with-trace", str(self.output_dir / f"{self.name}_trace.txt")]
            )

        if dag:
            run_opts.extend(
                ["-with-dag", str(self.output_dir / f"{self.name}_flowchart.html")]
            )

        if report:
            run_opts.extend(
                ["-with-report", str(self.output_dir / f"{self.name}_report.html")]
            )

        if timeline:
            run_opts.extend(
                ["-with-timeline", str(self.output_dir / f"{self.name}_timeline.txt")]
            )

        # Add additional arguments
        if nextflow_args:
            run_opts.extend(nextflow_args)

        # Assemble the final command list:
        cmd = nxf_exe + global_opts + ["run"] + run_opts

        self.logger.debug(f"Running command: {' '.join(map(str, cmd))}")

        try:
            # Run without capturing output to show live progress
            result = subprocess.run(cmd, check=True, cwd=self.output_dir, env=env)
            return result
        except subprocess.CalledProcessError as e:
            self.warn_temp_files()
            self.logger.error(
                f"Workflow execution failed with exit code: {e.returncode}"
            )
            raise RuntimeError(
                f"Workflow execution failed with exit code: {e.returncode}"
            ) from e
        finally:
            if params_tfile:
                os.unlink(params_tfile)


class VCFDatabase:
    """Base class for performing operations on VCF (Variant Call Format) databases.

    Provides tools for managing, validating, and manipulating VCF/BCF files including
    index validation, header validation, INFO field parsing, and genome reference checks.
    Also includes utilities for logging and configuration setup designed for bioinformatics
    workflow management.

    Attributes:
        TRANSCRIPT_KEYS: List of predefined transcript consequence keys used in
                         INFO annotations of VCF files.
    """

    TRANSCRIPT_KEYS = [
        "SYMBOL",
        "Feature",
        "Consequence",
        "HGVS_OFFSET",
        "HGVSc",
        "HGVSp",
        "IMPACT",
        "DISTANCE",
        "PICK",
        "VARIANT_CLASS",
    ]
    """Base class for VCF database operations"""

    def __init__(self, db_path: Path, verbosity: int, debug: bool):
        self.debug = debug
        self.stashed_output = StashOutput(str(db_path))
        self.stash_path = self.stashed_output.root_dir
        self.stash_name = self.stash_path.name
        self.blueprint_dir = self.stash_path / "blueprint"
        self.workflow_dir = self.stash_path / "workflow"
        self.stash_dir = self.stash_path / "stash"
        self.workflow_dir_src = self.stashed_output.workflow_src_dir
        self.blueprint_bcf = self.blueprint_dir / "vcfstash.bcf"

        self.info_file = self.blueprint_dir / "sources.info"
        self.verbosity = verbosity if not debug else 2

        if self.info_file.exists():
            with open(self.info_file) as f:
                self.db_info = json.load(f)
        else:
            self.db_info = {"input_files": []}
        self.logger: Optional[Logger] = None

    def connect_loggers(self, logger_name: str = "vcfdb") -> Logger:
        # Set up central logging
        log_file = self.stash_path / f"{logger_name}.log"
        return setup_logging(verbosity=self.verbosity, log_file=log_file)

    def setup_config(self, config_file: Path | str, config_name: str) -> Path:
        config_path = Path(config_file) if isinstance(config_file, str) else config_file
        config_path = config_path.expanduser().resolve()
        shutil.copyfile(config_path, self.workflow_dir / config_name)
        return self.workflow_dir / config_name

    def ensure_indexed(self, file_path: Path) -> None:
        """Ensure input file has an index file (CSI or TBI)"""
        if self.logger:
            self.logger.debug(f"Checking index for: {file_path}")
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        csi_index = Path(str(file_path) + ".csi")
        tbi_index = Path(str(file_path) + ".tbi")

        if not (csi_index.exists() or tbi_index.exists()):
            raise RuntimeError(
                f"No index found for {file_path}. Use bcftools index for BCF/compressed VCF "
                "or tabix for compressed VCF files."
            )

    def validate_bcf_header(
        self, bcf_path: Path, norm: bool = True
    ) -> Tuple[bool, Optional[str]]:
        """Validate BCF header format"""
        if self.logger:
            self.logger.debug(f"Validating BCF header: {bcf_path}")
        result = validate_bcf_header(bcf_path, norm)
        if not result[0] and self.logger:
            self.logger.error(f"BCF header validation failed: {result[1]}")
        return result

    def parse_vcf_info(self, vcf_data: list) -> list:
        """Parses VCF INFO field and expands transcript consequences."""

        def convert_vcfstr(value):
            if value is None or value == "":
                return None
            try:
                if isinstance(value, str) and (
                    "." in value or "e" in value or "E" in value
                ):
                    return float(value)
                return int(value)
            except ValueError:
                return value

        expanded_data: List[Dict[str, Any]] = []
        for tn, transcript in enumerate(vcf_data):
            expanded_data.append({})
            for entry in self.TRANSCRIPT_KEYS:
                if entry not in vcf_data[tn]:
                    raise ValueError(
                        f"Did not find key={entry} {transcript=} in CSQ INFO tag"
                    )
                if entry == "PICK":
                    expanded_data[tn][entry] = vcf_data[tn][entry] == "1"
                else:
                    expanded_data[tn][entry] = convert_vcfstr(vcf_data[tn][entry])

        return expanded_data

    def validate_vcf_reference(
        self,
        vcf_path: Path,
        ref_fasta: Path,
        chr_add_path: Optional[Path] = None,
        variants_to_check: int = 3,
    ) -> Tuple[bool, Optional[str]]:
        r"""Validates a VCF/BCF file against a reference genome with optional chromosome name mapping.

        This method performs two checks:
        1. Verifies that all contigs in the VCF/BCF index are found in the reference FASTA index
           (with optional mapping through chr_add file)
        2. Checks a configurable number of variants from different chromosomes to ensure their reference
           alleles match the reference genome

        Args:
            vcf_path (Path): Path to the VCF/BCF file to validate
            ref_fasta (Path): Path to the reference genome FASTA file
            chr_add_path (Path, optional): Path to chromosome mapping file. Format: two columns
                                           separated by tab, e.g., "1\tchr1". Default is None.
            variants_to_check (int, optional): Number of variants to check. Defaults to 3.

        Returns:
            Tuple[bool, Optional[str]]: A tuple containing:
                - bool: True if validation passed, False otherwise
                - Optional[str]: Error message if validation failed, None otherwise
        """
        if self.logger:
            self.logger.debug(f"Validating VCF/BCF file against reference: {vcf_path}")

        # Check if files exist
        if not vcf_path.exists():
            return False, f"VCF/BCF file not found: {vcf_path}"

        if not ref_fasta.exists():
            return False, f"Reference genome not found: {ref_fasta}"

        # Check for index files
        vcf_index_csi = Path(str(vcf_path) + ".csi")
        vcf_index_tbi = Path(str(vcf_path) + ".tbi")
        ref_index = Path(str(ref_fasta) + ".fai")

        if not (vcf_index_csi.exists() or vcf_index_tbi.exists()):
            return (
                False,
                f"VCF/BCF index not found for {vcf_path}. Use bcftools index or tabix.",
            )

        if not ref_index.exists():
            return (
                False,
                f"Reference index not found for {ref_fasta}. Use samtools faidx.",
            )

        try:
            # Load reference contigs from .fai
            ref_contigs = set()
            with open(ref_index, "r") as f:
                for line in f:
                    contig = line.split("\t")[0]
                    ref_contigs.add(contig)

            if self.logger:
                self.logger.debug(
                    f"Found {len(ref_contigs)} contigs in reference index"
                )

            # Get VCF/BCF contigs from index instead of header
            # This is important as the actual data may contain chromosomes not declared in the header
            vcf_contigs = set()
            vcf_index_file = vcf_index_csi if vcf_index_csi.exists() else vcf_index_tbi

            # Use bcftools to extract chromosomes from the index
            cmd = ["bcftools", "index", "--stats", str(vcf_index_file)]
            result = subprocess.run(cmd, capture_output=True, text=True)

            # Parse chromosomes from output
            for line in result.stdout.splitlines():
                if "[" not in line and "]" not in line:  # Skip header/footer lines
                    parts = line.strip().split()
                    if parts and parts[0]:  # Ensure we have a chromosome name
                        vcf_contigs.add(parts[0])

            if self.logger:
                self.logger.debug(
                    f"Found {len(vcf_contigs)} unique chromosomes in VCF index"
                )

            # Load chromosome mapping if provided
            chrom_map = {}
            if chr_add_path and Path(chr_add_path).exists():
                if self.logger:
                    self.logger.debug(f"Loading chromosome mapping from {chr_add_path}")
                with open(chr_add_path, "r") as f:
                    for line in f:
                        if line.strip() and not line.startswith("#"):
                            cols = line.strip().split("\t")
                            if len(cols) >= 2:
                                src, dest = cols[0], cols[1]
                                chrom_map[src] = dest
                                # Also add reverse mapping
                                chrom_map[dest] = src
                if self.logger:
                    self.logger.debug(f"Loaded {len(chrom_map)} chromosome mappings")

            # Check if all VCF contigs are in reference (with mapping)
            missing_contigs = set()
            for contig in vcf_contigs:
                # Check direct match
                if contig in ref_contigs:
                    continue

                # Check mapped contig name
                mapped_contig = chrom_map.get(contig)
                if mapped_contig and mapped_contig in ref_contigs:
                    continue

                # Try common transformations if no mapping found
                if not mapped_contig:
                    # Try with/without 'chr' prefix
                    if contig.startswith("chr") and contig[3:] in ref_contigs:
                        continue
                    if f"chr{contig}" in ref_contigs:
                        continue

                    # Try MT/M variations for mitochondrial DNA
                    if contig in ("MT", "M", "chrM", "chrMT") and (
                        "chrM" in ref_contigs
                        or "M" in ref_contigs
                        or "MT" in ref_contigs
                    ):
                        continue

                # If we get here, the contig is missing
                missing_contigs.add(contig)

            if missing_contigs:
                return (
                    False,
                    f"Contigs in VCF not found in reference (even after mapping): {', '.join(missing_contigs)}",
                )

            # Check reference alleles for a subset of variants
            checked_variants = 0
            checked_chroms = set()
            mismatches = []

            # Open reference FASTA
            ref_fasta_file = pysam.FastaFile(str(ref_fasta))

            # Open VCF file
            vcf = pysam.VariantFile(str(vcf_path))

            # Iterate through variants
            for variant in vcf.fetch():
                # Skip if we've already checked this chromosome
                if variant.chrom in checked_chroms:
                    continue

                # Add chromosome to checked set
                checked_chroms.add(variant.chrom)

                # Get reference sequence - try direct match first
                try:
                    chr_to_use = variant.chrom

                    # If direct access fails, try mapped version
                    if chr_to_use not in ref_fasta_file.references:
                        if (
                            chr_to_use in chrom_map
                            and chrom_map[chr_to_use] in ref_fasta_file.references
                        ):
                            chr_to_use = chrom_map[chr_to_use]
                        # Try with/without 'chr' prefix
                        elif (
                            chr_to_use.startswith("chr")
                            and chr_to_use[3:] in ref_fasta_file.references
                        ):
                            chr_to_use = chr_to_use[3:]
                        elif f"chr{chr_to_use}" in ref_fasta_file.references:
                            chr_to_use = f"chr{chr_to_use}"

                    # Skip if variant.ref is None
                    if variant.ref is None:
                        continue

                    ref_seq = ref_fasta_file.fetch(
                        chr_to_use, variant.pos - 1, variant.pos - 1 + len(variant.ref)
                    )

                    # Compare reference allele (ensuring both are not None)
                    if ref_seq is not None and variant.ref.upper() != ref_seq.upper():
                        mismatches.append(
                            f"{variant.chrom}:{variant.pos} (VCF: {variant.ref}, FASTA: {ref_seq})"
                        )

                    checked_variants += 1

                    # Stop if we've checked enough variants
                    if checked_variants >= variants_to_check:
                        break

                except Exception as e:
                    if self.logger:
                        self.logger.warning(
                            f"Error checking variant at {variant.chrom}:{variant.pos}: {e}"
                        )

            # Check if we found any variants to check
            if checked_variants == 0:
                return False, "Could not find any variants to check reference alleles"

            # Report mismatches
            if mismatches:
                return (
                    False,
                    f"Reference allele mismatches found: {', '.join(mismatches)}",
                )

            if self.logger:
                self.logger.info(
                    f"Successfully validated {checked_variants} variants against reference"
                )
            return True, None

        except Exception as e:
            return False, f"Error validating VCF against reference: {e}"

    @staticmethod
    def _copy_workflow_srcfiles(
        source: Path, destination: Path, skip_config: bool = False
    ) -> None:
        """Copy workflow files from source to destination, optionally skipping nextflow.config.

        Parameters:
            destination (Path): Destination directory for copying.
            skip_config (bool): If True, does not copy nextflow.config. Default is False.
        """
        try:
            # Ensure the source directory exists
            if not source.exists():
                raise RuntimeError(
                    f"Source workflow directory does not exist: {source}"
                )

            # Ensure the destination directory exists
            if not destination.exists():
                raise RuntimeError(
                    f"Destination directory does not exist: {destination}"
                )

            # Copy entire directory structure except *.config files
            if skip_config:
                shutil.copytree(
                    source,
                    destination,
                    ignore=shutil.ignore_patterns("*.config"),
                    dirs_exist_ok=True,
                )
            else:
                shutil.copytree(source, destination, dirs_exist_ok=True)

            # Verify copy worked by ensuring destination isn't empty
            if not list(destination.glob("*")):
                raise RuntimeError(f"Failed to copy workflow files to {destination}")

        except Exception as e:
            raise RuntimeError(f"Error copying workflow files: {str(e)}") from e
