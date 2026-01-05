# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""Abstract base class for workflow execution backends.

This module defines the interface that all workflow backends must implement.
"""

from abc import ABC, abstractmethod
from logging import Logger
from pathlib import Path
from typing import List, Optional, Union
import subprocess


class WorkflowBase(ABC):
    """Abstract base class for workflow execution backends.

    All workflow implementations (WorkflowManager) must
    inherit from this class and implement the required methods.

    Attributes:
        workflow_file (Path): The path to the workflow file.
        workflow_dir (Path): The directory of the workflow file.
        input_file (Path): The input file for the workflow.
        output_dir (Path): The directory where workflow output is stored.
        name (str): The unique identifier or name for the workflow instance.
        verbosity (int): Verbosity level for logging (0=quiet, 1=info, 2=debug).
        work_dir (Optional[Path]): The working directory for the workflow.
        logger (Logger): The logger instance for handling workflow logs.
        nf_config (Optional[Path]): The path to the main configuration file.
        nf_config_content (Optional[dict]): Parsed content of the main config.
        nfa_config (Optional[Path]): The path to the annotation configuration file.
        nfa_config_content (Optional[dict]): Parsed content of annotation config.
        params_file (Optional[Path]): The path to the YAML parameters file.
        params_file_content (Optional[dict]): Parsed content of the params file.
    """

    def __init__(
        self,
        input_file: Path,
        output_dir: Path,
        name: str,
        workflow: Path | None = None,
        anno_config_file: Optional[Path] = None,
        params_file: Optional[Path] = None,
        verbosity: int = 0,
    ):
        """Initialize the workflow backend.

        Args:
            workflow: Path to the workflow file (e.g., main.nf)
            input_file: Path to the input VCF/BCF file
            output_dir: Directory where output files will be stored
            name: Unique name for this workflow instance
            anno_config_file: Optional path to annotation configuration file
            params_file: Optional path to YAML parameters file
            verbosity: Verbosity level (0=quiet, 1=info, 2=debug)
        """
        if workflow is None:
            # Pure Python backend doesn't need a workflow script; use a placeholder
            self.workflow_file = Path(output_dir) / "workflow.stub"
        else:
            self.workflow_file = Path(workflow).expanduser()
        self.workflow_dir = self.workflow_file.parent
        self.input_file = Path(input_file).expanduser()
        self.output_dir = Path(output_dir).expanduser()
        self.name = name
        self.verbosity = verbosity
        self.work_dir: Optional[Path] = None

        # Configuration attributes (initialized by subclasses)
        self.nf_config: Optional[Path] = None
        self.nf_config_content: Optional[dict] = None
        self.nfa_config: Optional[Path] = None
        self.nfa_config_content: Optional[dict] = None
        self.params_file: Optional[Path] = None
        self.params_file_content: Optional[dict] = None
        self.logger: Optional[Logger] = None

    @abstractmethod
    def run(
        self,
        db_mode: str,
        trace: bool = False,
        db_bcf: Optional[Path] = None,
        dag: bool = False,
        timeline: bool = False,
        report: bool = False,
        temp: Union[Path, str] = "/tmp",
    ) -> subprocess.CompletedProcess:
        """Execute the workflow.

        This is the main entry point for running the workflow. Different
        db_mode values trigger different workflow logic:
        - 'blueprint-init': Create initial cache blueprint
        - 'blueprint-extend': Add variants to existing blueprint
        - 'cache-build': Annotate the blueprint
        - 'annotate': Annotate sample using cache
        - 'annotate-nocache': Annotate sample without cache

        Args:
            db_mode: The workflow mode (see above)
            trace: Whether to generate trace file
            db_bcf: Path to database BCF file (for blueprint-extend, annotate modes)
            dag: Whether to generate DAG visualization (may not be supported)
            timeline: Whether to generate timeline (may not be supported)
            report: Whether to generate report (may not be supported)
            temp: Temporary directory for intermediate files

        Returns:
            subprocess.CompletedProcess with execution results

        Raises:
            RuntimeError: If workflow execution fails
            ValueError: If invalid db_mode or missing required parameters
        """
        pass

    @abstractmethod
    def cleanup_work_dir(self) -> None:
        """Remove temporary work directory.

        This method is called to clean up temporary files after workflow
        execution. It should safely remove the work_dir if it exists.
        """
        pass

    def _get_temp_files(self) -> List[Path]:
        """Get list of temporary work directories.

        Returns:
            List of Path objects for temporary directories
        """
        if not self.work_dir:
            return []
        return [self.work_dir]

    def warn_temp_files(self) -> None:
        """Print warning about existing temporary files.

        This is called when workflow execution fails and temporary files
        were left behind. It informs the user about manual cleanup.
        """
        temp_files = self._get_temp_files()
        if temp_files and self.logger:
            self.logger.warning("Temporary files from failed run exist:")
            for path in temp_files:
                if path.exists():
                    self.logger.warning(f"- {path}")
            self.logger.warning("You may want to remove these files manually")
