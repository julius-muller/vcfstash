import json
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

import yaml

from vcfstash.database.base import NextflowWorkflow, VCFDatabase
from vcfstash.utils.validation import compute_md5


class DatabaseInitializer(VCFDatabase):
    """Class for initializing and managing a VCF database.

    The `DatabaseInitializer` class is designed to set up and initialize a VCF database by
    processing input VCF/BCF files, managing configurations, and running associated workflows.
    It provides tools for validating inputs, ensuring output directory structure, and creating
    the database while handling workflow execution and logging.

    Attributes:
        input_file: Path to the input VCF/BCF file.
        config_file: Path to the Nextflow configuration file.
        output_dir: Path to the output directory for database initialization.
        verbosity: Logging verbosity level.
        force: Boolean to indicate whether to overwrite existing database files.
        debug: Boolean to enable or disable debug mode.
        logger: Logger instance for logging workflow steps and errors.
        config_yaml: Path to the YAML configuration file.
        nx_workflow: Instance of the associated workflow manager.
    """

    def __init__(
        self,
        input_file: Path | str,
        params_file: Path | str,
        config_file: Optional[Path | str] = None,
        output_dir: Path | str = Path("."),
        verbosity: int = 0,
        force: bool = False,
        debug: bool = False,
    ) -> None:
        """Initialize the database creator.

        Args:
            input_file: Path to input BCF/VCF file (required)
            config_file: Path to Nextflow config file (required)
            output_dir: Output directory (default: current directory)
            verbosity: Logging verbosity level (0=WARNING, 1=INFO, 2=DEBUG)
            force: Force overwrite of existing database (default: False)

        """
        # Initialize the parent class
        super().__init__(
            Path(output_dir) if isinstance(output_dir, str) else output_dir,
            verbosity,
            debug,
        )
        self._setup_stash(force=force)
        self.logger = self.connect_loggers()

        # self.validate_label(name)
        self.input_file = Path(input_file).expanduser().resolve()

        self._copy_workflow_srcfiles(
            source=self.workflow_dir_src,
            destination=self.workflow_dir,
            skip_config=True,
        )

        self.config_file = None
        if config_file:
            self.config_file = self.workflow_dir / "init.config"
            config_path = (
                Path(config_file) if isinstance(config_file, str) else config_file
            )
            shutil.copyfile(config_path.expanduser().resolve(), self.config_file)

        self.config_yaml = self.workflow_dir / "init.yaml"
        shutil.copyfile(Path(params_file).expanduser().resolve(), self.config_yaml)

        # Initialize NextflowWorkflow
        if self.logger:
            self.logger.info("Initializing Nextflow workflow...")
        params_path = Path(params_file) if isinstance(params_file, str) else params_file
        self.nx_workflow = NextflowWorkflow(
            input_file=self.input_file,
            output_dir=self.blueprint_dir,
            name="init",
            workflow=self.workflow_dir / "main.nf",
            config_file=self.config_file,
            params_file=params_path,
            verbosity=self.verbosity,
        )

        self._validate_inputs()

        # Log initialization parameters
        if self.logger:
            self.logger.info(f"Initializing database: {self.stash_name}")
            self.logger.debug(f"Input file: {self.input_file}")
            self.logger.debug(f"Output directory: {self.blueprint_dir}")
            self.logger.debug(f"Config file: {self.config_file}")

    def _setup_stash(self, force: bool) -> None:
        # Remove destination directory if it exists to ensure clean copy
        if self.stashed_output.root_dir.exists():
            if (
                self.stashed_output.validate_structure()
            ):  # we dont want to remove a random dir....
                if force:
                    print(
                        f"Stash directory already exists, removing: {self.stashed_output.root_dir}"
                    )
                    shutil.rmtree(self.stashed_output.root_dir)
                else:
                    raise FileExistsError(
                        f"Output directory already exists: {self.stashed_output.root_dir}\nIf intended, use --force to overwrite."
                    )
            else:
                raise FileExistsError(
                    f"Output directory with an invalid structure detected: {self.stashed_output.root_dir}"
                )

        print(f"Creating stash structure: {self.stashed_output.root_dir}")
        self.stashed_output.create_structure()

    def initialize(self) -> None:
        """Initialize new VCF database
        self = DatabaseInitializer(name='nftest', input_file=Path('tests/data/nodata/gnomad_test.bcf'),
        config_file=Path('workflow/nextflow.config'), output_dir=Path('~/tmp/test'), verbosity=2, force=True)
        self.workflow_dir
        self.output_dir
        self.input_file
        self._validate_inputs()
        self._create_database()
        """
        if not self.input_file.exists():
            if self.logger:
                self.logger.error(f"Input BCF file does not exist: {self.input_file}")
            raise FileNotFoundError("Input BCF file does not exist.")

        if self.blueprint_bcf.exists():
            if self.logger:
                self.logger.error(
                    f"Output database already exists: {self.blueprint_bcf}"
                )
            raise FileExistsError("Output database already exists.")

        self._create_database()

    def _validate_inputs(self) -> None:
        """Validate input files, directories, and YAML parameters"""
        if self.logger:
            self.logger.debug("Validating inputs")

        # Check if output database already exists
        if self.blueprint_bcf and self.blueprint_bcf.exists():
            msg = f"Database file already exists: {self.blueprint_bcf}"
            if self.logger:
                self.logger.error(msg)
            raise FileExistsError(msg)

        # Check input VCF/BCF if provided
        if self.input_file:
            if not self.input_file.exists():
                msg = f"Input VCF/BCF file not found: {self.input_file}"
                if self.logger:
                    self.logger.error(msg)
                raise FileNotFoundError(msg)
            self.ensure_indexed(self.input_file)

        # Validate VCF reference
        try:
            # Load YAML file to get reference path
            params_yaml = yaml.safe_load(
                Path(self.config_yaml).expanduser().resolve().read_text()
            )

            # Get reference path from YAML
            reference_path_str = params_yaml.get("reference")
            if not reference_path_str:
                if self.logger:
                    self.logger.error("Reference path not found in YAML file")
                raise ValueError("Reference path not found in YAML file")

            # Expand environment variables in reference path
            if "${VCFSTASH_ROOT}" in reference_path_str:
                vcfstash_root = os.environ.get("VCFSTASH_ROOT")
                if not vcfstash_root:
                    if self.logger:
                        self.logger.error("VCFSTASH_ROOT environment variable not set")
                    raise ValueError("VCFSTASH_ROOT environment variable not set")
                reference_path_str = reference_path_str.replace(
                    "${VCFSTASH_ROOT}", vcfstash_root
                )

            reference_path = Path(reference_path_str).expanduser().resolve()

            # Validate VCF reference
            if self.logger:
                self.logger.info(f"Validating VCF reference: {reference_path}")
            valid, error_msg = self.validate_vcf_reference(
                self.input_file, reference_path
            )

            if not valid:
                if self.logger:
                    self.logger.error(f"VCF reference validation failed: {error_msg}")
                raise ValueError(f"VCF reference validation failed: {error_msg}")

            if self.logger:
                self.logger.info("VCF reference validation successful")
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during VCF reference validation: {e}")
            raise RuntimeError(f"Error during VCF reference validation: {e}") from e

        if self.logger:
            self.logger.debug("Input validation successful")

    def _create_database(self) -> None:
        """Create and initialize the database"""
        if self.logger:
            self.logger.info(
                "Creating database from normalized and annotated variants..."
            )

        self.stashed_output.validate_structure()

        try:
            db_info: Dict[str, Any] = {
                "name": self.stash_name,
                "created": datetime.now().isoformat(),
                "input_files": [],  # List of input file info dictionaries
            }

            if not self.input_file:
                raise ValueError("Input file is required for database initialization")

            try:
                input_md5 = compute_md5(self.input_file)
                db_info["input_files"].append(
                    {
                        "path": str(self.input_file),
                        "md5": input_md5,
                        "added": datetime.now().isoformat(),
                    }
                )
            except Exception as e:
                if self.logger:
                    self.logger.error("Failed to compute MD5 checksum for input file.")
                raise RuntimeError(f"MD5 computation failed: {e}") from e

            # Run the workflow in database mode
            start_time = datetime.now()
            if self.logger:
                self.logger.info("Starting the workflow execution...")
            self.nx_workflow.run(
                db_mode="stash-init", trace=True, dag=True, report=True
            )
            if self.logger:
                self.logger.info("Workflow execution completed.")
            duration = datetime.now() - start_time

            # Save database info
            with open(self.info_file, "w") as f:
                json.dump(db_info, f, indent=2)

            if self.logger:
                self.logger.info("Database creation completed successfully.")

            # Verify that required files in the database are present
            if not self.info_file.exists():
                msg = f"Database info file missing: {self.info_file}"
                if self.logger:
                    self.logger.error(msg)
                raise FileNotFoundError(msg)
            if self.logger:
                self.logger.info(f"- Created at: {db_info['created']}")
                self.logger.info(f"- Input file: {self.input_file}")
                self.logger.info(f"- Output file: {self.blueprint_bcf}")
                self.logger.info(f"- Input MD5: {input_md5}")
                self.logger.info(f"- Processing time: {duration.total_seconds():.2f}s")

            if not self.debug:
                self.nx_workflow.cleanup_work_dir()

        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during database creation: {e}")
            raise RuntimeError(f"Error during database creation: {e}") from e
