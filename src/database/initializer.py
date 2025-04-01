import json
import shutil
from pathlib import Path
from datetime import datetime
from src.database.base import VEPDatabase, NextflowWorkflow
from src.database.outputs import StashOutput
from src.utils.validation import compute_md5


class DatabaseInitializer(VEPDatabase):
    """Handles database initialization"""
    def __init__(self, input_file: Path | str, config_file: Path | str, output_dir: Path | str = Path("."),
                 verbosity: int = 0, force: bool = False) -> None:
        """Initialize the database creator.

        Args:
            input_file: Path to input BCF/VCF file (required)
            config_file: Path to Nextflow config file (required)
            output_dir: Output directory (default: current directory)
            verbosity: Logging verbosity level (0=WARNING, 1=INFO, 2=DEBUG)
            force: Force overwrite of existing database (default: False)

        """

        # Initialize the parent class
        super().__init__(output_dir,  verbosity)
        self._setup_stash(force=force)
        self.logger = self.connect_loggers()

        # self.validate_label(name)
        self.input_file = Path(input_file).expanduser().resolve()

        self._copy_workflow_srcfiles(source=self.workflow_dir_src, destination=self.workflow_dir, skip_config=True)

        self.config_file = self.workflow_dir / 'init_nextflow.config'
        shutil.copyfile(config_file.expanduser().resolve(), self.config_file)

        # Initialize NextflowWorkflow
        self.logger.info("Initializing Nextflow workflow...")
        self.nx_workflow = NextflowWorkflow(
            input_file=self.input_file,
            output_dir=self.blueprint_dir,
            name='init',
            workflow=self.workflow_dir / "main.nf",
            config_file=self.config_file,
            verbosity=self.verbosity
        )

        # Log initialization parameters
        self.logger.info(f"Initializing database: {self.stash_name}")
        self.logger.debug(f"Input file: {self.input_file}")
        self.logger.debug(f"Output directory: {self.blueprint_dir}")
        self.logger.debug(f"Config file: {self.config_file}")


    def _setup_stash(self, force:bool) -> None:
        # Remove destination directory if it exists to ensure clean copy
        if self.stashed_output.root_dir.exists():
            if self.stashed_output.validate_structure():  # we dont want to remove a random dir....
                if force:
                    print(f"Stash directory already exists, removing: {self.stashed_output.root_dir}")
                    shutil.rmtree(self.stashed_output.root_dir)
                else:
                    raise FileExistsError(
                        f"Output directory already exists: {self.stashed_output.root_dir}\nIf intended, use --force to overwrite.")
            else:
                raise FileNotFoundError(
                    f"Output directory must not exist if --force is not set and a valid stash directory: {self.stashed_output.root_dir}")

        print(f"Creating stash structure: {self.stashed_output.root_dir}")
        self.stashed_output.create_structure()

    def initialize(self) -> None:
        """Initialize new VEP database
        self = DatabaseInitializer(name='nftest', input_file=Path('tests/data/nodata/gnomad_test.bcf'),
        config_file=Path('workflow/nextflow.config'), output_dir=Path('~/tmp/test'), verbosity=2, force=True)
        self.workflow_dir
        self.output_dir
        self.input_file
        self._validate_inputs()
        self._create_database()
        """

        if not self.input_file.exists():
            self.logger.error(f"Input BCF file does not exist: {self.input_file}")
            raise FileNotFoundError("Input BCF file does not exist.")

        if self.blueprint_bcf.exists():
            self.logger.error(f"Output database already exists: {self.blueprint_bcf}")
            raise FileExistsError("Output database already exists.")

        self.logger.info("Starting database initialization")
        self._validate_inputs()
        self._create_database()
        self.logger.info("Database initialization completed successfully")

    def _validate_inputs(self) -> None:
        """Validate input files and directories"""
        self.logger.debug("Validating inputs")

        # Check if output database already exists
        if self.blueprint_bcf and self.blueprint_bcf.exists():
            msg = f"Database file already exists: {self.blueprint_bcf}"
            self.logger.error(msg)
            raise FileExistsError(msg)

        # Check input VCF/BCF if provided
        if self.input_file:
            if not self.input_file.exists():
                msg = f"Input VCF/BCF file not found: {self.input_file}"
                self.logger.error(msg)
                raise FileNotFoundError(msg)
            self.ensure_indexed(self.input_file)

        self.logger.debug("Input validation successful")

    def _create_database(self) -> None:
        """Create and initialize the database"""
        self.logger.info("Creating database from normalized and annotated variants...")


        self.logger.info(f"Workflow files copied to: {self.workflow_dir}")

        try:
            db_info = {
                "name": self.stash_name,
                "created": datetime.now().isoformat(),
                "input_files": []
            }

            if not self.input_file:
                raise ValueError("Input file is required for database initialization")

            try:
                input_md5 = compute_md5(self.input_file)
                db_info["input_files"].append({
                    "path": str(self.input_file),
                    "md5": input_md5,
                    "added": datetime.now().isoformat()
                })
            except Exception as e:
                self.logger.error("Failed to compute MD5 checksum for input file.")
                raise RuntimeError(f"MD5 computation failed: {e}")

            # Run the workflow in database mode
            start_time = datetime.now()
            self.logger.info("Starting the workflow execution...")
            self.nx_workflow.run(
                db_mode="stash-init",
                trace=True,
                dag=True,
                report=True
            )
            self.logger.info("Workflow execution completed.")
            duration = datetime.now() - start_time

            # Save database info
            with open(self.info_file, "w") as f:
                json.dump(db_info, f, indent=2)

            self.logger.info("Database creation completed successfully.")

            # Verify that required files in the database are present
            if not self.info_file.exists():
                msg = f"Database info file missing: {self.info_file}"
                self.logger.error(msg)
                raise FileNotFoundError(msg)
            self.logger.info(f"- Created at: {db_info['created']}")
            self.logger.info(f"- Input file: {self.input_file}")
            self.logger.info(f"- Output file: {self.blueprint_bcf}")
            self.logger.info(f"- Input MD5: {input_md5}")
            self.logger.info(f"- Processing time: {duration.total_seconds():.2f}s")

            self.nx_workflow.cleanup_work_dir()

        except Exception as e:
            self.logger.error(f"Error during database creation: {e}")
            raise RuntimeError(f"Error during database creation: {e}")

