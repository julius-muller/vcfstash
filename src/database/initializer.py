import json
import shutil
from pathlib import Path
from datetime import datetime
from src.database.base import VEPDatabase, NextflowWorkflow
from src.utils.validation import compute_md5


class DatabaseInitializer(VEPDatabase):
    """Handles database initialization"""
    def __init__(self, name: str, input_file: Path | str, config_file: Path | str = None, output_dir: Path | str = Path("."),
                 verbosity: int = 0, force: bool = False, test_mode: bool = False) -> None:
        """Initialize the database creator.

        Args:
            name: Name of the new database
            input_file: Path to input BCF/VCF file (required)
            output_dir: Output directory (default: current directory)
            verbosity: Logging verbosity level (0=WARNING, 1=INFO, 2=DEBUG)        """
        # Create the database directory first
        self.input_file = Path(input_file).expanduser().resolve()
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.stash_path = self.output_dir / name
        self.test_mode = test_mode

        # Remove destination directory if it exists to ensure clean copy
        if force and self.stash_path.exists():
            print(f"Stash directory already exists, removing: {self.stash_path}")
            shutil.rmtree(self.stash_path)

        # Now initialize the parent class
        super().__init__(self.stash_path, test_mode, verbosity)
        self.validate_labels(name)

        # Create the database directory
        self._setup_cachedir(force = force)

        self._copy_workflow_srcfiles(destination=self.workflow_dir, skip_config=True)

        self.config_file = self.setup_config(config_file=config_file, config_name='init_nextflow.config')

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
        self.logger.debug(f"Output directory: {self.output_dir}")
        self.logger.debug(f"Config file: {self.config_file}")


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

        # Check output directory permissions
        if not self.output_dir.exists():
            try:
                self.output_dir.mkdir(parents=True)
            except PermissionError:
                msg = f"Cannot create output directory: {self.output_dir}"
                self.logger.error(msg)
                raise PermissionError(msg)

        self.logger.debug("Input validation successful")

    def _create_database(self) -> None:
        """Create and initialize the database"""
        self.logger.info("Creating database from normalized and annotated variants...")

        self.blueprint_dir.mkdir(parents=True, exist_ok=True)

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
            self.logger.info(f"- Output file: {self.output_dir / self.blueprint_bcf}")
            self.logger.info(f"- Input MD5: {input_md5}")
            self.logger.info(f"- Processing time: {duration.total_seconds():.2f}s")

            self.nx_workflow.cleanup_work_dir()

        except Exception as e:
            self.logger.error(f"Error during database creation: {e}")
            raise RuntimeError(f"Error during database creation: {e}")

