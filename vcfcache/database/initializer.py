import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

from vcfcache.database.base import VCFDatabase
from vcfcache.utils.validation import compute_md5


class DatabaseInitializer(VCFDatabase):
    """Class for initializing and managing a VCF database.

    The `DatabaseInitializer` class is designed to set up and initialize a VCF database by
    processing input VCF/BCF files, managing configurations, and running associated workflows.
    It provides tools for validating inputs, ensuring output directory structure, and creating
    the database while handling workflow execution and logging.

    """

    def __init__(
        self,
        input_file: Path | str,
        bcftools_path: Path | str | None,
        params_file: Optional[Path | str] = None,
        output_dir: Path | str = Path("."),
        verbosity: int = 0,
        force: bool = False,
        debug: bool = False,
        threads: int = 1,
        normalize: bool = False,
    ) -> None:
        """Initialize the database creator.

        Args:
            input_file: Path to input BCF/VCF file (required)
            bcftools_path: Path to bcftools binary
            params_file: Path to params YAML file (optional, will auto-generate if not provided)
            output_dir: Output directory (default: current directory)
            verbosity: Logging verbosity level (0=WARNING, 1=INFO, 2=DEBUG)
            force: Force overwrite of existing database (default: False)
            debug: Boolean to enable or disable debug mode (default: False)
            threads: Number of threads for bcftools (default: 1)
            normalize: If True, split multiallelic variants using bcftools norm -m-

        """
        # Initialize the parent class
        super().__init__(
            Path(output_dir) if isinstance(output_dir, str) else output_dir,
            verbosity,
            debug,
            bcftools_path,
        )
        self._setup_cache(force=force)
        self.logger = self.connect_loggers()

        # self.validate_label(name)
        self.input_file = Path(input_file).expanduser().resolve()

        self._copy_workflow_srcfiles(
            source=self.workflow_dir_src,
            destination=self.workflow_dir,
            skip_config=True,
        )
        self.normalize = normalize

        # Handle params file - create minimal one if not provided
        self.config_yaml = self.workflow_dir / "init.yaml"
        if params_file:
            shutil.copyfile(Path(params_file).expanduser().resolve(), self.config_yaml)
            params_path = Path(params_file) if isinstance(params_file, str) else params_file
        else:
            # Create minimal params file with just bcftools and threads
            import yaml
            minimal_params = {
                "bcftools_cmd": str(bcftools_path),
                "threads": threads,
            }
            self.config_yaml.write_text(yaml.dump(minimal_params))
            params_path = self.config_yaml

        # Initialize workflow backend (pure Python)
        if self.logger:
            self.logger.info("Initializing pure Python workflow...")

        from vcfcache.database.base import create_workflow

        self.nx_workflow = create_workflow(
            input_file=self.input_file,
            output_dir=self.blueprint_dir,
            name="init",
            params_file=params_path,
            verbosity=self.verbosity,
        )

        self._validate_inputs()

        # Log initialization parameters
        if self.logger:
            self.logger.info(f"Initializing database: {self.cache_name}")
            self.logger.debug(f"Input file: {self.input_file}")
            self.logger.debug(f"Output directory: {self.blueprint_dir}")

    def _log_contigs(self) -> None:
        """Log a preview of contigs present in the input BCF/VCF (top 30)."""

        try:
            result = subprocess.run(
                [str(self.bcftools_path), "index", "-s", str(self.input_file)],
                capture_output=True,
                text=True,
                check=True,
            )
            contigs = [line.split("\t", 1)[0] for line in result.stdout.splitlines()]

            # Sort first by string length, then lexicographically
            contigs.sort(key=lambda c: (len(c), c))

            preview = ", ".join(contigs[:30]) if contigs else "(none)"
            msg = f"Cache will be set up including the following contigs (top 30): {preview}"
            if self.logger:
                self.logger.info(msg)
            else:
                print(msg)
        except Exception as exc:
            if self.logger:
                self.logger.warning(
                    f"Could not list contigs via bcftools index -s: {exc}"
                )
            else:
                print(f"Warning: Could not list contigs via bcftools index -s: {exc}")

    def _setup_cache(self, force: bool) -> None:
        # Remove destination directory if it exists to ensure clean copy
        if self.cached_output.root_dir.exists():
            if (
                self.cached_output.validate_structure()
            ):  # we dont want to remove a random dir....
                if force:
                    print(
                        f"Cache directory already exists, removing: {self.cached_output.root_dir}"
                    )
                    shutil.rmtree(self.cached_output.root_dir)
                else:
                    raise FileExistsError(
                        f"Output directory already exists: {self.cached_output.root_dir}\nIf intended, use --force to overwrite."
                    )
            else:
                raise FileExistsError(
                    f"Output directory with an invalid structure detected: {self.cached_output.root_dir}"
                )

        print(f"Creating cache structure: {self.cached_output.root_dir}")
        self.cached_output.create_structure()

    def initialize(self) -> None:
        """Initialize new VCF database
        self = DatabaseInitializer(name='nftest', input_file=Path('tests/data/nodata/gnomad_test.bcf'),
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

        if self.logger:
            self.logger.debug("Input validation successful")

    def _create_database(self) -> None:
        """Create and initialize the database"""
        if self.logger:
            self.logger.info(
                "Creating database from normalized and annotated variants..."
            )

        self.cached_output.validate_structure()

        try:
            db_info: Dict[str, Any] = {
                "name": self.cache_name,
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
            # Log contig preview
            self._log_contigs()
            self.nx_workflow.run(
                db_mode="blueprint-init",
                trace=True,
                dag=True,
                report=True,
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
