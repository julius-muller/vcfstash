import json
import shutil
from pathlib import Path
import subprocess
from datetime import datetime
from vcfstash.database.base import VCFDatabase, NextflowWorkflow
from vcfstash.utils.validation import compute_md5, get_bcf_stats, check_duplicate_md5

class DatabaseUpdater(VCFDatabase):
    """Handles adding new variants to database"""
    def __init__(self, db_path: Path | str, input_file: Path | str, config_file: Path | str = None,
                 params_file: Path | str = None, verbosity: int = 0, debug: bool = False):
        super().__init__(Path(db_path), verbosity, debug)
        self.stashed_output.validate_structure()
        self.logger = self.connect_loggers()
        self.input_file = Path(input_file).expanduser().resolve()
        self.input_md5 = compute_md5(self.input_file) # might take too long to do here
        if config_file:
            self.config_file = self.blueprint_dir / f'add_{self.input_md5}.config'
            shutil.copyfile(config_file.expanduser().resolve(), self.config_file)
        else:
            wfini = self.workflow_dir / "init.config"
            self.config_file = wfini if wfini.exists() else None

        if params_file:
            self.params_file = self.blueprint_dir / f'add_{self.input_md5}.yaml'
            shutil.copyfile(params_file.expanduser().resolve(), self.params_file)
        else:
            wfini = self.workflow_dir / "init.yaml"
            self.params_file = wfini if wfini.exists() else None

        # Initialize NextflowWorkflow
        self.nx_workflow = NextflowWorkflow(
            input_file=self.input_file,
            output_dir=self.blueprint_dir,
            name=f'add_{self.input_md5}',
            workflow=self.workflow_dir / "main.nf",
            config_file=self.config_file,
            params_file=self.params_file,
            verbosity=self.verbosity
        )

        # Log initialization parameters
        self.logger.info("Initializing database update")
        self.logger.debug(f"Input file: {input_file}")

    def add(self) -> None:
        """
        Add new variants to existing database

        self = DatabaseUpdater(db_path=Path('~/projects/vcfstash/tests/data/test_out/nftest'),
        input_file=Path('~/projects/vcfstash/tests/data/nodata/dbsnp_test.bcf'),
        verbosity=2)
        profile='test'
        """
        self.logger.info("Starting database update")

        try:
            # Check for duplicate before validation
            if check_duplicate_md5(db_info=self.db_info, new_md5=self.input_md5):
                self.logger.warning(f"Skipping duplicate file (MD5 match): {self.input_file}")
                return

            self._validate_inputs()
            self._merge_variants()
            self.logger.info("Database update completed successfully")

        except FileNotFoundError as e:
            self.logger.error(str(e))
            raise
        except ValueError as e:
            self.logger.error(str(e))
            raise

    def _validate_inputs(self) -> None:
        """Validate input files and check for duplicates"""
        self.logger.debug("Validating inputs")

        # First check database BCF
        if not self.blueprint_bcf.exists():
            self.logger.error("Database BCF file does not exist")
            raise FileNotFoundError("Database BCF file does not exist")

        # Check input VCF/BCF
        if not self.input_file.exists():
            self.logger.error("Input VCF/BCF file does not exist")
            raise FileNotFoundError("Input VCF/BCF file does not exist")

        # Validate remaining inputs
        self.ensure_indexed(self.blueprint_bcf)
        self.ensure_indexed(self.input_file)
        self.logger.debug("Input validation successful")

    def _merge_variants(self) -> None:
        """Merge new variants into the database"""

        self.logger.info("Starting variant merge")

        try:
            # Get statistics before merge
            pre_stats = get_bcf_stats(self.blueprint_bcf)

            # Run the workflow in database mode
            start_time = datetime.now()
            self.nx_workflow.run(
                db_mode="stash-add",
                trace=True,
                dag=True,
                report=True,
                db_bcf=self.blueprint_bcf
            )
            duration = datetime.now() - start_time

            # Get statistics after merge
            post_stats = get_bcf_stats(self.blueprint_bcf)

            # Calculate differences
            diff_stats = {}
            for key in set(pre_stats.keys()) | set(post_stats.keys()):
                try:
                    pre_val = int(pre_stats.get(key, 0))
                    post_val = int(post_stats.get(key, 0))
                    diff_stats[key] = post_val - pre_val
                except ValueError:
                    continue

            self.logger.debug("Merge completed, updating database files")

            self.db_info["input_files"].append({
                "path": str(self.input_file),
                "md5": self.input_md5,
                "added": datetime.now().isoformat()
            })

            with open(self.info_file, "w") as f:
                json.dump(self.db_info, f, indent=2)

            # Log update details
            self.logger.info("Database update summary:")
            self.logger.info(f"- Added file: {self.input_file}")
            self.logger.info(f"- Input MD5: {self.input_md5}")
            self.logger.info("- Database statistics changes:")
            for key, diff in diff_stats.items():
                prefix = "+" if diff > 0 else ""
                self.logger.info(f"  {key}: {prefix}{diff:d} ({pre_stats[key]} -> {post_stats[key]})")
            self.logger.info(f"- Processing time: {duration.total_seconds():.2f}s")
            if not self.debug:
                self.nx_workflow.cleanup_work_dir()

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to merge variants: {e}")
            raise RuntimeError(f"Failed to add variants: {e}")