# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional

from vcfcache.database.base import VCFDatabase
from vcfcache.utils.validation import check_duplicate_md5, compute_md5, get_bcf_stats


class DatabaseUpdater(VCFDatabase):
    """Represents a specialized utility for updating a VCF/BCF database.

    This class extends the functionality of the VCFDatabase to facilitate the addition
    of new variants into the existing database. It validates input files, manages input
    configuration and parameters, and implements workflows to integrate new data into
    the database. The primary purpose is to ensure VCF/BCF files are merged seamlessly
    into the database while maintaining consistency and log tracking. It also keeps
    blueprints around input files and database statistics up-to-date. Validation of input
    VCF reference files against YAML configuration is supported as well.

    Attributes:
        db_path (Path): Path to the database directory.
        input_file (Path): File path to the input VCF/BCF file.
        params_file (Optional[Path]): Optional parameters file path (auto-generated if not provided).
        verbosity (int): Verbosity level for logging. Higher numbers increase detail.
        debug (bool): Debug flag indicating whether debugging is enabled.
        threads (int): Number of threads for bcftools operations (default: 1).
    """

    def __init__(
        self,
        db_path: Path | str,
        input_file: Path | str,
        bcftools_path: Path | str,
        params_file: Optional[Path] | Optional[str] = None,
        verbosity: int = 0,
        debug: bool = False,
        normalize: bool = False,
    ):
        super().__init__(Path(db_path), verbosity, debug, bcftools_path)
        self.cached_output.validate_structure()
        self.logger = self.connect_loggers()
        self.input_file = Path(input_file).expanduser().resolve()
        self.input_md5 = compute_md5(self.input_file)  # might take too long to do here
        self.normalize = normalize

        # Handle params file - use existing or auto-generate
        if params_file:
            self.params_file = self.blueprint_dir / f"add_{self.input_md5}.yaml"
            params_path = (
                Path(params_file) if isinstance(params_file, str) else params_file
            )
            shutil.copyfile(params_path.expanduser().resolve(), self.params_file)
        else:
            # Check if init.yaml exists from blueprint-init
            wfini = self.workflow_dir / "init.yaml"
            if wfini.exists():
                self.params_file = wfini
            else:
                # Auto-generate minimal params file with all required fields
                import yaml
                self.params_file = self.workflow_dir / "extend.yaml"
                minimal_params = {
                    "annotation_tool_cmd": str(bcftools_path),  # Use bcftools as default
                    "bcftools_cmd": str(bcftools_path),
                    "temp_dir": "/tmp",
                    "threads": 1,  # Default to 1 thread
                    "genome_build": "UNKNOWN",
                    "optional_checks": {},
                }
                self.params_file.write_text(yaml.dump(minimal_params))

        # Initialize workflow backend (pure Python)
        from vcfcache.database.base import create_workflow
        self.nx_workflow = create_workflow(
            input_file=self.input_file,
            output_dir=self.blueprint_dir,
            name=f"add_{self.input_md5}",
            params_file=self.params_file,
            verbosity=self.verbosity,
        )

        self._validate_inputs()
        # Log initialization parameters
        if self.logger:
            self.logger.info("Initializing database update")
            self.logger.debug(f"Input file: {input_file}")

    def _validate_inputs(self) -> None:
        """Validate input files, directories, and YAML parameters"""
        if self.logger:
            self.logger.debug("Validating inputs")

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

    def add(self) -> None:
        """Add new variants to existing database

        self = DatabaseUpdater(db_path=Path('~/projects/vcfcache/tests/data/test_out/nftest'),
        input_file=Path('~/projects/vcfcache/tests/data/nodata/dbsnp_test.bcf'),
        verbosity=2)
        profile='test'
        """
        if self.logger:
            self.logger.info("Starting database update")

        try:
            # Check for duplicate before validation
            if check_duplicate_md5(db_info=self.db_info, new_md5=self.input_md5):
                if self.logger:
                    self.logger.warning(
                        f"Skipping duplicate file (MD5 match): {self.input_file}"
                    )
                return

            self._validate_input_files()
            self._merge_variants()
            if self.logger:
                self.logger.info("Database update completed successfully")

        except FileNotFoundError as e:
            if self.logger:
                self.logger.error(str(e))
            raise
        except ValueError as e:
            if self.logger:
                self.logger.error(str(e))
            raise

    def _validate_input_files(self) -> None:
        """Validate input files and check for duplicates"""
        if self.logger:
            self.logger.debug("Validating inputs")

        # First check database BCF
        if not self.blueprint_bcf.exists():
            if self.logger:
                self.logger.error("Database BCF file does not exist")
            raise FileNotFoundError("Database BCF file does not exist")

        # Check input VCF/BCF
        if not self.input_file.exists():
            if self.logger:
                self.logger.error("Input VCF/BCF file does not exist")
            raise FileNotFoundError("Input VCF/BCF file does not exist")

        # Validate remaining inputs
        self.ensure_indexed(self.blueprint_bcf)
        self.ensure_indexed(self.input_file)
        if self.logger:
            self.logger.debug("Input validation successful")

    def _merge_variants(self) -> None:
        """Merge new variants into the database"""
        if self.logger:
            self.logger.info("Starting variant merge")

        try:
            # Get statistics before merge
            pre_stats = get_bcf_stats(self.blueprint_bcf, bcftools_path=self.bcftools_path)

            # Run the workflow in database mode
            start_time = datetime.now()
            # Pass normalize flag via nextflow_args
            nextflow_args = ["--normalize"] if getattr(self, "normalize", False) else None
            self.nx_workflow.run(
                db_mode="blueprint-extend",
                trace=True,
                dag=True,
                report=True,
                db_bcf=self.blueprint_bcf,
                nextflow_args=nextflow_args,
            )
            duration = datetime.now() - start_time

            # Get statistics after merge
            post_stats = get_bcf_stats(self.blueprint_bcf, bcftools_path=self.bcftools_path)

            # Calculate differences
            diff_stats = {}
            for key in set(pre_stats.keys()) | set(post_stats.keys()):
                try:
                    pre_val = int(pre_stats.get(key, 0))
                    post_val = int(post_stats.get(key, 0))
                    diff_stats[key] = post_val - pre_val
                except ValueError:
                    continue

            if self.logger:
                self.logger.debug("Merge completed, updating database files")

            self.db_info["input_files"].append(
                {
                    "path": str(self.input_file),
                    "md5": self.input_md5,
                    "added": datetime.now().isoformat(),
                }
            )

            with open(self.info_file, "w") as f:
                json.dump(self.db_info, f, indent=2)

            # Log update details
            if self.logger:
                self.logger.info("Database update summary:")
                self.logger.info(f"- Added file: {self.input_file}")
                self.logger.info(f"- Input MD5: {self.input_md5}")
                self.logger.info("- Database statistics changes:")
                for key, diff in diff_stats.items():
                    prefix = "+" if diff > 0 else ""
                    self.logger.info(
                        f"  {key}: {prefix}{diff:d} ({pre_stats[key]} -> {post_stats[key]})"
                    )
                self.logger.info(f"- Processing time: {duration.total_seconds():.2f}s")

            # Write completion flag
            from vcfcache.utils.completion import write_completion_flag
            write_completion_flag(
                output_dir=self.cached_output.root_dir,
                command="blueprint-extend",
                mode="blueprint-extend"
            )

            if not self.debug:
                self.nx_workflow.cleanup_work_dir()

        except subprocess.CalledProcessError as e:
            if self.logger:
                self.logger.error(f"Failed to merge variants: {e}")
            raise RuntimeError(f"Failed to add variants: {e}") from e
