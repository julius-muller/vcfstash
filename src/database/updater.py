from pathlib import Path
import subprocess
import os
from datetime import datetime
from .base import VEPDatabase
from ..utils.validation import compute_md5, get_bcf_stats, check_duplicate_md5

class DatabaseUpdater(VEPDatabase):
    """Handles adding new variants to database"""
    def __init__(self, db_path: Path, input_file: Path, fasta_ref: Path, threads: int):
        super().__init__(db_path)
        self.input_file = Path(input_file)
        self.fasta_ref = Path(fasta_ref)
        self.threads = max(int(threads), 1)
        self.logger.info("Initializing database update")
        self.logger.debug(f"Input file: {input_file}")
        self.logger.debug(f"Reference: {fasta_ref}")
        self.logger.debug(f"Threads: {threads}")

    def add(self) -> None:
        """Add new variants to existing database"""
        self.logger.info("Starting database update")
        self._validate_inputs()
        self._merge_variants()
        self.logger.info("Database update completed successfully")

    def _validate_inputs(self) -> None:
        """Validate input files and check for duplicates"""
        self.logger.debug("Validating inputs")

        if not self.variants_bcf.exists():
            self.logger.error("Database BCF file does not exist")
            raise FileNotFoundError("Database BCF file does not exist")

        if not self.input_file.exists():
            self.logger.error("Input VCF/BCF file does not exist")
            raise FileNotFoundError("Input VCF/BCF file does not exist")

        self.ensure_indexed(self.variants_bcf)
        self.ensure_indexed(self.input_file)

        input_md5 = compute_md5(self.input_file)
        self.logger.debug(f"Input file MD5: {input_md5}")

        if check_duplicate_md5(self.info_file, input_md5):
            self.logger.error("Duplicate file detected (MD5 match)")
            raise ValueError("This file was already added to the database (MD5 match)")

        self.logger.debug("Input validation successful")

    def _merge_variants(self) -> None:
        """Merge new variants into the database"""
        temp_merged = Path(str(self.variants_bcf) + ".tmp.bcf")
        self.logger.info("Starting variant merge")
        self.logger.debug(f"Using temporary file: {temp_merged}")

        try:
            start_time = datetime.now()
            cmd = [
                "bcftools", "concat",
                "--allow-overlaps",
                "--rm-dup", "all",
                "-Ob",
                "--write-index",
                "-o", str(temp_merged),
                "--threads", str(self.threads),
                str(self.variants_bcf),
                str(self.input_file)
            ]
            self.logger.debug(f"Executing command: {' '.join(str(x) for x in cmd)}")

            subprocess.run(cmd, check=True)
            duration = datetime.now() - start_time

            stats = get_bcf_stats(temp_merged)
            self.logger.debug("Merge completed, updating database files")

            # Replace old database with merged file
            os.replace(temp_merged, self.variants_bcf)
            os.replace(str(temp_merged) + ".csi", str(self.variants_bcf) + ".csi")

            # Log update details
            input_md5 = compute_md5(self.input_file)
            self.logger.info("Database update summary:")
            self.logger.info(f"- Added file: {self.input_file}")
            self.logger.info(f"- Input MD5: {input_md5}")
            self.logger.info("- Database statistics after update:")
            for key, value in stats.items():
                self.logger.info(f"  {key}: {value}")
            self.logger.info(f"- Processing time: {duration.total_seconds():.2f}s")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to merge variants: {e}")
            self.logger.debug("Cleaning up temporary files")
            temp_merged.unlink(missing_ok=True)
            Path(str(temp_merged) + ".csi").unlink(missing_ok=True)
            raise RuntimeError(f"Failed to add variants: {e}")