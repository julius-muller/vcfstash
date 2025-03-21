from pathlib import Path
import subprocess
from datetime import datetime
from .base import VEPDatabase
from ..utils.validation import compute_md5


class DatabaseInitializer(VEPDatabase):
    """Handles database initialization"""
    def __init__(self, name: str, input_file: Path, fasta_ref: Path, output_dir: Path, threads: int):
        self.output_dir = Path(output_dir)
        super().__init__(self.output_dir / name)
        self.input_file = Path(input_file)
        self.fasta_ref = Path(fasta_ref)
        self.threads = max(int(threads), 1)
        self.logger.info(f"Initializing database: {name}")
        self.logger.debug(f"Input file: {input_file}")
        self.logger.debug(f"Reference: {fasta_ref}")
        self.logger.debug(f"Output directory: {output_dir}")
        self.logger.debug(f"Threads: {threads}")

    def initialize(self) -> None:
        """Initialize new VEP database"""
        if not self.input_file.exists():
            self.logger.error(f"Input BCF file does not exist: {self.input_file}")
            raise FileNotFoundError("Input BCF file does not exist.")

        if self.variants_bcf.exists():
            self.logger.error(f"Output database already exists: {self.variants_bcf}")
            raise FileExistsError("Output database already exists.")

        self.logger.info("Starting database initialization")
        self.ensure_indexed(self.input_file)
        self._validate_inputs()
        self._create_database()
        self.logger.info("Database initialization completed successfully")

    def _validate_inputs(self) -> None:
        """Validate the input BCF file format"""
        self.logger.debug("Validating input BCF format")
        is_valid, error = self.validate_bcf_header(self.input_file, norm=False)
        if not is_valid:
            self.logger.error(f"Invalid input BCF: {error}")
            raise ValueError(f"Invalid input BCF: {error}")
        self.logger.debug("Input BCF validation successful")

    def _create_database(self) -> None:
        """Create and initialize the database"""
        self.logger.info("Creating database structure")
        self.blueprint_dir.mkdir(parents=True, exist_ok=True)
        input_md5 = compute_md5(self.input_file)
        self.logger.debug(f"Input file MD5: {input_md5}")

        try:
            cmd = [
                "bcftools", "view", "-Ou", "--threads", str(self.threads), str(self.input_file),
                "|", "bcftools", "annotate", "-x", "INFO", "--threads", str(self.threads),
                "|", "bcftools", "norm", "-m-", "-f", str(self.fasta_ref), "-c", "x",
                "--threads", str(self.threads), "--rm-dup", "all",
                "-Ob", "--write-index", "-o", str(self.variants_bcf)
            ]

            self.logger.debug(f"Executing command: {' '.join(cmd)}")
            start_time = datetime.now()
            subprocess.run(" ".join(cmd), shell=True, check=True)
            duration = datetime.now() - start_time
            self.logger.info(f"Database creation completed in {duration.total_seconds():.2f} seconds")

            # Store initialization details
            self.info_file.parent.mkdir(parents=True, exist_ok=True)
            self.info_file.write_text("")  # Create empty info file
            self.logger.debug(f"Created info file: {self.info_file}")

            # Log initialization details
            self.logger.info("Database initialization summary:")
            self.logger.info(f"- Database name: {self.db_path.name}")
            self.logger.info(f"- Input file: {self.input_file}")
            self.logger.info(f"- Input MD5: {input_md5}")
            self.logger.info(f"- Processing time: {duration.total_seconds():.2f}s")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error during bcftools operation: {e}")
            self.logger.error(f"Command output: {e.output if hasattr(e, 'output') else 'No output available'}")
            raise RuntimeError(f"Error during bcftools operation: {e}")