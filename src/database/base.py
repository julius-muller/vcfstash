import subprocess
import sys
from pathlib import Path
from typing import Tuple, List, Optional

import yaml
import hashlib
from src.utils.validation import ensure_indexed, validate_bcf_header
from src.utils.logging import setup_logger

class VEPDatabase:
    """Base class for VEP database operations"""
    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        self.blueprint_dir = self.db_path / "blueprint"
        self.variants_bcf = self.blueprint_dir / "variants.bcf"
        self.info_file = self.blueprint_dir / "variants.info"

        # Set up logging
        log_file = self.db_path / "vepdb.log"
        self.logger = setup_logger(log_file)
        self.logger.info(f"Initializing VEP database: {db_path}")

    def ensure_indexed(self, file_path: Path) -> None:
        """Validate that BCF/VCF file has an index"""
        self.logger.debug(f"Checking index for: {file_path}")
        ensure_indexed(file_path)
        self.lo

    def validate_bcf_header(self, bcf_path: Path, norm: bool = True) -> Tuple[bool, str]:
        """Validate BCF header format"""
        self.logger.debug(f"Validating BCF header: {bcf_path}")
        result = validate_bcf_header(bcf_path, norm)
        if not result[0]:
            self.logger.error(f"BCF header validation failed: {result[1]}")
        return result

class NextflowWorkflow:
    """Base class for Nextflow workflow operations"""
    def __init__(self, workflow_dir: Path, filename: str = "main", params_file: Optional[Path] = None):
        # Set up logging
        log_file = Path(workflow_dir) / "workflow.log"
        self.logger = setup_logger(log_file)
        self.logger.info(f"Initializing Nextflow workflow in: {workflow_dir}")

        self.workflow_dir = Path(workflow_dir)
        self.workflow_path = self.workflow_dir / f"{filename}.nf"
        if not self.workflow_path.exists():
            self.logger.error(f"Workflow file not found: {self.workflow_path}")
            raise FileNotFoundError(f"Workflow file not found: {self.workflow_path}")

        self.workflow_config_path = self.workflow_dir / "nextflow.config"
        if not self.workflow_config_path.exists():
            self.logger.error(f"Nextflow config file not found: {self.workflow_config_path}")
            raise FileNotFoundError(f"Nextflow config file not found: {self.workflow_config_path}")

        # Handle optional params_file
        self.params_file = Path(params_file) if params_file else None
        if self.params_file and not self.params_file.exists():
            self.logger.error(f"Parameters file not found: {self.params_file}")
            raise FileNotFoundError(f"Parameters file not found: {self.params_file}")


        self.workflow_hash = self.get_workflow_hash(self.workflow_dir)
        self.config = self.load_nextflow_config()
        self.run_dir = self.workflow_dir / "run"
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.logger.debug(f"Workflow hash: {self.workflow_hash}")

    def load_nextflow_config(self, test_mode=False):
        """Load Nextflow configuration from YAML files."""
        if self.params_file:
            self.logger.debug(f"Loading params from: {self.params_file}")
            with open(self.params_file) as f:
                config = yaml.safe_load(f)
            self.logger.debug("Parameters loaded successfully")
            return config
        return {}

    def get_workflow_hash(self, workflow_dir: Path) -> str:
        """Get combined hash of workflow files"""
        workflow_files = ['main.nf', 'nextflow.config']
        combined_content = b''

        self.logger.debug("Computing workflow hash")
        for file in workflow_files:
            path = Path(workflow_dir) / file
            if not path.exists():
                self.logger.error(f"Required workflow file not found: {path}")
                raise FileNotFoundError(f"Required workflow file not found: {path}")
            combined_content += path.read_bytes()

        workflow_hash = hashlib.md5(combined_content).hexdigest()
        self.logger.debug(f"Workflow hash computed: {workflow_hash}")
        return workflow_hash

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


    def run(self, input_file: Path, output_dir: Path, db_mode: bool = False, **kwargs) -> subprocess.CompletedProcess:
        """Run Nextflow pipeline with the appropriate configuration."""
        if not self.workflow_path.is_file():
            self.logger.error(f"Required workflow file not found: {self.workflow_path}")
            raise FileNotFoundError(f"Required workflow file not found: {self.workflow_path}")

        cmd = [
            "nextflow", "run", str(self.workflow_path),
            "-with-trace",
            "--input", str(input_file),
            "--output", str(output_dir)
        ]

        if db_mode:
            cmd.extend(["--db_mode", "true"])

        if self.params_file:
            cmd.extend(["-params-file", str(self.params_file)])

        for key, value in kwargs.items():
            if key == "nextflow_args" and value:
                args = [arg for arg in value if arg and arg != '[]']
                if args:
                    cmd.extend(args)
            elif value:
                cmd.extend([f"--{key}", str(value)])

        self.logger.info("Starting Nextflow workflow")
        self.logger.debug(f"Command: {' '.join(str(x) for x in cmd)}")

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.info("Workflow completed successfully")
            self.logger.debug(result.stdout)
            self.store_workflow_dag(output_dir, cmd)
            return result

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Workflow execution failed (exit code: {e.returncode})")
            self.logger.error(f"STDOUT: {e.stdout}")
            self.logger.error(f"STDERR: {e.stderr}")
            raise

