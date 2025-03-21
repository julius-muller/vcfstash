
import subprocess
from pathlib import Path
from typing import Dict, List
import yaml
import hashlib
from src.utils.logging import setup_logger


class WorkflowManager:
    """Handles workflow operations and configuration"""
    def __init__(self, project_root: Path):
        self.project_root = Path(project_root)
        log_file = self.project_root / "workflow.log"
        self.logger = setup_logger(log_file)
        self.logger.info("Initializing workflow manager")
        self.logger.debug(f"Project root: {project_root}")

    def load_nextflow_config(self, test_mode: bool = False) -> Dict:
        """Load Nextflow configuration from YAML files"""
        config_file = 'nextflow_test.yml' if test_mode else 'nextflow.yml'
        config_path = self.project_root / 'config' / config_file
        self.logger.debug(f"Loading config from: {config_path}")

        try:
            with open(config_path) as f:
                config = yaml.safe_load(f)
            self.logger.debug("Configuration loaded successfully")
            return config

        except FileNotFoundError:
            self.logger.error(f"Configuration file not found: {config_path}")
            raise
        except yaml.YAMLError as e:
            self.logger.error(f"Error parsing configuration: {e}")
            raise

    def get_workflow_hash(self, workflow_dir: Path) -> str:
        """Get combined hash of workflow files"""
        workflow_files = ['main.nf', 'nextflow.config']
        combined_content = b''

        self.logger.debug(f"Computing workflow hash for: {workflow_dir}")
        for file in workflow_files:
            path = Path(workflow_dir) / file
            if not path.exists():
                self.logger.error(f"Required workflow file not found: {path}")
                raise FileNotFoundError(f"Required workflow file not found: {path}")
            combined_content += path.read_bytes()

        workflow_hash = hashlib.md5(combined_content).hexdigest()
        self.logger.debug(f"Workflow hash: {workflow_hash}")
        return workflow_hash

    def store_workflow_dag(self, run_dir: Path, cmd: List[str]) -> None:
        """Generate and store workflow DAG visualization"""
        self.logger.info("Generating workflow DAG")
        self.logger.debug(f"Run directory: {run_dir}")

        try:
            try:
                del cmd[cmd.index("-with-trace")]
            except ValueError:
                pass

            dag_cmd = cmd + ["-preview", "-with-dag", str(run_dir / "flowchart.html")]
            self.logger.debug(f"DAG command: {' '.join(dag_cmd)}")

            subprocess.run(dag_cmd, check=True, capture_output=True)
            self.logger.info(f"Workflow DAG saved to: {run_dir}/flowchart.html")

        except subprocess.CalledProcessError as e:
            self.logger.warning(f"Failed to generate workflow DAG: {e}")
            self.logger.debug(f"Command output: {e.output if hasattr(e, 'output') else 'No output available'}")