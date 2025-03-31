##########################################################################
#                                                                        #
# This file contains abstract base classes and implementations           #
# for managing output file structures in VEP-based workflows.
# NOT CURRENTLY USED.                                            #
#                                                                        #
##########################################################################

# vepstash/outputs.py
from abc import ABC, abstractmethod
from pathlib import Path


class BaseOutput(ABC):
    """
    Partially abstract base for output structures. 
    Subclasses define the exact dirs/files to create or check.
    """

    def __init__(self, root_dir: str):
        self.root_dir = Path(root_dir)

    @abstractmethod
    def create_structure(self) -> None:
        """Create needed subdirs and files as placeholders."""
        pass

    @abstractmethod
    def validate_structure(self) -> bool:
        """Ensure required dirs/files exist; return True/False."""
        pass


class StashOutput(BaseOutput):
    """
    Encapsulates the structure for stash-init / stash-add:

      <root_dir>/
      ├── blueprint/
      ├── vepdb.log
      └── workflow/
          ├── main.nf
          ├── modules/
          │   ├── annotate.nf
          │   ├── intersect.nf
          │   ├── merge.nf
          │   ├── merge_variants.nf
          │   ├── normalize.nf
          │   └── utils.nf
          ├── add_ddb497f91944bd880b65655594538a2f_nextflow.config (example)
          ├── annotation.config
          ├── init_nextflow.config
          └── workflow.log
    """

    def create_structure(self) -> None:
        # Top-level elements
        (self.root_dir / "blueprint").mkdir(parents=True, exist_ok=True)

        # workflow directory & sub-structure
        workflow_dir = self.root_dir / "workflow"
        workflow_dir.mkdir(parents=True, exist_ok=True)

        # modules
        modules_dir = workflow_dir / "modules"
        modules_dir.mkdir(exist_ok=True)


    def validate_structure(self) -> bool:
        # Minimal existence checks
        required_paths = [
            self.root_dir / "blueprint",
            self.root_dir / "vepdb.log",
            self.root_dir / "workflow" / "main.nf",
            self.root_dir / "workflow" / "modules" / "annotate.nf",
            self.root_dir / "workflow" / "modules" / "intersect.nf",
            self.root_dir / "workflow" / "modules" / "merge.nf",
            self.root_dir / "workflow" / "modules" / "merge_variants.nf",
            self.root_dir / "workflow" / "modules" / "normalize.nf",
            self.root_dir / "workflow" / "modules" / "utils.nf",
            self.root_dir / "workflow" / "annotation.config",
            self.root_dir / "workflow" / "init_nextflow.config",
            self.root_dir / "workflow" / "workflow.log"
        ]
        return all(p.exists() for p in required_paths)


class AnnotateOutput(BaseOutput):
    """
    Encapsulates the structure for annotation workflows. Example:

      <root_dir>/
      ├── annotations/
      │   └── <any subfolders, e.g. 'testor'>
      └── workflow/
          ├── main.nf
          ├── modules/
          │   ├── annotate.nf
          │   ├── intersect.nf
          │   ├── merge.nf
          │   ├── merge_variants.nf
          │   ├── normalize.nf
          │   └── utils.nf
          ├── testor_annotation.config
          └── testor_nextflow.config
    """

    def create_structure(self) -> None:
        # Create base annotations dir
        annotations_dir = self.root_dir / "annotations"
        annotations_dir.mkdir(parents=True, exist_ok=True)

        # workflow directory & sub-structure
        workflow_dir = self.root_dir / "workflow"
        workflow_dir.mkdir(parents=True, exist_ok=True)

        # modules
        modules_dir = workflow_dir / "modules"
        modules_dir.mkdir(exist_ok=True)


    def validate_structure(self) -> bool:
        required_paths = [
            self.root_dir / "annotations",
            self.root_dir / "workflow" / "main.nf",
            self.root_dir / "workflow" / "modules" / "annotate.nf",
            self.root_dir / "workflow" / "modules" / "intersect.nf",
            self.root_dir / "workflow" / "modules" / "merge.nf",
            self.root_dir / "workflow" / "modules" / "merge_variants.nf",
            self.root_dir / "workflow" / "modules" / "normalize.nf",
            self.root_dir / "workflow" / "modules" / "utils.nf",
            self.root_dir / "workflow" / "*_annotation.config",
            self.root_dir / "workflow" / "*_nextflow.config",
        ]
        return all(p.exists() for p in required_paths)
