##########################################################################
#                                                                        #
# This file contains abstract base classes and implementations           #
# for managing output file structures in VEP-based workflows.
# NOT CURRENTLY USED.                                            #
#                                                                        #
##########################################################################
import os
import warnings
from abc import ABC, abstractmethod
from pathlib import Path


class BaseOutput(ABC):
    """
    Partially abstract base for output structures. 
    Subclasses define the exact dirs/files to create or check.
    """

    def __init__(self, root_dir: str):
        self.root_dir = Path(root_dir).expanduser().resolve()
        # define the base directory of the module
        self.module_src_dir = Path(os.getenv('VEPSTASH_HOME', Path('.').resolve())) if '__file__' not in globals() else Path(__file__).parent.parent.parent

    @abstractmethod
    def required_paths(self) -> dir:
        """Dict with structure {'label':Path(),...} of required paths to check for existence."""
        pass

    @abstractmethod
    def create_structure(self) -> None:
        """Create needed subdirs and files as placeholders."""
        pass

    @abstractmethod
    def validate_structure(self) -> bool:
        """Ensure required dirs/files exist; return True/False."""
        pass

    @staticmethod
    def validate_label(label: str) -> None:
        """
        Validates that the label is valid in this context.
        """
        if len(label) > 30:
            raise ValueError(f"Annotation name must be less than 30 characters, but has {len(label)}: {label}")
        if " " in label:
            raise ValueError(f"Annotation name must not contain white spaces: {label}")
        if not all(c.isalnum() or c in "_-." for c in label):
            raise ValueError(f"Annotation name must only contain alphanumeric characters, underscores, dots, or dashes: {label}")


    @staticmethod
    def create_directories(dirs_to_create: dict) -> None:
        # Create directories and verify they exist
        for name, dir_path in dirs_to_create.items():
            if dir_path.is_dir:
                try:
                    # Create directory with parents if it doesn't exist
                    dir_path.mkdir(parents=True, exist_ok=True)

                    # Verify the directory exists after creation
                    if not dir_path.exists():
                        raise RuntimeError(f"Failed to create {name} directory: {dir_path}")

                    # Verify it's actually a directory
                    if not dir_path.is_dir():
                        raise RuntimeError(f"Path exists but is not a directory: {dir_path}")

                    # Verify we have write access by creating and removing a test file
                    test_file = dir_path / ".write_test"
                    try:
                        test_file.touch()
                        test_file.unlink()
                    except (IOError, PermissionError) as e:
                        raise RuntimeError(f"No write permission in {name} directory {dir_path}: {e}")

                except Exception as e:
                    # Catch any other exceptions that might occur during directory setup
                    raise RuntimeError(f"Error setting up {name} directory {dir_path}: {e}")


class StashOutput(BaseOutput):
    """
    Encapsulates the structure for stash-init / stash-add:

      <stash_root_dir>/
      ├── blueprint/
      ├── annotations/
      └── workflow/
          ├── ... parse from src
          ├── modules/
          │   ├── ... parse from src


    self = StashOutput(stash_root_dir='.')
    """

    def __init__(self, stash_root_dir: str):
        super().__init__(stash_root_dir)
        self.stash_root_dir = self.root_dir
        self.workflow_dir = self.stash_root_dir / "workflow"
        self.workflow_src_dir = self.module_src_dir / "workflow"

    def required_paths(self) -> dict:
        """
        Returns a dictionary with the required paths for the stash output structure.
        """

        return  {
            "blueprint": self.stash_root_dir / "blueprint",
            "annotations": self.stash_root_dir / "annotations",
            "workflow": self.workflow_dir,
            "workflow_src": self.module_src_dir / "workflow",
            "modules": self.workflow_dir / "modules",
        }
    
    def create_structure(self) -> None:
        req_dirs = {k:v for k,v in self.required_paths().items() if k != "workflow_src"}
        self.create_directories(req_dirs)

    def validate_structure(self) -> bool:

        # Minimal existence checks
        required_paths = self.required_paths()

        # for path in self.workflow_src_dir.rglob("*"):  # Recursively find all files and dirs
        #     if not path.name.endswith(".config"):  # Exclude .config files
        #         required_paths[f"{path.parent.stem}>{path.name}"] = self.workflow_dir / path.name

        for pname, path in required_paths.items():
            if not path.exists():
                warnings.warn(f"Missing required path {pname}: {path}")
                return False
        return True

class AnnotatedStashOutput(BaseOutput):
    """
    Encapsulates the structure for annotation stash from stash-annotate. Example:

      <stash_root_dir>/
      ├── annotations/
      │   └── <any subfolders, e.g. 'test'>  <- annotation_dir

    """

    def __init__(self, annotation_dir: str):
        super().__init__(annotation_dir)
        self.annotation_dir = self.root_dir
        self.annotations_dir = self.root_dir.parent
        self.stash_root_dir = self.root_dir.parent.parent
        self.stash_output = StashOutput(str(self.stash_root_dir))
        self.name = self.annotation_dir.name


    def required_paths(self) -> dict:
        """
        Returns a dictionary with the required paths for the stash output structure.
        These come on top of self.stash_ouptput.required_paths()
        """
        return {# we don't really need blueprint at this stage anymore
            "annotation": self.annotation_dir,
            "initial_config": self.stash_output.workflow_dir / 'init_nextflow.config'
        }

    def create_structure(self) -> None:
        self.create_directories({'annotation': self.annotation_dir})
    
    def validate_structure(self) -> bool:
        # this is valid if it sits inside annotations of a valid stash output
        valid_structure = self.stash_output.validate_structure()
        required_paths = self.required_paths()
        for pname, path in required_paths.items():
            if not path.exists():
                warnings.warn(f"Missing required path {pname}: {path}")
                valid_structure = False
                break

        try:
            self.validate_label(label=self.name)
        except ValueError as e:
            warnings.warn(f"Invalid annotation name {self.name}: {e}")
            valid_structure = False

        return valid_structure

class AnnotatedUserOutput(BaseOutput):
    """
    Encapsulates the structure for annotation workflows. Example:

      <stash_root_dir>/
      ├── annotations/
      │   └── <any subfolders, e.g. 'testor'>

    """

    def __init__(self, output_dir: str):
        super().__init__(output_dir)
        self.workflow_dir = self.root_dir / "workflow"
        self.workflow_src_dir = self.module_src_dir / "workflow"
        self.name = self.root_dir.name

    def required_paths(self) -> dict:
        """
        Returns a dictionary with the required paths for the stash output structure.
        """
        required_paths = {
            "workflow": self.workflow_dir,
            "workflow_modules": self.workflow_dir / "modules"
        }
        # for path in self.workflow_src_dir.rglob("*"):  # Recursively find all files and dirs
        #     if not path.name.endswith(".config"):  # Exclude .config files
        #         required_paths[f"{path.parent.stem}>{path.name}"] = self.workflow_dir / path.name
        return required_paths

    def create_structure(self) -> None:
        dirs_to_create = {
            "workflow": self.workflow_dir # remaining sub dirs are created by the copytree in VEPDatabase._copy_workflow_srcfiles()
        }
        self.create_directories(dirs_to_create)


    def validate_structure(self) -> bool:
        valid_structure = True

        for pname, path in self.required_paths().items():
            if not path.exists():
                warnings.warn(f"Missing required path {pname}: {path}")
                valid_structure = False
                break

        try:
            self.validate_label(label=self.name)
        except ValueError as e:
            warnings.warn(f"Invalid annotation name {self.name}: {e}")
            valid_structure = False

        return valid_structure
