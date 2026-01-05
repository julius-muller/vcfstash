# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

##########################################################################
#                                                                        #
# This file contains abstract base classes and implementations           #
# for managing output file structures in VCF-based workflows.
# NOT CURRENTLY USED.                                            #
#                                                                        #
##########################################################################
import importlib.resources
import os
import warnings
from abc import ABC, abstractmethod
from pathlib import Path

MAXCHAR_CACHENAME = 50 # maximum amount of characters allowed for a cachename, to keep organization somewhat tidy

class BaseOutput(ABC):
    """Partially abstract base for output structures.
    Subclasses define the exact dirs/files to create or check.
    """

    def __init__(self, root_dir: str):
        self.root_dir = Path(root_dir).expanduser().resolve()
        # define the base directory of the module
        # Use importlib.resources instead of environment variable
        try:
            self.module_src_dir = Path(str(importlib.resources.files("vcfcache")))
        except ModuleNotFoundError:
            self.module_src_dir = Path(os.getenv("VCFCACHE_ROOT", Path(".").resolve()))

    @abstractmethod
    def required_paths(self) -> dict:
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
        """Validates that the label is valid in this context."""
        if len(label) > MAXCHAR_CACHENAME:
            raise ValueError(
                f"Annotation name must be less than {MAXCHAR_CACHENAME} characters, but has {len(label)}: {label}"
            )
        if " " in label:
            raise ValueError(f"Annotation name must not contain white spaces: {label}")
        if not all(c.isalnum() or c in "_-." for c in label):
            raise ValueError(
                f"Annotation name must only contain alphanumeric characters, underscores, dots, or dashes: {label}"
            )

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
                        raise RuntimeError(
                            f"Failed to create {name} directory: {dir_path}"
                        )

                    # Verify it's actually a directory
                    if not dir_path.is_dir():
                        raise RuntimeError(
                            f"Path exists but is not a directory: {dir_path}"
                        )

                    # Verify we have write access by creating and removing a test file
                    test_file = dir_path / ".write_test"
                    try:
                        test_file.touch()
                        test_file.unlink()
                    except (IOError, PermissionError) as e:
                        raise RuntimeError(
                            f"No write permission in {name} directory {dir_path}: {e}"
                        ) from e

                except Exception as e:
                    # Catch any other exceptions that might occur during directory setup
                    raise RuntimeError(
                        f"Error setting up {name} directory {dir_path}: {e}"
                    ) from e


class CacheOutput(BaseOutput):
    """Encapsulates the structure for blueprint-init / blueprint-extend:

      <cache_root_dir>/
      ├── blueprint/
      ├── cache/

    self = CacheOutput(cache_root_dir='.')
    """

    def __init__(self, cache_root_dir: str):
        super().__init__(cache_root_dir)
        self.cache_root_dir = self.root_dir
        self.workflow_dir = self.cache_root_dir / "workflow"
        self.workflow_src_dir = self.module_src_dir / "workflow"

        # We keep the workflow_dir for storing config snapshots

    def required_paths(self) -> dict:
        """Returns a dictionary with the required paths for the cache output structure."""
        return {
            "blueprint": self.cache_root_dir / "blueprint",
            "cache": self.cache_root_dir / "cache",
            "workflow": self.workflow_dir,
            "workflow_src": self.workflow_src_dir
        }

    def create_structure(self) -> None:
        req_dirs = {
            k: v for k, v in self.required_paths().items() if k != "workflow_src"
        }
        self.create_directories(req_dirs)

    def validate_structure(self) -> bool:

        # Minimal existence checks
        required_paths = self.required_paths()

        # Skip it in validation since it's excluded from create_structure
        paths_to_check = {k: v for k, v in required_paths.items() if k != "workflow_src"}

        # for path in self.workflow_src_dir.rglob("*"):  # Recursively find all files and dirs
        #     if not path.name.endswith(".config"):  # Exclude .config files
        #         required_paths[f"{path.parent.stem}>{path.name}"] = self.workflow_dir / path.name

        for pname, path in paths_to_check.items():
            if not path.exists():
                warnings.warn(f"Missing required path {pname}: {path}", stacklevel=2)
                return False
        return True


class AnnotatedCacheOutput(BaseOutput):
    """Encapsulates the structure for annotation cache from cache-build. Example:

    <cache_root_dir>/
    ├── cache/
    │   └── <any subfolders, e.g. 'test'>  <- annotation_dir

    """

    def __init__(self, annotation_dir: str):
        super().__init__(annotation_dir)
        self.annotation_dir = self.root_dir
        self.cache_dir = self.root_dir.parent
        self.cache_root_dir = self.root_dir.parent.parent
        self.cache_output = CacheOutput(str(self.cache_root_dir))
        self.name = self.annotation_dir.name

    def required_paths(self) -> dict:
        """Returns a dictionary with the required paths for the cache output structure.
        These come on top of self.cache_ouptput.required_paths()
        """
        return {  # we don't really need blueprint at this stage anymore
            "annotation": self.annotation_dir,
            "initial_config": self.cache_output.workflow_dir / "init.yaml",
        }

    def create_structure(self) -> None:
        self.create_directories({"annotation": self.annotation_dir})

    def validate_structure(self) -> bool:
        # this is valid if it sits inside cache of a valid cache output
        valid_structure = self.cache_output.validate_structure()
        required_paths = self.required_paths()
        for pname, path in required_paths.items():
            if not path.exists():
                warnings.warn(f"Missing required path {pname}: {path}", stacklevel=2)
                valid_structure = False
                break

        try:
            self.validate_label(label=self.name)
        except ValueError as e:
            warnings.warn(f"Invalid annotation name {self.name}: {e}", stacklevel=2)
            valid_structure = False

        return valid_structure


class AnnotatedUserOutput(BaseOutput):
    """Encapsulates the structure for annotation workflows. Example:

    <cache_root_dir>/
    ├── cache/
    │   └── <any subfolders, e.g. 'testor'>

    """

    def __init__(self, output_dir: str):
        super().__init__(output_dir)
        self.workflow_dir = self.root_dir / "workflow"
        # Dynamically locate the workflow directory in the installed package
        try:
            self.workflow_src_dir = (
                Path(str(importlib.resources.files("vcfcache"))) / "workflow"
            )

        except ModuleNotFoundError:
            raise RuntimeError(
                "Workflow directory not found in the installed package."
            ) from None
        self.name = self.root_dir.name

    def required_paths(self) -> dict:
        """Returns a dictionary with the required paths for the cache output structure."""
        required_paths = {"workflow": self.workflow_dir}
        # for path in self.workflow_src_dir.rglob("*"):  # Recursively find all files and dirs
        #     if not path.name.endswith(".config"):  # Exclude .config files
        #         required_paths[f"{path.parent.stem}>{path.name}"] = self.workflow_dir / path.name
        return required_paths

    def create_structure(self) -> None:
        """Create the required directories for the annotated user output structure."""
        dirs_to_create = {
            "workflow": self.workflow_dir  # remaining sub dirs are created by the copytree in VCFDatabase._copy_workflow_srcfiles()
        }
        self.create_directories(dirs_to_create)

    def validate_structure(self) -> bool:
        """Validate the structure of the annotated user output directory."""
        valid_structure = True

        for pname, path in self.required_paths().items():
            if not path.exists():
                warnings.warn(f"Missing required path {pname}: {path}", stacklevel=2)
                valid_structure = False
                break

        try:
            self.validate_label(label=self.name)
        except ValueError as e:
            warnings.warn(f"Invalid annotation name {self.name}: {e}", stacklevel=2)
            valid_structure = False

        return valid_structure
