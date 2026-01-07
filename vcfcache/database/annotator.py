# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""Database annotation module for the vcfcache package.

This module provides classes for annotating the variant database and for annotating
user VCF files using the annotated database.
"""

import hashlib
import shutil
import tempfile
import subprocess
import sys
import time
from datetime import datetime
from logging import Logger
from multiprocessing import Pool
from pathlib import Path
from typing import Optional, Union

import pysam
import yaml

from vcfcache.database.base import VCFDatabase
from vcfcache.database.outputs import AnnotatedCacheOutput, AnnotatedUserOutput
from vcfcache.utils.logging import setup_logging


class DatabaseAnnotator(VCFDatabase):
    """A class for annotating variant data using a VCF database.

    DatabaseAnnotator is built for managing the annotation workflow of variant data.
    It integrates with a VCF database and leverages configuration files for setting
    up the workflow parameters. The class supports validation of input files, creation
    of necessary directories, and setting up configurations for annotation.

    Attributes:
        annotation_name (str): Name of the annotation operation.
        cached_annotations (AnnotatedCacheOutput): Object managing annotation cached data.
        logger (Logger): Logging utility for the annotator.
        output_dir (Path): Output directory for annotation-related files.
        info_snapshot_file (Path): Path for storing snapshot information about the blueprint in use.
        anno_config_file (Path): Preprocessed annotation configuration file path.
        params_file (Path): Path to the annotation parameters YAML file.
        nx_workflow (WorkflowManager): Object for managing the actual workflow.

    Methods:
        _preprocess_annotation_config(user_config: Path) -> Path
            Preprocesses the given annotation configuration file to address variable substitution issues.

        _validate_inputs() -> None
            Validates the presence and structure of input files, annotation YAML parameters,
            and other critical resources.

        _setup_annotation_cache(force: bool) -> None
            Ensures the cached directory structure is properly set up for the annotation
            process. Deletes existing directories if `force` is set.
    """

    def __init__(
        self,
        annotation_name: str,
        db_path: Path | str,
        anno_config_file: Path | str,
        bcftools_path: Path | str,
        params_file: Optional[Path | str] = None,
        verbosity: int = 0,
        force: bool = False,
        debug: bool = False,
        read_only: bool = False,
    ):
        """Initialize database annotator.

        self = DatabaseAnnotator(annotation_name="testor", anno_config_file=Path('~/projects/vcfcache/tests/config/example_annotation.config'),
         db_path=Path('~/tmp/vcfcache/test_cache'),force=True)

        Args:
            db_path: Path to the database
            params_file: Optional parameters file
            verbosity: Logging verbosity level (0=WARNING, 1=INFO, 2=DEBUG)
            read_only: If True, skip cache setup (for operating on existing cache)
        """
        super().__init__(
            Path(db_path) if isinstance(db_path, str) else db_path, verbosity, debug, bcftools_path
        )

        self.cached_annotations = AnnotatedCacheOutput(
            str(self.cache_dir / annotation_name)
        )
        self.cached_annotations.validate_label(annotation_name)
        self.annotation_name = annotation_name
        self.logger: Logger = self.connect_loggers()

        # Only setup cache if not read-only mode
        if not read_only:
            self._setup_annotation_cache(force)
        self.output_dir = self.cached_annotations.annotation_dir

        self.info_snapshot_file = self.output_dir / "blueprint_snapshot.info"

        self.anno_config_file = self._preprocess_annotation_config(
            user_config=Path(anno_config_file).expanduser().resolve()
        )

        self.params_file = self.output_dir / "params.snapshot.yaml"
        if params_file:
            params_path = (
                Path(params_file) if isinstance(params_file, str) else params_file
            )
            shutil.copyfile(params_path.expanduser().resolve(), self.params_file)
        else:
            wfi = self.workflow_dir / "init.yaml"
            assert wfi.exists(), f"Workflow init params file not found: {wfi}"
            shutil.copyfile(wfi, self.params_file)
            assert (


                self.params_file.exists()
            ), f"Workflow params file not found: {self.params_file}"

        # Initialize workflow backend (pure Python)
        from vcfcache.database.base import create_workflow
        self.nx_workflow = create_workflow(
            input_file=self.blueprint_bcf,
            output_dir=self.output_dir,
            name=self.annotation_name,
            anno_config_file=self.anno_config_file,
            params_file=self.params_file,
            verbosity=self.verbosity,
        )

        self._validate_inputs()

        # Log initialization parameters
        if self.logger:
            self.logger.info("Initializing database annotation")
            self.logger.debug(f"Annotation directory: {self.output_dir}")

    def _preprocess_annotation_config(self, user_config: Path) -> Path:
        """Preprocess annotation.yaml to fix variable substitution issues.
        Replaces problematic variable references with their escaped versions.
        This only has to be done once, as the config is copied to the output directory and used for all subsequent vcfcache annotate runs

        Args:
            anno_config_file (Path): Original annotation config file
            output_config_file (Path): Where to save the processed config
        """
        assert user_config.exists(), f"Annotation config file not found: {user_config}"

        with open(user_config, "r") as f:
            content = f.read()

        # Use regex to ensure exactly one backslash before variables
        # This fixes files from old buggy code that may have double backslashes
        import re

        # First, fix any double (or more) backslashes before variables (from old buggy code)
        cleanup_patterns = [
            (r'\\+\$\{INPUT_BCF\}', '\\${INPUT_BCF}'),  # Replace \\+${...} with \${...}
            (r'\\+\$INPUT_BCF(?![_{])', '\\$INPUT_BCF'),
            (r'\\+\$\{OUTPUT_BCF\}', '\\${OUTPUT_BCF}'),
            (r'\\+\$OUTPUT_BCF(?![_{])', '\\$OUTPUT_BCF'),
            (r'\\+\$\{AUXILIARY_DIR\}', '\\${AUXILIARY_DIR}'),
            (r'\\+\$AUXILIARY_DIR(?![_{])', '\\$AUXILIARY_DIR'),
        ]

        modified_content = content
        for pattern, replacement in cleanup_patterns:
            before = modified_content
            modified_content = re.sub(pattern, replacement, modified_content)
            if before != modified_content and self.logger:
                self.logger.debug(f"Cleanup: Fixed double backslashes for pattern: {pattern}")

        # Then, add backslashes where missing (for unescaped variables)
        # Pattern explanation: negative lookbehind (?<!\\) ensures we don't match if backslash already present
        # Note: Use regular strings for replacements, not raw strings, so \\$ = one backslash
        add_escape_patterns = [
            (r'(?<!\\)\$\{INPUT_BCF\}', '\\${INPUT_BCF}'),
            (r'(?<!\\)\$INPUT_BCF(?![_{])', '\\$INPUT_BCF'),  # Don't match if followed by _ or {
            (r'(?<!\\)\$\{OUTPUT_BCF\}', '\\${OUTPUT_BCF}'),
            (r'(?<!\\)\$OUTPUT_BCF(?![_{])', '\\$OUTPUT_BCF'),
            (r'(?<!\\)\$\{AUXILIARY_DIR\}', '\\${AUXILIARY_DIR}'),
            (r'(?<!\\)\$AUXILIARY_DIR(?![_{])', '\\$AUXILIARY_DIR'),
        ]

        # Apply each regex replacement
        for pattern, replacement in add_escape_patterns:
            before = modified_content
            modified_content = re.sub(pattern, replacement, modified_content)
            if before != modified_content and self.logger:
                self.logger.debug(f"Add escape: Added backslash for pattern: {pattern}")

        if self.logger:
            # Debug: Check final state of variables in the preprocessed content
            import_bcf_count = modified_content.count('\\${INPUT_BCF}')
            double_backslash_count = modified_content.count('\\\\${INPUT_BCF}')
            self.logger.debug(f"Preprocessing result: \\${{INPUT_BCF}} count={import_bcf_count}, \\\\${{INPUT_BCF}} count={double_backslash_count}")

        output_cfg = self.output_dir / "annotation.yaml"
        with open(output_cfg, "w") as f:
            f.write(modified_content)

        return output_cfg

    def _validate_inputs(self) -> None:
        """Validate input files, directories, and YAML parameters."""
        if self.logger:
            self.logger.debug("Validating inputs")

        if self.logger:
            self.logger.debug("Input validation successful")

    def _setup_annotation_cache(self, force: bool) -> None:
        # Remove destination directory if it exists to ensure clean copy
        if self.cached_annotations.annotation_dir.exists():
            if (
                self.cached_output.validate_structure()
            ):  # we dont want to remove a random dir....
                if force:
                    print(
                        f"Cache directory already exists, removing: {self.cached_annotations.annotation_dir}"
                    )
                    shutil.rmtree(self.cached_annotations.annotation_dir)
                else:
                    raise FileExistsError(
                        f"Output directory already exists: {self.cached_annotations.annotation_dir}\nIf intended, use --force to overwrite."
                    )
            else:
                if not force:
                    raise FileNotFoundError(
                        f"Output directory must not exist if --force is not set and a valid cache directory: {self.cached_annotations.annotation_dir}"
                    )

        if self.logger:
            self.logger.debug(f"Creating cache structure: {self.cached_annotations.annotation_dir}")
        self.cached_annotations.create_structure()

    def annotate(self, extra_files: bool = True) -> None:
        """Run annotation workflow on database"""
        # Store blueprint snapshot and workflow files
        shutil.copy2(self.info_file, self.info_snapshot_file)

        try:
            if self.logger:
                self.logger.info("Starting annotation workflow")

            start_time = datetime.now()
            self.nx_workflow.run(
                db_mode="cache-build",
                db_bcf=self.blueprint_bcf,
                trace=extra_files,
                dag=extra_files,
                report=extra_files,
            )
            if not self.debug:
                self.nx_workflow.cleanup_work_dir()

            duration = datetime.now() - start_time
            if self.logger:
                self.logger.info(
                    f"Annotation to {self.output_dir} completed in {duration.total_seconds():.2f} seconds"
                )

            # Write completion flag
            from vcfcache.utils.completion import write_completion_flag
            write_completion_flag(
                output_dir=self.cached_annotations.root_dir,
                command="cache-build",
                mode="cache-build"
            )

        except subprocess.CalledProcessError as e:
            if self.logger:
                self.logger.error(f"Workflow execution failed: {e.stderr}")
                self.logger.warning(f"Removing output directory: {self.output_dir}")
            shutil.rmtree(self.output_dir)
            sys.exit(1)


class VCFAnnotator(VCFDatabase):
    """Provides functionality to annotate VCF/BCF files using a predefined annotation
    database while managing output directories and workflows.

    This class is designed for handling large-scale genomic data annotations.
    It requires an annotation database and an input VCF/BCF
    file, both of which should be prepared and structured correctly beforehand.

    Attributes:
        input_vcf (Path): Path to the input VCF/BCF file.
        annotation_db_path (Path): Path to the structured annotation database directory.
        annotation_name (str): Name derived from the annotation database.
        output_dir (Path): Directory where annotation-related files are stored.
        params_file (Path): Path to the parameters YAML file.
        logger (Logger): Logging instance for the class.
        nx_workflow (WorkflowManager): Instance of the workflow manager.

    Args:
        input_vcf (Path | str): Path to the input VCF/BCF file, which must be indexed.
        annotation_db (Path | str): Path to the annotation database.
        output_file (Path | str): Path where annotated results will be written (or '-'/'stdout').
        params_file (Optional[Path | str]): Path to a custom parameters file; if
            not provided, defaults to "annotation.yaml" from the annotation database.
        verbosity (int): Logging verbosity level; 0 (WARNING), 1 (INFO), or 2 (DEBUG).
        force (bool): Indicates whether to overwrite the existing output directory.
        debug (bool): Enables extensive logging and debug mode for troubleshooting.

    Raises:
        FileNotFoundError: Raised if input VCF/BCF, annotation database, or
            required reference/configuration files are not found.
        FileExistsError: Raised if output directory exists but force flag is not set.
        ValueError: Raised if invalid parameters or VCF reference information
            is encountered during validation.
        RuntimeError: Raised for unexpected errors during validation or setup.
    """

    VALID_VCF_EXTENSIONS = (".bcf", ".vcf.gz", ".vcf")
    INFO_FIELDS = [
        "GT",
        "DP",
        "AF",
        "gnomadg_af",
        "gnomade_af",
        "gnomadg_ac",
        "gnomade_ac",
        "clinvar_clnsig",
        "deeprvat_score",
    ]
    BASES = {"A", "C", "G", "T"}

    def __init__(
        self,
        input_vcf: Path | str,
        annotation_db: Path | str,
        output_file: Path | str,
        bcftools_path: Path,
        params_file: Optional[Path | str] = None,
        stats_dir: Optional[Path | str] = None,
        no_stats: bool = False,
        verbosity: int = 0,
        force: bool = False,
        debug: bool = False,
    ):
        """Initialize database annotator.

        Args:
            input_vcf: Path to the input BCF/VCF file, needs to be indexed!
            annotation_db: Path to the annotation database
            output_file: Path to the output BCF file, or '-'/'stdout' to stream
            force: Whether to overwrite existing output directory
            debug: Whether to enable debug mode
            verbosity: Logging verbosity level (0=WARNING, 1=INFO, 2=DEBUG)
        """
        self.input_vcf = Path(input_vcf).expanduser().resolve()
        self.no_stats = no_stats
        self.vcf_name, fext = self._validate_and_extract_sample_name()

        if not self.input_vcf.exists():
            raise FileNotFoundError(f"Input VCF file not found: {self.input_vcf}")

        self.cached_annotations = AnnotatedCacheOutput(str(annotation_db))
        if not self.cached_annotations.validate_structure():
            raise FileNotFoundError(
                f"Annotation database annotation_db not valid: {self.cached_annotations.annotation_dir}"
            )
        self.annotation_db_path = self.cached_annotations.annotation_dir
        self.annotation_name = self.cached_annotations.name
        super().__init__(
            self.cached_annotations.cache_output.root_dir, verbosity, debug, bcftools_path
        )

        self.output_to_stdout = str(output_file) in {"-", "stdout"}
        self.output_vcf: Optional[Path]
        if self.output_to_stdout:
            self.output_vcf = None
            output_name = "stdout"
        else:
            output_path = Path(output_file).expanduser()
            if output_path.exists() and output_path.is_dir():
                raise ValueError(
                    f"Output path is a directory: {output_path}\n"
                    "Provide a file path (e.g., sample_vc.bcf) and use --stats-dir for logs."
                )
            if output_path.suffix == "":
                output_path = output_path.with_suffix(".bcf")
            if output_path.suffix != ".bcf":
                raise ValueError(
                    f"Output file must end with .bcf or be '-'/'stdout': {output_file}"
                )
            output_path = output_path.resolve()
            if output_path.exists():
                if force:
                    output_path.unlink()
                    csi_path = output_path.with_suffix(output_path.suffix + ".csi")
                    if csi_path.exists():
                        csi_path.unlink()
                else:
                    raise FileExistsError(
                        f"Output file already exists: {output_path}\nIf intended, use --force to overwrite."
                    )
            output_path.parent.mkdir(parents=True, exist_ok=True)
            self.output_vcf = output_path
            output_name = output_path.name

        if self.no_stats and stats_dir:
            raise ValueError("--no-stats cannot be used together with --stats-dir.")

        input_stats_name = f"{self._input_basename()}_vcstats"
        if self.no_stats:
            stats_output_dir = Path(tempfile.mkdtemp(prefix="vcfcache_stats_")).resolve()
        elif stats_dir:
            stats_output_dir = Path(stats_dir).expanduser().resolve() / input_stats_name
        else:
            stats_output_dir = Path.cwd() / input_stats_name

        self.stats_dir = stats_output_dir
        self.output_annotations = AnnotatedUserOutput(str(stats_output_dir))
        self.output_annotations.validate_label(self.output_annotations.name)
        self._setup_output(force=force)
        self.output_dir = self.output_annotations.root_dir

        self.logger: Logger = setup_logging(
            verbosity=self.verbosity,
            log_file=self.output_dir / "annotation.log",
        )

        self.annotation_wfl_path = self.output_annotations.workflow_dir
        self.annotation_wfl_path.mkdir(parents=True, exist_ok=True)

        # now also import the mandatory annotation file, that cannot be provided by the user at this stage
        self.anno_config_file = self.annotation_db_path / "annotation.yaml"
        if not self.anno_config_file.exists():
            raise FileNotFoundError(
                f"Annotation config file not found: {self.anno_config_file}"
            )
        self.anno_snapshot_file = self.annotation_wfl_path / "annotation.snapshot.yaml"
        shutil.copyfile(self.anno_config_file, self.anno_snapshot_file)

        self.params_file = self.annotation_wfl_path / "params.snapshot.yaml"
        if params_file:
            params_path = (
                Path(params_file) if isinstance(params_file, str) else params_file
            )
            shutil.copyfile(params_path.expanduser().resolve(), self.params_file)
        else:
            wfi = self.annotation_db_path / "params.snapshot.yaml"
            assert (
                wfi.exists()
            ), f"Workflow annotation params file not found: {wfi}, required if no yaml provided!"
            shutil.copyfile(wfi, self.params_file)
            print(f"Using default params from cache: {wfi.relative_to(wfi.parent.parent.parent)}")
        assert (
            self.params_file.exists()
        ), f"Workflow params file not found: {self.params_file}"

        self.cache_file = self.annotation_db_path / "vcfcache_annotated.bcf"
        if not self.cache_file.exists():
            raise FileNotFoundError(f"Cache file not found: {self.cache_file}")

        # Initialize workflow backend (pure Python)
        from vcfcache.database.base import create_workflow
        self.nx_workflow = create_workflow(
            input_file=self.input_vcf,
            output_dir=self.output_dir,
            name=self.annotation_name,
            anno_config_file=self.anno_config_file,
            params_file=self.params_file,
            verbosity=self.verbosity,
            output_file=self.output_vcf if not self.output_to_stdout else "-",
        )

        self._validate_inputs()

        # Log initialization parameters
        if self.logger:
            self.logger.info(f"Initializing annotation of {self.input_vcf.name}")
            self.logger.debug(f"Cache file: {self.cache_file}")

    def _log_contig_overlap(self) -> None:
        """Log contig overlap between cache and input, and error on no overlap."""
        cache_contigs = self._list_contigs(self.cache_file)
        input_contigs = self._list_contigs(self.input_vcf)
        cache_set = set(cache_contigs)
        input_set = set(input_contigs)
        overlap = sorted(cache_set & input_set)

        msg = (
            f"Contig overlap: cache={len(cache_set)} input={len(input_set)} "
            f"overlap={len(overlap)}"
        )
        if self.logger:
            self.logger.info(msg)
        else:
            print(msg)

        if overlap:
            preview = ", ".join(overlap[:10])
            preview_msg = f"Overlapping contigs (first 10): {preview}"
            if self.logger:
                self.logger.info(preview_msg)
            else:
                print(preview_msg)
            return

        raise RuntimeError(
            "No contigs in common between cache and input. "
            "Check that both use the same reference genome and contig naming."
        )

    def _validate_inputs(self) -> None:
        """Validate input files, directories, and YAML parameters."""
        if self.logger:
            self.logger.debug("Validating inputs")

        # Check input VCF/BCF if provided
        if self.input_vcf:
            if not self.input_vcf.exists():
                msg = f"Input VCF/BCF file not found: {self.input_vcf}"
                if self.logger:
                    self.logger.error(msg)
                raise FileNotFoundError(msg)
            self.ensure_indexed(self.input_vcf)

        if self.logger:
            self.logger.debug("Input validation successful")

    def _setup_output(self, force: bool) -> None:
        # Remove destination directory if it exists to ensure clean copy
        if self.output_annotations.root_dir.exists():
            if (
                self.output_annotations.validate_structure()
            ):  # we dont want to remove a random dir....
                if force:
                    print(
                        f"Output directory already exists, removing: {self.output_annotations.root_dir}"
                    )
                    shutil.rmtree(self.output_annotations.root_dir)
                else:
                    raise FileExistsError(
                        f"Output directory already exists: {self.output_annotations.root_dir}\nIf intended, use --force to overwrite."
                    )
            else:
                raise FileNotFoundError(
                    f"Output directory must not exist if --force is not set and a valid output directory: {self.output_annotations.root_dir}"
                )

        if self.logger:
            self.logger.debug(f"Creating output structure: {self.output_annotations.root_dir}")
        self.output_annotations.create_structure()

    def _validate_and_extract_sample_name(self) -> tuple[str, str]:
        """Validates the input VCF file has an acceptable extension
        ('.bcf', '.vcf.gz', '.vcf') and extracts the sample name
        (filename without directory path and extension). Also checks that the file is indexed.

        Returns:
            str: the extracted sample name

        Raises:
            ValueError: if the input file has an invalid extension
            FileNotFoundError: if the index file is missing
        """
        input_vcf_path = self.input_vcf

        # Validate file extension
        if not input_vcf_path.suffixes:
            raise ValueError(
                f"Input VCF file '{input_vcf_path}' lacks a file extension."
            )

        # Check for valid extensions, considering multi-part extensions
        # vcfcache internally uses BCF format, so we prefer .csi indices
        if input_vcf_path.name.endswith(".vcf.gz"):
            extension = ".vcf.gz"
            sample_name = input_vcf_path.name[:-7]  # Removes '.vcf.gz'
            # Check for .csi first (BCF-compatible), then .tbi (legacy VCF)
            index_file_csi = Path(str(input_vcf_path) + ".csi")
            index_file_tbi = input_vcf_path.with_suffix(input_vcf_path.suffix + ".tbi")
            if index_file_csi.exists():
                index_file = index_file_csi
            elif index_file_tbi.exists():
                index_file = index_file_tbi
            else:
                index_file = index_file_csi  # For error message
        elif input_vcf_path.name.endswith(".bcf"):
            extension = ".bcf"
            sample_name = input_vcf_path.name[:-4]  # Removes '.bcf'
            index_file = input_vcf_path.with_suffix(".bcf.csi")
        elif input_vcf_path.name.endswith(".vcf"):
            extension = ".vcf"
            sample_name = input_vcf_path.name[:-4]  # Removes '.vcf'
            # Check for .csi first (BCF-compatible), then .tbi (legacy VCF)
            index_file_csi = input_vcf_path.with_suffix(".vcf.csi")
            index_file_tbi = input_vcf_path.with_suffix(".vcf.tbi")
            if index_file_csi.exists():
                index_file = index_file_csi
            elif index_file_tbi.exists():
                index_file = index_file_tbi
            else:
                index_file = index_file_csi  # For error message
        else:
            raise ValueError(
                f"Input VCF file '{input_vcf_path}' must end with one of {self.VALID_VCF_EXTENSIONS}"
            )

        if not index_file.exists():
            raise FileNotFoundError(
                f"Index file for '{input_vcf_path}' not found: '{index_file}' (vcfcache uses .csi indices for BCF-compatible operations)"
            )

        return sample_name, extension

    def _process_region(self, args: tuple) -> "pd.DataFrame":
        """Process a single genomic region from BCF file."""
        try:
            import pandas as pd

            bcf_path, region = args
            if self.logger:
                self.logger.debug(f"Processing region: {region}")

            vcf = pysam.VariantFile(str(bcf_path))
            records = []
            variant_count = 0
            excluded_count = 0

            try:
                for record in vcf.fetch(region=region):
                    variant_count += 1
                    try:
                        # Extract basic variant fields
                        chrom = record.chrom
                        if chrom[:3] != "chr":
                            continue

                        pos = record.pos
                        ref = record.ref
                        alt = (
                            record.alts[0] if record.alts else ""
                        )  # Assuming single ALT
                        if not all([x in self.BASES for x in alt]):
                            continue

                        # Process INFO fields
                        info = {key: None for key in self.INFO_FIELDS}
                        info |= {
                            key: record.info.get(key, None)
                            for key in self.INFO_FIELDS
                            if key in record.info
                        }

                        # Extract FORMAT fields
                        if len(record.samples):

                            sample = record.samples[0]
                            ad = sample.get("AD", None)
                            dp = sample.get("DP", None)

                            # Calculate AF
                            af = None
                            if ad and len(ad) >= 2:
                                ref_depth = ad[0]
                                alt_depth = ad[1]
                                af = (
                                    alt_depth / (ref_depth + alt_depth)
                                    if (ref_depth + alt_depth) > 0
                                    else None
                                )

                            info |= {
                                "GT": sample.get("GT", None),
                                "AD": ad[1] if ad and len(ad) > 1 else None,
                                "DP": dp,
                                "AF": af,
                            }

                        # Process clinvar and gnomad fields
                        self._process_variant_annotations(info)

                        # Process VCF annotations
                        # TODO: currently does not work for any Tag apart from CSQ, need to pull from .config via ${params.must_contain_info_tag}
                        if "CSQ" in record.info:
                            csq_values = record.info["CSQ"]
                            if csq_values:
                                description = vcf.header.info["CSQ"].description
                                if description:
                                    fields = (
                                        description.split(" ")[-1].split("|")
                                        if " " in description
                                        else []
                                    )
                                    vcf_csqs = [
                                        dict(zip(fields, x.split("|"), strict=False))
                                        for x in csq_values
                                    ]
                                    # TODO: This currently only works with VEP annotations, need to be more flexible
                                    expanded_transcripts = self.parse_vcf_info(vcf_csqs)

                                    for transcript in expanded_transcripts:
                                        row = {
                                            "CHROM": chrom,
                                            "POS": pos,
                                            "REF": ref,
                                            "ALT": alt,
                                            **info,
                                            **transcript,
                                        }
                                        records.append(row)

                    except Exception as e:
                        excluded_count += 1
                        if self.logger:
                            self.logger.error(
                                f"Error processing variant at {record.chrom}:{record.pos}: {e}"
                            )
                        continue

                if self.logger:
                    self.logger.debug(
                        f"Processed {variant_count} variants, excluded {excluded_count} in {region}"
                    )
                return pd.DataFrame(records)

            except Exception as e:
                if self.logger:
                    self.logger.error(f"Error processing region {region}: {e}")
                raise
        except ImportError:
            raise ImportError(
                "This feature requires pandas and pyarrow. "
                "Install with 'pip install pandas pyarrow'"
            ) from None

    def _process_variant_annotations(self, info: dict) -> None:
        """Process clinvar and gnomad annotations."""
        # Process clinvar
        clinvar_clnsig = info.get("clinvar_clnsig", None)
        if clinvar_clnsig and clinvar_clnsig[0] != "null":
            info["clinvar_clnsig"] = ", ".join(clinvar_clnsig)
        else:
            info["clinvar_clnsig"] = None

        # Process gnomad fields
        gnomad_fields = ["gnomadg_ac", "gnomade_ac", "gnomadg_af", "gnomade_af"]
        for field in gnomad_fields:
            value = info.get(field, None)
            if isinstance(value, (int, float)) and value < 0:
                value = None
            info[field] = value

        # Calculate weighted average
        gnomadg_ac = float(info["gnomadg_af"]) if info.get("gnomadg_ac", None) else 0
        gnomade_ac = (
            int(info.get("gnomade_ac", 0)) if info.get("gnomade_ac", None) else 0
        )
        gnomadg_af = float(info["gnomadg_af"]) if info.get("gnomadg_af", None) else None
        gnomade_af = float(info["gnomade_af"]) if info.get("gnomade_af", None) else None

        info["gnomad_af"] = self.wavg(gnomadg_af, gnomade_af, gnomadg_ac, gnomade_ac)

        # Remove individual gnomad fields
        for field in gnomad_fields:
            info.pop(field, None)

    def annotate(self, uncached: bool = False, convert_parquet: bool = False, preserve_unannotated: bool = False, skip_split_multiallelic: bool = False, md5_all: bool = False) -> None:
        """Run annotation workflow on input VCF file.

        Args:
            uncached: Whether to run the workflow in uncached mode
            convert_parquet: Whether to convert output to Parquet format
            preserve_unannotated: Whether to preserve variants without annotation in output
            skip_split_multiallelic: Skip splitting multiallelic variants (use only if certain input has none)
            md5_all: Compute full MD5 of all variants (no header) and store in stats

        Returns:
            Path to output file (BCF or Parquet)
            self = VCFAnnotator(input_vcf="~/projects/vcfcache/tests/data/nodata/sample4.bcf",
             annotation_db="~/tmp/test/test_out/cache/testor", output_file="~/tmp/test/aout.bcf", stats_dir="~/tmp/test/aout_stats", force=True)

        """
        start_time = time.time()
        if self.logger:
            self.logger.debug("Starting VCF annotation")

        try:
            if not uncached:
                self._log_contig_overlap()

            # Run the workflow in database mode
            self.nx_workflow.run(
                db_mode="annotate" if not uncached else "annotate-nocache",
                db_bcf=self.cache_file,
                trace=True,
                dag=True,
                report=True,
                preserve_unannotated=preserve_unannotated,
                skip_split_multiallelic=skip_split_multiallelic,
            )
            duration = time.time() - start_time
            # Always show completion (even in default mode)
            print(f"Annotation completed in {duration:.1f}s")

            if not self.no_stats:
                from vcfcache.utils.completion import write_completion_flag
                mode = "uncached" if uncached else "cached"
                write_completion_flag(
                    output_dir=self.output_annotations.root_dir,
                    command="annotate",
                    mode=mode,
                    output_file=str(self.output_vcf) if self.output_vcf else "stdout",
                )
                self._write_compare_stats(mode=mode, md5_all=md5_all)

            if convert_parquet:
                if self.output_vcf is None:
                    raise ValueError("Parquet conversion is not supported when output is stdout.")
                # threads = self.nx_workflow.nf_config_content['params'].get('vep_max_forks',1) * self.nx_workflow.nf_config_content['params'].get('vep_max_chr_parallel', 1)
                self._convert_to_parquet(self.output_vcf)  # , threads=threads)

            if not self.debug:
                self.nx_workflow.cleanup_work_dir()

        except Exception:
            if self.logger:
                self.logger.error("Annotation failed", exc_info=True)
            raise
        finally:
            if self.no_stats:
                shutil.rmtree(self.output_dir, ignore_errors=True)

    def _input_basename(self) -> str:
        name = self.input_vcf.name
        suffixes = "".join(self.input_vcf.suffixes)
        if suffixes:
            return name[: -len(suffixes)]
        return self.input_vcf.stem

    def _write_compare_stats(self, mode: str, md5_all: bool = False) -> None:
        stats_file = self.output_dir / "compare_stats.yaml"
        output_file = str(self.output_vcf) if self.output_vcf else "stdout"
        anno_md5 = None
        if self.anno_snapshot_file.exists():
            from vcfcache.utils.validation import compute_md5
            anno_md5 = compute_md5(self.anno_snapshot_file)

        threads = None
        if hasattr(self.nx_workflow, "params_file_content"):
            threads = self.nx_workflow.params_file_content.get("threads")

        stats = {
            "command": "annotate",
            "mode": mode,
            "run_timestamp": datetime.now().isoformat(timespec="seconds"),
            "output_file": output_file,
            "input_file": str(self.input_vcf),
            "input_name": self.input_vcf.name,
            "cache_name": self.annotation_name,
            "cache_path": str(self.annotation_db_path),
            "annotation_yaml_md5": anno_md5,
            "genome_build_params": self.nx_workflow.params_file_content.get("genome_build"),
            "genome_build_annotation": self.nx_workflow.nfa_config_content.get("genome_build"),
            "vcfcache_version": getattr(sys.modules.get("vcfcache"), "__version__", "unknown"),
            "threads": threads,
            "variant_counts": {
                "total_output": None,
                "annotated_output": None,
                "tool_annotated": None,
                "input_variants": None,
                "dropped_variants": None,
                "missing_variants": None,
                "missing_annotated": None,
            },
            "variant_md5": {
                "top10": None,
                "bottom10": None,
                "all": None,
            },
        }

        last_stats = getattr(self.nx_workflow, "last_run_stats", {}) or {}
        stats["variant_counts"]["dropped_variants"] = last_stats.get("dropped_variants")
        stats["variant_counts"]["missing_variants"] = last_stats.get("missing_variants")
        stats["variant_counts"]["missing_annotated"] = last_stats.get("missing_annotated")
        stats["variant_counts"]["input_variants"] = last_stats.get("input_variants")

        if self.output_vcf and self.output_vcf.exists():
            stats["variant_counts"]["total_output"] = self._count_variants(self.output_vcf)
            tag = self.nx_workflow.nfa_config_content.get("must_contain_info_tag")
            stats["variant_counts"]["annotated_output"] = self._count_annotated_variants(self.output_vcf, tag)
            top_md5, bottom_md5 = self._compute_top_bottom_md5(self.output_vcf)
            stats["variant_md5"]["top10"] = top_md5
            stats["variant_md5"]["bottom10"] = bottom_md5
            if md5_all:
                stats["variant_md5"]["all"] = self._compute_all_md5(self.output_vcf)

        tool_annotated = None
        if mode == "cached":
            tool_annotated = last_stats.get("missing_annotated")
        elif mode == "uncached":
            tool_annotated = last_stats.get("output_variants")
        if tool_annotated is None:
            tool_annotated = stats["variant_counts"]["annotated_output"]
        stats["variant_counts"]["tool_annotated"] = tool_annotated

        stats_file.write_text(yaml.safe_dump(stats, sort_keys=False))

    def _count_variants(self, bcf_path: Path) -> Optional[int]:
        result = subprocess.run(
            [str(self.bcftools_path), "index", "-n", str(bcf_path)],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            return None
        try:
            return int(result.stdout.strip())
        except ValueError:
            return None

    def _count_annotated_variants(self, bcf_path: Path, tag: Optional[str]) -> Optional[int]:
        if not tag:
            return None
        result = subprocess.run(
            [str(self.bcftools_path), "view", "-H", "-i", f"INFO/{tag}!=\"\"", str(bcf_path)],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            return None
        return len([line for line in result.stdout.splitlines() if line.strip()])

    def _compute_top_bottom_md5(self, bcf_path: Path) -> tuple[Optional[str], Optional[str]]:
        try:
            proc = subprocess.Popen(
                [str(self.bcftools_path), "view", "-H", str(bcf_path)],
                stdout=subprocess.PIPE,
                text=True,
            )
        except Exception:
            return None, None

        top_lines = []
        bottom_lines = []
        for line in proc.stdout:  # type: ignore[union-attr]
            line = line.rstrip("\n")
            if len(top_lines) < 10:
                top_lines.append(line)
            bottom_lines.append(line)
            if len(bottom_lines) > 10:
                bottom_lines.pop(0)
        proc.wait()

        if not top_lines:
            return None, None

        def _md5(lines):
            h = hashlib.md5()
            for l in lines:
                h.update((l + "\n").encode())
            return h.hexdigest()

        top_md5 = _md5(top_lines)
        bottom_md5 = _md5(bottom_lines)
        if len(top_lines) <= 10 and len(bottom_lines) <= 10 and len(top_lines) == len(bottom_lines):
            bottom_md5 = top_md5
        return top_md5, bottom_md5

    def _compute_all_md5(self, bcf_path: Path) -> Optional[str]:
        try:
            proc = subprocess.Popen(
                [str(self.bcftools_path), "view", "-H", str(bcf_path)],
                stdout=subprocess.PIPE,
                text=True,
            )
        except Exception:
            return None

        md5_hash = hashlib.md5()
        for line in proc.stdout:  # type: ignore[union-attr]
            md5_hash.update(line.encode())
        proc.wait()
        if proc.returncode != 0:
            return None
        return md5_hash.hexdigest()

    def _convert_to_parquet(self, bcf_path: Path, threads: int = 2) -> Path:
        """Convert annotated BCF to optimized Parquet format"""
        try:
            import pandas as pd
            import pyarrow as pa  # type: ignore
            import pyarrow.parquet as pq  # type: ignore

            vcf = pysam.VariantFile(str(bcf_path))
            regions = list(vcf.header.contigs.keys())
            args_list = [(str(bcf_path), region) for region in regions]

            if self.logger:
                self.logger.info(f"Converting BCF to Parquet: {bcf_path}")
                self.logger.debug(
                    f"Processing {len(regions)} regions using {threads} threads"
                )

            with Pool(threads) as pool:
                dataframes = pool.map(self._process_region, args_list)

            # Filter and combine dataframes
            dataframes = [df for df in dataframes if not df.empty]
            if not dataframes:
                if self.logger:
                    self.logger.error("No valid variants found in annotated file")
                raise ValueError("No valid variants found in annotated file")

            combined_df = pd.concat(dataframes, ignore_index=True)
            output_vcf = getattr(self, "output_vcf", None)
            output_base = output_vcf.parent if output_vcf else self.output_dir
            output_file = output_base / f"{self.input_vcf.stem}.parquet"
            if self.logger:
                self.logger.info(f"Writing Parquet file: {output_file}")
                self.logger.debug(f"Total variants: {len(combined_df)}")

            # Write optimized parquet
            table = pa.Table.from_pandas(combined_df)
            pq.write_table(
                table,
                output_file,
                compression="snappy",
                use_dictionary=True,
                row_group_size=100000,
                data_page_size=65536,
                write_statistics=True,
            )
            if self.logger:
                self.logger.info("Parquet conversion completed")
            return output_file

        except ImportError:
            raise ImportError(
                "Converting to Parquet requires additional dependencies. "
                "Please install them with: pip install pandas pyarrow"
            ) from None

    @staticmethod
    def wavg(
        f1: float | None, f2: float | None, n1: Union[int, float], n2: Union[int, float]
    ) -> float | None:
        """Weighted average for Allele Frequencies."""
        total_weight = n1 + n2
        if total_weight == 0:
            return None
        if f1 is not None and f2 is not None:
            return (f1 * n1 + f2 * n2) / total_weight
        elif f1 is None and f2 is None:
            return None
        elif f1 is None:
            return f2
        else:  # f2 is None
            return f1
