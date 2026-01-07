# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

"""Pure Python workflow manager for VCFcache operations.

This module implements a pure Python alternative to the Nextflow-based workflow
system. It executes the same bcftools commands directly via subprocess, eliminating
the need for Java and Nextflow dependencies.
"""

import datetime
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Union

import yaml

from vcfcache.database.workflow_base import WorkflowBase
from vcfcache.utils.logging import setup_logging


class BcftoolsCommand:
    """Helper class to execute and log bcftools commands."""

    def __init__(
        self, cmd: str, logger, work_dir: Path, log_file: Optional[Path] = None
    ):
        """Initialize a bcftools command.

        Args:
            cmd: The command to execute (bash script)
            logger: Logger instance for output
            work_dir: Working directory for execution
            log_file: Optional path to log file
        """
        self.cmd = cmd
        self.logger = logger
        self.work_dir = work_dir
        self.log_file = log_file or (work_dir / "command.log")

    def run(self, check: bool = True) -> subprocess.CompletedProcess:
        """Execute the command with logging.

        Args:
            check: If True, raise exception on non-zero exit code

        Returns:
            subprocess.CompletedProcess with results

        Raises:
            subprocess.CalledProcessError: If command fails and check=True
        """
        import time

        self.logger.debug(f"Running command in {self.work_dir}")
        self.logger.debug(f"Command: {self.cmd}")

        # Write command to file for debugging
        script_file = self.work_dir / "command.sh"
        script_file.write_text(f"#!/bin/bash\nset -euo pipefail\n\n{self.cmd}\n")
        script_file.chmod(0o755)

        # Start timing
        start_time = time.time()

        # Execute command
        result = subprocess.run(
            self.cmd,
            shell=True,
            capture_output=True,
            text=True,
            cwd=self.work_dir,
            executable="/bin/bash",
        )

        # Calculate duration
        duration = time.time() - start_time

        # Extract command name (first word after bcftools)
        cmd_name = "unknown"
        if "bcftools" in self.cmd:
            parts = self.cmd.split()
            for i, part in enumerate(parts):
                if "bcftools" in part and i + 1 < len(parts):
                    cmd_name = f"bcftools {parts[i+1].split()[0]}"
                    break

        # Log timing
        self.logger.info(f"Command completed in {duration:.3f}s: {cmd_name}")

        # Save timing info
        timing_file = self.work_dir / "timing.txt"
        with timing_file.open("a") as f:
            f.write(f"{cmd_name}\t{duration:.3f}\n")

        # Save stdout/stderr
        (self.work_dir / "stdout.txt").write_text(result.stdout)
        (self.work_dir / "stderr.txt").write_text(result.stderr)

        # Log output
        if result.stdout:
            self.logger.debug(f"STDOUT: {result.stdout}")
        if result.stderr:
            self.logger.debug(f"STDERR: {result.stderr}")

        if result.returncode != 0:
            self.logger.error(f"Command failed with exit code {result.returncode}")
            self.logger.error(f"STDERR: {result.stderr}")
            if check:
                raise subprocess.CalledProcessError(
                    result.returncode, self.cmd, result.stdout, result.stderr
                )

        return result


class WorkflowManager(WorkflowBase):
    """Pure Python workflow manager for VCFcache operations.

    This class implements the same workflow logic as Nextflow but using
    pure Python subprocess calls. It eliminates the Java/Nextflow dependency
    while maintaining identical functionality and output.
    """

    def __init__(
        self,
        input_file: Path,
        output_dir: Path,
        name: str,
        workflow: Path | None = None,
        anno_config_file: Optional[Path] = None,
        params_file: Optional[Path] = None,
        output_file: Optional[Path | str] = None,
        verbosity: int = 0,
    ):
        """Initialize the pure Python workflow manager.

        Args:
            workflow: Path to workflow file (for compatibility, not used)
            input_file: Path to input VCF/BCF file
            output_dir: Directory for output files
            name: Unique name for this workflow instance
            anno_config_file: Optional annotation configuration YAML file
            params_file: Required YAML parameters file
            output_file: Optional output BCF path for annotate modes (or '-' for stdout)
            verbosity: Verbosity level (0=quiet, 1=info, 2=debug)
        """
        super().__init__(
            workflow=workflow,
            input_file=input_file,
            output_dir=output_dir,
            name=name,
            anno_config_file=anno_config_file,
            params_file=params_file,
            verbosity=verbosity,
        )

        # Set up logging
        log_file = self.output_dir / "workflow.log" if self.output_dir else None
        self.logger = setup_logging(
            verbosity=verbosity,
            log_file=log_file if log_file and self.output_dir.exists() else None,
        )

        self.logger.debug(f"Initializing Pure Python workflow in: {self.output_dir}")
        if output_file is not None:
            if str(output_file) in {"-", "stdout"}:
                self.output_file = "-"
            else:
                self.output_file = Path(output_file).expanduser().resolve()
        else:
            self.output_file = None

        # Load and validate params file (required for all modes)
        if params_file:
            self.params_file = Path(params_file).expanduser().resolve()
            if not self.params_file.exists():
                self.logger.error(f"Parameters file not found: {self.params_file}")
                raise FileNotFoundError(
                    f"Parameters file not found: {self.params_file}"
                )

            # Load params with environment variable expansion
            params_content = self.params_file.read_text()
            params_expanded = os.path.expandvars(params_content)
            raw_params = yaml.safe_load(params_expanded) or {}

            # Validate params.yaml schema
            from vcfcache.utils.schemas import ParamsYAMLSchema
            is_valid, error = ParamsYAMLSchema.validate(raw_params, self.params_file)
            if not is_valid:
                self.logger.error(error)
                raise ValueError(error)

            def _expand(val):
                if isinstance(val, str):
                    return os.path.expanduser(os.path.expandvars(val))
                if isinstance(val, dict):
                    return {k: _expand(v) for k, v in val.items()}
                if isinstance(val, list):
                    return [_expand(v) for v in val]
                return val

            self.params_file_content = _expand(raw_params)

            self.logger.debug(f"Loaded parameters from: {self.params_file}")
        else:
            self.params_file_content = {}

        # Load annotation config if provided (YAML format)
        if anno_config_file:
            self.nfa_config = Path(anno_config_file).expanduser().resolve()
            if not self.nfa_config.exists():
                self.logger.error(
                    f"Annotation config file not found: {self.nfa_config}"
                )
                raise FileNotFoundError(
                    f"Annotation config file not found: {self.nfa_config}"
                )

            # Load annotation config (YAML format, not Groovy)
            anno_content = self.nfa_config.read_text()
            anno_expanded = os.path.expandvars(anno_content)
            self.nfa_config_content = yaml.safe_load(anno_expanded)

            # Validate annotation.yaml schema
            from vcfcache.utils.schemas import AnnotationYAMLSchema
            is_valid, error = AnnotationYAMLSchema.validate(self.nfa_config_content, self.nfa_config)
            if not is_valid:
                self.logger.error(error)
                raise ValueError(error)

            self.logger.debug(f"Loaded annotation config from: {self.nfa_config}")
        else:
            self.nfa_config_content = {}

        params_genome = None
        anno_genome = None
        if self.params_file_content:
            params_genome = self.params_file_content.get("genome_build")
            if params_genome:
                self.logger.info(f"Genome build (params.yaml): {params_genome}")
        if self.nfa_config_content:
            anno_genome = self.nfa_config_content.get("genome_build")
            if anno_genome:
                self.logger.info(f"Genome build (annotation.yaml): {anno_genome}")
        if params_genome and anno_genome and params_genome != anno_genome:
            raise ValueError(
                f"Genome build mismatch between params.yaml ({params_genome}) "
                f"and annotation.yaml ({anno_genome})."
            )

    def run(
        self,
        db_mode: str,
        nextflow_args: Optional[List[str]] = None,
        trace: bool = False,
        db_bcf: Optional[Path] = None,
        dag: bool = False,
        timeline: bool = False,
        report: bool = False,
        temp: Union[Path, str] = "/tmp",
        preserve_unannotated: bool = False,
        skip_split_multiallelic: bool = False,
    ) -> subprocess.CompletedProcess:
        """Execute the workflow using pure Python.

        Args:
            db_mode: Workflow mode ('blueprint-init', 'blueprint-extend', 'cache-build', 'annotate')
            nextflow_args: Additional arguments (repurposed for normalize flag)
            trace: Whether to generate trace file
            db_bcf: Path to database BCF file (required for blueprint-extend and annotate)
            dag: Ignored (not supported in pure Python mode)
            timeline: Ignored (not supported in pure Python mode)
            report: Ignored (not supported in pure Python mode)
            temp: Temporary directory for intermediate files
            preserve_unannotated: Whether to preserve variants without annotation in output
            skip_split_multiallelic: Skip splitting multiallelic variants (use only if certain input has none)

        Returns:
            subprocess.CompletedProcess with execution results

        Raises:
            ValueError: If invalid db_mode or missing required parameters
            RuntimeError: If workflow execution fails
        """
        start_time = datetime.datetime.now()
        normalize = bool(
            nextflow_args
            and any(a in ("-n", "--normalize", "normalize") for a in nextflow_args)
        )

        # Create work directory
        self._create_work_dir(self.output_dir, dirname="work")

        self.logger.info(f"Running workflow in mode: {db_mode}")

        try:
            # Route to appropriate workflow method
            if db_mode == "blueprint-init":
                result = self._run_blueprint_init(normalize=normalize)
            elif db_mode == "blueprint-extend":
                if not db_bcf:
                    raise ValueError("db_bcf is required for blueprint-extend mode")
                result = self._run_blueprint_extend(db_bcf, normalize=normalize)
            elif db_mode == "cache-build":
                if not db_bcf:
                    raise ValueError("db_bcf is required for cache-build mode")
                result = self._run_cache_build(db_bcf)
            elif db_mode == "annotate":
                if not db_bcf:
                    raise ValueError("db_bcf is required for annotate mode")
                result = self._run_annotate(db_bcf, preserve_unannotated=preserve_unannotated, skip_split_multiallelic=skip_split_multiallelic)
            elif db_mode == "annotate-nocache":
                result = self._run_annotate_nocache(skip_split_multiallelic=skip_split_multiallelic)
            else:
                raise ValueError(
                    f"Invalid db_mode: {db_mode}. Must be one of: "
                    "blueprint-init, blueprint-extend, cache-build, annotate, annotate-nocache"
                )

            end_time = datetime.datetime.now()

            # Write trace file if requested
            if trace:
                self._write_trace_file(db_mode, start_time, end_time)

            self.logger.info(
                f"Workflow completed successfully in {(end_time - start_time).total_seconds():.1f}s"
            )

            return result

        except subprocess.CalledProcessError as e:
            self.warn_temp_files()
            self.logger.error(f"Workflow execution failed: {e}")
            raise RuntimeError(f"Workflow execution failed: {e}") from e
        except Exception as e:
            self.warn_temp_files()
            self.logger.error(f"Unexpected error: {e}")
            raise

    def _run_blueprint_init(self, normalize: bool = False) -> subprocess.CompletedProcess:
        """Initialize database blueprint.

        Returns:
            subprocess.CompletedProcess with results
        """
        self.logger.info("Initializing database blueprint")

        work_task = self.work_dir / "blueprint-init"
        work_task.mkdir(parents=True, exist_ok=True)

        input_bcf = self.input_file
        output_bcf = self.output_dir / "vcfcache.bcf"
        bcftools = self.params_file_content["bcftools_cmd"]

        # Get thread count from params
        threads = self.params_file_content.get("threads", 1)

        if normalize:
            self.logger.info("Removing GT/INFO fields and splitting multiallelic sites")
            cmd = f"""{bcftools} view -G -Ou --threads {threads} {input_bcf} | \\
                    {bcftools} annotate -x INFO -Ou --threads {threads} | \\
                    {bcftools} norm -m- -o {output_bcf} -Ob --write-index --threads {threads}
                """
        else:
            self.logger.info("Removing GT/INFO fields (no multiallelic splitting)")
            cmd = f"""{bcftools} view -G -Ou --threads {threads} {input_bcf} | \\
                    {bcftools} annotate -x INFO -Ou --threads {threads} | \\
                    {bcftools} view -o {output_bcf} -Ob --write-index --threads {threads}
                """

        return BcftoolsCommand(cmd, self.logger, work_task).run()

    def _run_blueprint_extend(
        self, db_bcf: Path, normalize: bool = False
    ) -> subprocess.CompletedProcess:
        """Add variants to existing blueprint.

        Args:
            db_bcf: Path to existing blueprint BCF

        Returns:
            subprocess.CompletedProcess with results
        """
        self.logger.info("Adding variants to existing blueprint")

        work_task = self.work_dir / "blueprint-extend"
        work_task.mkdir(parents=True, exist_ok=True)

        bcftools = self.params_file_content["bcftools_cmd"]

        # Get thread count from params
        threads = self.params_file_content.get("threads", 1)

        # Step 1: Filter new input (drop INFO), optionally split multiallelics
        filtered = work_task / "filtered.bcf"
        input_bcf = self.input_file

        if normalize:
            self.logger.info("Filtering new input (drop INFO) and splitting multiallelic sites")
            step_cmd = f"""{bcftools} view -G -Ou --threads {threads} {input_bcf} | \\
                {bcftools} annotate -x INFO -Ou --threads {threads} | \\
                {bcftools} norm -m- -o {filtered} -Ob --write-index --threads {threads}
            """
        else:
            self.logger.info("Filtering new input (drop INFO) (no multiallelic splitting)")
            step_cmd = f"""{bcftools} view -G -Ou --threads {threads} {input_bcf} | \\
                {bcftools} annotate -x INFO -Ou --threads {threads} | \\
                {bcftools} view -o {filtered} -Ob --write-index --threads {threads}
            """

        BcftoolsCommand(step_cmd, self.logger, work_task).run()

        # Step 2: Merge with existing blueprint
        output_bcf = self.output_dir / "vcfcache.bcf"
        self.logger.info(f"Merging with existing blueprint: {db_bcf}")

        merge_cmd = (
            f"{bcftools} merge -m none --threads {threads} {db_bcf} {filtered} -o {output_bcf} -Ob --write-index"
        )

        return BcftoolsCommand(merge_cmd, self.logger, work_task).run()

    def _run_cache_build(self, db_bcf: Path) -> subprocess.CompletedProcess:
        """Run annotation on blueprint to create cache.

        Args:
            db_bcf: Path to blueprint BCF

        Returns:
            subprocess.CompletedProcess with results
        """
        self.logger.info("Annotating blueprint to create cache")

        work_task = self.work_dir / "cache-build"
        work_task.mkdir(parents=True, exist_ok=True)

        # Prepare auxiliary directory
        aux_dir = self.output_dir / "auxiliary"
        aux_dir.mkdir(exist_ok=True)

        output_bcf = self.output_dir / "vcfcache_annotated.bcf"

        # Substitute variables in annotation command
        anno_cmd = self._substitute_variables(
            self.nfa_config_content["annotation_cmd"],
            extra_vars={
                "INPUT_BCF": str(db_bcf),
                "OUTPUT_BCF": str(output_bcf),
                "AUXILIARY_DIR": str(aux_dir),
            },
        )

        self.logger.info("Running annotation command on blueprint")

        # Execute annotation command
        result = BcftoolsCommand(anno_cmd, self.logger, work_task).run()

        # Copy annotation logs to cache directory for posterity
        # These logs contain the full output from the annotation tool (VEP, SnpEff, etc.)
        stdout_file = work_task / "stdout.txt"
        stderr_file = work_task / "stderr.txt"

        # Save stdout to annotation_tool.log
        annotation_log = self.output_dir / "annotation_tool.log"
        if stdout_file.exists() and stdout_file.stat().st_size > 0:
            annotation_log.write_text(stdout_file.read_text())
            self.logger.info(f"Annotation output saved to: {annotation_log}")

        # Save stderr to annotation_tool_err.log (only if non-empty)
        annotation_err_log = self.output_dir / "annotation_tool_err.log"
        if stderr_file.exists() and stderr_file.stat().st_size > 0:
            annotation_err_log.write_text(stderr_file.read_text())
            self.logger.info(f"Annotation errors saved to: {annotation_err_log}")

        # Remove auxiliary directory if empty (some annotation tools don't use it)
        if aux_dir.exists() and not any(aux_dir.iterdir()):
            aux_dir.rmdir()
            self.logger.debug("Removed empty auxiliary directory")

        # Validate output has required INFO tag
        tag = self.nfa_config_content["must_contain_info_tag"]
        self._validate_info_tag(output_bcf, tag)

        self.logger.info(f"Annotation complete, cache created at: {output_bcf}")

        return result

    def _resolve_output_bcf(self, sample_name: str) -> Path | str:
        """Resolve output target for annotate modes."""
        if self.output_file is None:
            return self.output_dir / f"{sample_name}_vc.bcf"
        return self.output_file

    def _run_annotate(self, db_bcf: Path, preserve_unannotated: bool = False, skip_split_multiallelic: bool = False) -> subprocess.CompletedProcess:
        """Annotate sample using cache - 4-step caching process.

        This is the key performance feature that makes VCFcache fast.
        It uses pre-annotated common variants from the cache and only
        annotates novel variants.

        Args:
            db_bcf: Path to cache BCF (vcfcache_annotated.bcf)
            preserve_unannotated: Whether to preserve variants without annotation in output
            skip_split_multiallelic: Skip splitting multiallelic variants (use only if certain input has none)

        Returns:
            subprocess.CompletedProcess with results

        Example and Debug:
          ROOT = Path("/home/j380r/projects/vcfcache").resolve()  # or Path.cwd() if you're in repo root
          TEST_DATA = ROOT / "tests" / "data" / "nodata"

          sample_bcf   = TEST_DATA / "sample4.bcf"       # toy sample to annotate
          blueprint_bcf = TEST_DATA / "gnomad_test.bcf"  # toy “blueprint” to turn into an annotated cache

          params = ROOT / "tests" / "config" / "test_params.yaml"
          anno   = ROOT / "tests" / "config" / "test_annotation.yaml"  # adds INFO/MOCK_ANNO

          scratch = Path(tempfile.mkdtemp(prefix="vcfcache_dbg_"))
          cache_dir = scratch / "cache"
          out_dir   = scratch / "annot_out"
          cache_dir.mkdir(parents=True, exist_ok=True)
          out_dir.mkdir(parents=True, exist_ok=True)

          # 1) Build a tiny annotated cache BCF (creates cache_dir/vcfcache_annotated.bcf in tmp)
          self = WorkflowManager(
              input_file=sample_bcf,
              output_dir=cache_dir,
              name="dbg_cache",
              anno_config_file=anno,
              params_file=params,
              verbosity=2,
          )
                  # Create work directory
          self._create_work_dir(self.output_dir, dirname="work")
          db_bcf=blueprint_bcf
        """
        self.logger.info("Annotating sample using cache (4-step process)")

        work_task = self.work_dir / "annotate"
        work_task.mkdir(parents=True, exist_ok=True)

        sample_name = self.input_file.stem
        input_bcf = self.input_file
        bcftools = self.params_file_content["bcftools_cmd"]
        tag = self.nfa_config_content["must_contain_info_tag"]

        # Get thread count from params
        threads = self.params_file_content.get("threads", 1)

        # Step 0: Split multiallelic variants and remove spanning deletions (ALT=*)
        if skip_split_multiallelic:
            self.logger.info("Step 0/4: Skipping multiallelic variant splitting (--skip-split-multiallelic flag set)")
            # Still need to remove spanning deletions even when skipping split
            self.logger.info("Step 0/4: Removing spanning deletion alleles (ALT=*)")
            filtered_input = work_task / f"{sample_name}_no_span_del.bcf"
            filter_cmd = (
                f"{bcftools} view -e 'ALT=\"*\"' {input_bcf} -o {filtered_input} -Ob -W --threads {threads}"
            )
            BcftoolsCommand(filter_cmd, self.logger, work_task).run()
            input_bcf = filtered_input
        else:
            self.logger.info("Step 0/4: Splitting multiallelic variants and removing spanning deletions")
            normalized_input = work_task / f"{sample_name}_normalized.bcf"
            # Combine norm and filter: split multiallelic, then remove ALT=*
            norm_cmd = (
                f"{bcftools} norm -m- {input_bcf} -Ob | "
                f"{bcftools} view -e 'ALT=\"*\"' -o {normalized_input} -Ob -W --threads {threads}"
            )
            BcftoolsCommand(norm_cmd, self.logger, work_task).run()
            input_bcf = normalized_input  # Use normalized input for subsequent steps

        # Step 1: Add cache annotations
        self.logger.info("Step 1/4: Adding cache annotations")
        step1_bcf = work_task / f"{sample_name}_isecvst.bcf"
        cmd1 = (
            f"{bcftools} annotate -a {db_bcf} {input_bcf} -c INFO "
            f"-o {step1_bcf} -Ob -W --threads {threads}"
        )
        BcftoolsCommand(cmd1, self.logger, work_task).run()

        # Step 2: Filter to get missing annotations
        self.logger.info("Step 2/4: Identifying variants missing from cache")
        step2_bcf = work_task / f"{sample_name}_isecvst_miss.bcf"
        cmd2 = (
            f"{bcftools} filter -i 'INFO/{tag}==\"\"' {step1_bcf} "
            f"-o {step2_bcf} -Ob -W --threads {threads}"
        )
        BcftoolsCommand(cmd2, self.logger, work_task).run()

        # Count missing variants
        count_result = subprocess.run(
            f"{bcftools} index -n {step2_bcf}.csi",
            shell=True,
            capture_output=True,
            text=True,
        )
        missing_count = (
            int(count_result.stdout.strip()) if count_result.returncode == 0 else 0
        )
        self.logger.info(f"Found {missing_count} variants not in cache")

        # Step 3: Annotate only missing variants
        if missing_count > 0:
            self.logger.info(f"Step 3/4: Annotating {missing_count} missing variants")
            step3_bcf = work_task / f"{sample_name}_missing_annotated.bcf"
            aux_dir = self.output_dir / "auxiliary"
            aux_dir.mkdir(exist_ok=True)

            anno_cmd = self._substitute_variables(
                self.nfa_config_content["annotation_cmd"],
                extra_vars={
                    "INPUT_BCF": str(step2_bcf),
                    "OUTPUT_BCF": str(step3_bcf),
                    "AUXILIARY_DIR": str(aux_dir),
                },
            )
            BcftoolsCommand(anno_cmd, self.logger, work_task).run()

            # Save annotation tool logs for novel variants
            stdout_file = work_task / "stdout.txt"
            stderr_file = work_task / "stderr.txt"

            annotation_log = self.output_dir / "annotation_tool.log"
            if stdout_file.exists() and stdout_file.stat().st_size > 0:
                annotation_log.write_text(stdout_file.read_text())
                self.logger.info(f"Annotation output saved to: {annotation_log}")

            annotation_err_log = self.output_dir / "annotation_tool_err.log"
            if stderr_file.exists() and stderr_file.stat().st_size > 0:
                annotation_err_log.write_text(stderr_file.read_text())
                self.logger.info(f"Annotation errors saved to: {annotation_err_log}")

            # Remove auxiliary directory if empty
            if aux_dir.exists() and not any(aux_dir.iterdir()):
                aux_dir.rmdir()
                self.logger.debug("Removed empty auxiliary directory")
        else:
            self.logger.info("Step 3/4: Skipped (all variants found in cache)")
            # Create empty file for step 4
            step3_bcf = step2_bcf

        # Step 4: Merge newly annotated back into original
        self.logger.info("Step 4/4: Merging cache and newly annotated variants")

        # First merge annotations into a temporary file
        step4_bcf = work_task / f"{sample_name}_merged.bcf"
        output_bcf = self._resolve_output_bcf(sample_name)
        output_target = output_bcf
        if output_bcf == "-":
            output_target = work_task / f"{sample_name}_vc.bcf"

        if missing_count > 0:
            # Check if annotation tool dropped any variants
            step3_count_result = subprocess.run(
                f"{bcftools} index -n {step3_bcf}.csi",
                shell=True,
                capture_output=True,
                text=True,
            )
            step3_count = (
                int(step3_count_result.stdout.strip()) if step3_count_result.returncode == 0 else 0
            )

            dropped_count = 0
            if step3_count < missing_count:
                dropped_count = missing_count - step3_count
                self.logger.info(
                    f"Annotation tool dropped {dropped_count} variants from input. "
                    f"By default, these are also removed from cached output to match annotation tool behavior."
                )

            # Merge annotations from step3 into step1
            cmd4 = (
                f"{bcftools} annotate -a {step3_bcf} {step1_bcf} -c INFO "
                f"-o {step4_bcf} -Ob -W --threads {threads}"
            )
        else:
            # No new annotations, just copy step1
            cmd4 = f"cp {step1_bcf} {step4_bcf} && cp {step1_bcf}.csi {step4_bcf}.csi"
            step3_count = 0
            dropped_count = 0

        BcftoolsCommand(cmd4, self.logger, work_task).run()

        # Post-filter: Remove variants without annotation to match annotation tool behavior
        # This ensures cached output is identical to uncached output
        # Skip filtering if --preserve-unannotated flag is set
        if preserve_unannotated:
            self.logger.info(
                "Preserving unannotated variants in output (--preserve-unannotated flag set)"
            )
            if output_bcf == "-":
                copy_cmd = f"{bcftools} view {step4_bcf} -Ob -o -"
                result = BcftoolsCommand(copy_cmd, self.logger, work_task).run()
            else:
                copy_cmd = f"cp {step4_bcf} {output_bcf} && cp {step4_bcf}.csi {output_bcf}.csi"
                result = BcftoolsCommand(copy_cmd, self.logger, work_task).run()
        else:
            if output_bcf == "-":
                filter_cmd = (
                    f"{bcftools} view -i 'INFO/{tag}!=\"\"' {step4_bcf} "
                    f"-o - -Ob --threads {threads}"
                )
            else:
                filter_cmd = (
                    f"{bcftools} view -i 'INFO/{tag}!=\"\"' {step4_bcf} "
                    f"-o {output_bcf} -Ob -W --threads {threads}"
                )
            result = BcftoolsCommand(filter_cmd, self.logger, work_task).run()

        self.last_run_stats = {
            "missing_variants": missing_count,
            "missing_annotated": step3_count,
            "dropped_variants": dropped_count if missing_count > 0 else 0,
        }
        self.logger.info(f"Annotation complete: {output_bcf}")

        return result

    def _run_annotate_nocache(self, skip_split_multiallelic: bool = False) -> subprocess.CompletedProcess:
        """Annotate sample directly without using cache.

        This is for benchmarking or when cache is not available.

        Args:
            skip_split_multiallelic: Skip splitting multiallelic variants (use only if certain input has none)

        Returns:
            subprocess.CompletedProcess with results
        """
        self.logger.info("Annotating sample directly (no cache)")

        work_task = self.work_dir / "annotate-nocache"
        work_task.mkdir(parents=True, exist_ok=True)

        sample_name = self.input_file.stem
        input_bcf = self.input_file
        output_bcf = self._resolve_output_bcf(sample_name)
        output_target = output_bcf
        if output_bcf == "-":
            output_target = work_task / f"{sample_name}_vc.bcf"
        aux_dir = self.output_dir / "auxiliary"
        aux_dir.mkdir(exist_ok=True)
        bcftools = self.params_file_content["bcftools_cmd"]

        # Get thread count from params
        threads = self.params_file_content.get("threads", 1)

        # Step 0: Split multiallelic variants and remove spanning deletions (ALT=*)
        if skip_split_multiallelic:
            self.logger.info("Skipping multiallelic variant splitting (--skip-split-multiallelic flag set)")
            # Still need to remove spanning deletions even when skipping split
            self.logger.info("Removing spanning deletion alleles (ALT=*)")
            filtered_input = work_task / f"{sample_name}_no_span_del.bcf"
            filter_cmd = (
                f"{bcftools} view -e 'ALT=\"*\"' {input_bcf} -o {filtered_input} -Ob -W --threads {threads}"
            )
            BcftoolsCommand(filter_cmd, self.logger, work_task).run()
            input_bcf = filtered_input
        else:
            self.logger.info("Splitting multiallelic variants and removing spanning deletions")
            normalized_input = work_task / f"{sample_name}_normalized.bcf"
            # Combine norm and filter: split multiallelic, then remove ALT=*
            norm_cmd = (
                f"{bcftools} norm -m- {input_bcf} -Ob | "
                f"{bcftools} view -e 'ALT=\"*\"' -o {normalized_input} -Ob -W --threads {threads}"
            )
            BcftoolsCommand(norm_cmd, self.logger, work_task).run()
            input_bcf = normalized_input  # Use normalized input for annotation

        # Use a unique placeholder for stdin when piping
        STDIN_PLACEHOLDER = "__VCFCACHE_STDIN__"

        anno_cmd = self._substitute_variables(
            self.nfa_config_content["annotation_cmd"],
            extra_vars={
                "INPUT_BCF": str(input_bcf),
                "OUTPUT_BCF": str(output_target),
                "AUXILIARY_DIR": str(aux_dir),
            },
        )

        result = BcftoolsCommand(anno_cmd, self.logger, work_task).run()

        if output_bcf == "-":
            stream_cmd = f"{bcftools} view {output_target} -Ob -o -"
            result = BcftoolsCommand(stream_cmd, self.logger, work_task).run()

        # Validate output
        tag = self.nfa_config_content["must_contain_info_tag"]
        self._validate_info_tag(output_target, tag)

        input_count = self._count_variants(input_bcf, bcftools)
        output_count = self._count_variants(output_target, bcftools)
        dropped = 0
        if input_count is not None and output_count is not None:
            dropped = max(input_count - output_count, 0)
        self.last_run_stats = {
            "input_variants": input_count,
            "output_variants": output_count,
            "dropped_variants": dropped,
        }

        self.logger.info(f"Direct annotation complete: {output_bcf}")

        return result

    @staticmethod
    def _count_variants(path: Path, bcftools_cmd: str) -> Optional[int]:
        try:
            result = subprocess.run(
                f"{bcftools_cmd} index -n {path}",
                shell=True,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                return None
            return int(result.stdout.strip())
        except Exception:
            return None

    def _substitute_variables(
        self, text: str, extra_vars: Optional[Dict[str, str]] = None, skip_vars: Optional[list] = None
    ) -> str:
        """Replace variables in command strings.

        Order of substitution:
        1. Environment variables (${VCFCACHE_ROOT})
        2. Params file variables (${params.vep_cache})
        3. Special workflow variables (${INPUT_BCF}, ${OUTPUT_BCF}, ${AUXILIARY_DIR})

        Supports both escaped (\\${VAR}) and unescaped (${VAR}, $VAR) variable formats.
        Escaped variables are useful in YAML to prevent premature expansion.

        Args:
            text: Text with variables to substitute
            extra_vars: Additional variables to substitute (e.g., INPUT_BCF, OUTPUT_BCF)
            skip_vars: List of variable names to skip substitution for

        Returns:
            Text with variables substituted
        """
        skip_vars = skip_vars or []

        # 1. Environment variables
        text = os.path.expandvars(text)

        # 2. Params file variables
        if self.params_file_content:
            for key, value in self.params_file_content.items():
                if key != "optional_checks":
                    text = text.replace(f"${{params.{key}}}", str(value))

        # 3. Special workflow variables
        if extra_vars:
            for key, value in extra_vars.items():
                if key not in skip_vars:
                    # Handle escaped variables first (e.g., \${INPUT_BCF} from YAML)
                    text = text.replace(f"\\${{{key}}}", str(value))
                    # Then handle normal variables
                    text = text.replace(f"${{{key}}}", str(value))
                    text = text.replace(f"${key}", str(value))  # Also support $VAR format

        return text

    def _validate_info_tag(self, bcf_path: Path, tag: str) -> None:
        """Validate that BCF header contains required INFO tag.

        Args:
            bcf_path: Path to BCF file to validate
            tag: INFO tag that must be present

        Raises:
            RuntimeError: If tag is not found in header
        """
        bcftools = self.params_file_content["bcftools_cmd"]

        result = subprocess.run(
            f"{bcftools} view -h {bcf_path} | grep '##INFO=<ID={tag},'",
            shell=True,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"Annotation validation failed: Required INFO tag '{tag}' not found in {bcf_path}. "
                f"Check that your annotation command is working correctly."
            )

        self.logger.debug(f"Validated INFO tag '{tag}' present in {bcf_path}")

    def _write_trace_file(self, mode: str, start_time, end_time) -> None:
        """Write execution trace file.

        Args:
            mode: Workflow mode
            start_time: Start time
            end_time: End time
        """
        trace_file = self.output_dir / f"{self.name}_trace.txt"

        duration = (end_time - start_time).total_seconds()

        with open(trace_file, "w") as f:
            f.write(f"task_id\tname\tstatus\texit\tduration\n")
            f.write(f"1\t{mode}\tCOMPLETED\t0\t{duration:.1f}s\n")

        self.logger.debug(f"Trace file written to: {trace_file}")

    def _create_work_dir(self, parent: Path, dirname: str = "work") -> None:
        """Create a temporary work directory.

        Args:
            parent: Parent directory
            dirname: Name of work directory
        """
        self.work_dir = parent / dirname
        if self.work_dir is not None:
            if not self.work_dir.exists():
                self.work_dir.mkdir(parents=True, exist_ok=True)
                if self.logger:
                    self.logger.debug(f"Created work directory: {self.work_dir}")
            else:
                if self.logger:
                    self.logger.warning(
                        f"Work directory already exists: {self.work_dir}"
                    )

    def cleanup_work_dir(self) -> None:
        """Remove temporary work directory."""
        if not self.work_dir:
            return

        self.logger.debug("Cleaning up work directory")
        try:
            if self.work_dir.exists():
                self.logger.debug(f"Removing work directory: {self.work_dir}")
                shutil.rmtree(self.work_dir)
        except Exception as e:
            self.logger.warning(f"Failed to remove work directory {self.work_dir}: {e}")
        self.work_dir = None
