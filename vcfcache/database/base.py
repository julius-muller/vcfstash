# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2024-2026 Julius Müller

import importlib.resources
import json
import os
import re
import shutil
import subprocess
import tempfile
from logging import Logger
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pysam
import yaml

from vcfcache.database.outputs import CacheOutput
from vcfcache.database.workflow_base import WorkflowBase
from vcfcache.utils.logging import setup_logging
from vcfcache.utils.validation import validate_bcf_header


def create_workflow(**kwargs) -> WorkflowBase:
    """Factory function to create workflow backend.

    This function creates a WorkflowManager (pure Python) instance for running
    vcfcache workflows. The pure Python backend replaces the legacy Nextflow
    system and provides a simpler, faster workflow execution.

    Args:
        **kwargs: Arguments to pass to the workflow constructor:
                 - workflow: Path to workflow file (ignored, for compatibility)
                 - input_file: Path to input VCF/BCF
                 - output_dir: Output directory
                 - name: Workflow instance name
                 - anno_config_file: Optional annotation config
                 - params_file: Optional YAML params
                 - output_file: Optional output BCF path for annotate modes
                 - verbosity: Verbosity level

    Returns:
        WorkflowBase instance (WorkflowManager)

    Example:
        >>> workflow = create_workflow(
        ...     input_file=Path("input.bcf"),
        ...     output_dir=Path("/tmp/output"),
        ...     name="test",
        ...     anno_config_file=Path("annotation.yaml"),
        ...     params_file=Path("params.yaml"),
        ...     verbosity=1
        ... )
    """
    from vcfcache.database.workflow_manager import WorkflowManager

    # Ensure a placeholder workflow path exists for compatibility
    if "workflow" not in kwargs or kwargs.get("workflow") is None:
        output_dir = kwargs.get("output_dir")
        if output_dir is None:
            raise ValueError("output_dir is required to initialize workflow backend")
        kwargs["workflow"] = Path(output_dir) / "workflow.stub"

    return WorkflowManager(**kwargs)


class VCFDatabase:
    """Base class for performing operations on VCF (Variant Call Format) databases.

    Provides tools for managing, validating, and manipulating VCF/BCF files including
    index validation, header validation, INFO field parsing, and genome reference checks.
    Also includes utilities for logging and configuration setup designed for bioinformatics
    workflow management.

    Attributes:
        TRANSCRIPT_KEYS: List of predefined transcript consequence keys used in
                         INFO annotations of VCF files.
    """

    TRANSCRIPT_KEYS = [
        "SYMBOL",
        "Feature",
        "Consequence",
        "HGVS_OFFSET",
        "HGVSc",
        "HGVSp",
        "IMPACT",
        "DISTANCE",
        "PICK",
        "VARIANT_CLASS",
    ]
    """Base class for VCF database operations"""

    def __init__(self, db_path: Path, verbosity: int, debug: bool, bcftools_path: Path):
        self.debug = debug
        self.cached_output = CacheOutput(str(db_path))
        self.cache_path = self.cached_output.root_dir
        self.cache_name = self.cache_path.name
        self.blueprint_dir = self.cache_path / "blueprint"
        self.workflow_dir = self.cache_path / "workflow"
        self.cache_dir = self.cache_path / "cache"
        self.workflow_dir_src = self.cached_output.workflow_src_dir
        self.blueprint_bcf = self.blueprint_dir / "vcfcache.bcf"
        self.info_file = self.blueprint_dir / "sources.info"
        self.verbosity = verbosity if not debug else 2
        self.bcftools_path = bcftools_path

        if self.info_file.exists():
            with open(self.info_file) as f:
                self.db_info = json.load(f)
        else:
            self.db_info = {"input_files": []}
        self.logger: Optional[Logger] = None

    def connect_loggers(self, logger_name: str = "vcfdb") -> Logger:
        # Set up central logging
        log_file = self.cache_path / f"{logger_name}.log"
        return setup_logging(verbosity=self.verbosity, log_file=log_file)

    def ensure_indexed(self, file_path: Path) -> None:
        """Ensure input file has an index file (CSI or TBI)"""
        if self.logger:
            self.logger.debug(f"Checking index for: {file_path}")
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        csi_index = Path(str(file_path) + ".csi")
        tbi_index = Path(str(file_path) + ".tbi")

        if not (csi_index.exists() or tbi_index.exists()):
            raise RuntimeError(
                f"No index found for {file_path}. Use bcftools index for BCF/compressed VCF "
                "or tabix for compressed VCF files."
            )

    def _list_contigs(self, bcf: Path) -> list[str]:
        """List contigs using bcftools index -s."""
        self.ensure_indexed(bcf)
        res = subprocess.run(
            [str(self.bcftools_path), "index", "-s", str(bcf)],
            check=True,
            capture_output=True,
            text=True,
        )
        return [line.split()[0] for line in res.stdout.splitlines() if line.strip()]

    def validate_bcf_header(
        self, bcf_path: Path, norm: bool = True
    ) -> Tuple[bool, Optional[str]]:
        """Validate BCF header format"""
        if self.logger:
            self.logger.debug(f"Validating BCF header: {bcf_path}")
        result = validate_bcf_header(bcf_path, norm, bcftools_path=self.bcftools_path)
        if not result[0] and self.logger:
            self.logger.error(f"BCF header validation failed: {result[1]}")
        return result

    def parse_vcf_info(self, vcf_data: list) -> list:
        """Parses VCF INFO field and expands transcript consequences."""

        def convert_vcfstr(value):
            if value is None or value == "":
                return None
            try:
                if isinstance(value, str) and (
                    "." in value or "e" in value or "E" in value
                ):
                    return float(value)
                return int(value)
            except ValueError:
                return value

        expanded_data: List[Dict[str, Any]] = []
        for tn, transcript in enumerate(vcf_data):
            expanded_data.append({})
            for entry in self.TRANSCRIPT_KEYS:
                if entry not in vcf_data[tn]:
                    raise ValueError(
                        f"Did not find key={entry} {transcript=} in CSQ INFO tag"
                    )
                if entry == "PICK":
                    expanded_data[tn][entry] = vcf_data[tn][entry] == "1"
                else:
                    expanded_data[tn][entry] = convert_vcfstr(vcf_data[tn][entry])

        return expanded_data

    def validate_vcf_reference(
        self,
        vcf_path: Path,
        ref_fasta: Path,
        variants_to_check: int = 3,
    ) -> Tuple[bool, Optional[str]]:
        r"""Validates a VCF/BCF file against a reference genome with automatic chromosome name mapping.

        This method performs two checks:
        1. Verifies that all contigs in the VCF/BCF index are found in the reference FASTA index
           (with automatic chr prefix and mitochondrial chromosome handling)
        2. Checks a configurable number of variants from different chromosomes to ensure their reference
           alleles match the reference genome

        Args:
            vcf_path (Path): Path to the VCF/BCF file to validate
            ref_fasta (Path): Path to the reference genome FASTA file
            variants_to_check (int, optional): Number of variants to check. Defaults to 3.

        Returns:
            Tuple[bool, Optional[str]]: A tuple containing:
                - bool: True if validation passed, False otherwise
                - Optional[str]: Error message if validation failed, None otherwise
        """
        if self.logger:
            self.logger.debug(f"Validating VCF/BCF file against reference: {vcf_path}")

        # Check if files exist
        if not vcf_path.exists():
            return False, f"VCF/BCF file not found: {vcf_path}"

        if not ref_fasta.exists():
            return False, f"Reference genome not found: {ref_fasta}"

        # Check for index files
        vcf_index_csi = Path(str(vcf_path) + ".csi")
        vcf_index_tbi = Path(str(vcf_path) + ".tbi")
        ref_index = Path(str(ref_fasta) + ".fai")

        if not (vcf_index_csi.exists() or vcf_index_tbi.exists()):
            return (
                False,
                f"VCF/BCF index not found for {vcf_path}. Use bcftools index or tabix.",
            )

        if not ref_index.exists():
            return (
                False,
                f"Reference index not found for {ref_fasta}. Use samtools faidx.",
            )

        try:
            # Load reference contigs from .fai
            ref_contigs = set()
            with open(ref_index, "r") as f:
                for line in f:
                    contig = line.split("\t")[0]
                    ref_contigs.add(contig)

            if self.logger:
                self.logger.debug(
                    f"Found {len(ref_contigs)} contigs in reference index"
                )

            # Use bcftools to extract chromosomes from the index
            vcf_contigs = set()
            vcf_index_file = vcf_index_csi if vcf_index_csi.exists() else vcf_index_tbi

            # Use bcftools to extract chromosomes from the index
            cmd = [self.bcftools_path, "index", "--stats", str(vcf_index_file)]
            result = subprocess.run(cmd, capture_output=True, text=True)

            # Parse chromosomes from output
            for line in result.stdout.splitlines():
                if "[" not in line and "]" not in line:  # Skip header/footer lines
                    parts = line.strip().split()
                    if parts and parts[0]:  # Ensure we have a chromosome name
                        vcf_contigs.add(parts[0])

            if self.logger:
                self.logger.debug(
                    f"Found {len(vcf_contigs)} unique chromosomes in VCF index"
                )

            # Check if all VCF contigs are in reference (exact match)
            missing_contigs = {c for c in vcf_contigs if c not in ref_contigs}

            if missing_contigs:
                return (
                    False,
                    f"Contigs in VCF not found in reference: {', '.join(missing_contigs)}",
                )

            # Check reference alleles for a subset of variants
            checked_variants = 0
            checked_chroms = set()
            mismatches = []

            # Open reference FASTA
            ref_fasta_file = pysam.FastaFile(str(ref_fasta))

            # Open VCF file
            vcf = pysam.VariantFile(str(vcf_path))

            # Iterate through variants
            for variant in vcf.fetch():
                # Skip if we've already checked this chromosome
                if variant.chrom in checked_chroms:
                    continue

                # Add chromosome to checked set
                checked_chroms.add(variant.chrom)

                # Get reference sequence - try direct match first
                try:
                    chr_to_use = variant.chrom

                    # If direct access fails, try mapped version
                    if chr_to_use not in ref_fasta_file.references:
                        if (
                            chr_to_use in chrom_map
                            and chrom_map[chr_to_use] in ref_fasta_file.references
                        ):
                            chr_to_use = chrom_map[chr_to_use]
                        # Try with/without 'chr' prefix
                        elif (
                            chr_to_use.startswith("chr")
                            and chr_to_use[3:] in ref_fasta_file.references
                        ):
                            chr_to_use = chr_to_use[3:]
                        elif f"chr{chr_to_use}" in ref_fasta_file.references:
                            chr_to_use = f"chr{chr_to_use}"

                    # Skip if variant.ref is None
                    if variant.ref is None:
                        continue

                    ref_seq = ref_fasta_file.fetch(
                        chr_to_use, variant.pos - 1, variant.pos - 1 + len(variant.ref)
                    )

                    # Compare reference allele (ensuring both are not None)
                    if ref_seq is not None and variant.ref.upper() != ref_seq.upper():
                        mismatches.append(
                            f"{variant.chrom}:{variant.pos} (VCF: {variant.ref}, FASTA: {ref_seq})"
                        )

                    checked_variants += 1

                    # Stop if we've checked enough variants
                    if checked_variants >= variants_to_check:
                        break

                except Exception as e:
                    if self.logger:
                        self.logger.warning(
                            f"Error checking variant at {variant.chrom}:{variant.pos}: {e}"
                        )

            # Check if we found any variants to check
            if checked_variants == 0:
                return False, "Could not find any variants to check reference alleles"

            # Report mismatches
            if mismatches:
                return (
                    False,
                    f"Reference allele mismatches found: {', '.join(mismatches)}",
                )

            if self.logger:
                self.logger.info(
                    f"Successfully validated {checked_variants} variants against reference"
                )
            return True, None

        except Exception as e:
            return False, f"Error validating VCF against reference: {e}"

    @staticmethod
    def _copy_workflow_srcfiles(
        source: Path, destination: Path, skip_config: bool = False
    ) -> None:
        """Copy workflow files from source to destination,

        In this case, we just ensure the destination directory exists.

        Parameters:
            destination (Path): Destination directory for copying.
        """
        try:
            # Ensure the source directory exists
            if not source.exists():
                # Source doesn't exist - just create destination
                destination.mkdir(parents=True, exist_ok=True)
                return

            # Ensure the destination directory exists
            destination.mkdir(parents=True, exist_ok=True)

            # Check if source has any files to copy
            source_files = list(source.glob("*"))
            if not source_files:
                return

            # Copy entire directory structure except *.config files
            if skip_config:
                shutil.copytree(
                    source,
                    destination,
                    ignore=shutil.ignore_patterns("*.config"),
                    dirs_exist_ok=True,
                )
            else:
                shutil.copytree(source, destination, dirs_exist_ok=True)

        except Exception as e:
            raise RuntimeError(f"Error copying workflow files: {str(e)}") from e
