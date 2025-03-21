import json
import shutil
import sys
import tarfile
import tempfile
import time
from src.utils.logging import setup_logger
from typing import List, Optional, Dict
from ..utils.validation import compute_md5
from pathlib import Path
import subprocess
from datetime import datetime
from multiprocessing import Pool
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

from .base import VEPDatabase, NextflowWorkflow
from ..utils.validation import validate_vcf_format, get_bcf_stats

# Add at the top of annotator.py
INFO_FIELDS = ['GT', 'DP', 'AF', "gnomadg_af", "gnomade_af", "gnomadg_ac", "gnomade_ac", 'clinvar_clnsig', "deeprvat_score"]
TRANSCRIPT_KEYS = [
    'SYMBOL', 'Feature', 'Consequence', 'HGVS_OFFSET', 'HGVSc', 'HGVSp',
    'IMPACT', 'DISTANCE', 'PICK', 'VARIANT_CLASS'
]
BASES = {"A", "C", "G", "T"}

def wavg(f1: float | None, f2: float | None, n1: int, n2: int) -> float | None:
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
    elif f2 is None:
        return f1

def parse_vep_info(vep_data: list) -> list:
    """Parses VEP INFO field and expands transcript consequences."""
    def convert_vepstr(value):
        if value is None or value == "":
            return None
        try:
            if isinstance(value, str) and ("." in value or "e" in value or "E" in value):
                return float(value)
            return int(value)
        except ValueError:
            return value

    expanded_data = []
    for tn, transcript in enumerate(vep_data):
        expanded_data.append({})
        for entry in TRANSCRIPT_KEYS:
            if entry not in vep_data[tn]:
                raise ValueError(f"Did not find key={entry} in CSQ INFO tag")
            if entry == "PICK":
                expanded_data[tn][entry] = vep_data[tn][entry] == "1"
            else:
                expanded_data[tn][entry] = convert_vepstr(vep_data[tn][entry])

    return expanded_data


class DatabaseAnnotator(VEPDatabase):
    """Handles database annotation workflows"""
    def __init__(self, db_path: Path, workflow_dir: Path, params_file: Optional[Path], nextflow_args: List[str]):
        super().__init__(db_path)

        self._nf_workflow = NextflowWorkflow(workflow_dir, params_file=params_file)
        # raise ValueError('here!!')
        self.workflow_dir = self._nf_workflow.workflow_dir
        self.workflow_path = self._nf_workflow.workflow_path
        self.workflow_config_path = self._nf_workflow.workflow_config_path
        self.params_file = self._nf_workflow.params_file
        self.nextflow_args = nextflow_args
        self.run_dir = self._create_run_dir(db_dir=self.db_path , workflow_hash=self._nf_workflow.workflow_hash)
        self.run_info = self.run_dir / "annotation.info"

        # Set up logging
        log_file = self.run_dir / "annotation.log"
        self.logger = setup_logger(log_file)
        self.logger.info("Initializing database annotation")
        self.logger.debug(f"Workflow dir: {workflow_dir}")
        self.logger.debug(f"Params file: {params_file}")

    def annotate(self) -> None:
        """Run annotation workflow on database"""
        if not self.variants_bcf.exists():
            self.logger.error("Database BCF file does not exist")
            raise FileNotFoundError("Database BCF file does not exist.")

        # Copy workflow files to database
        workflow_dir = Path(self.workflow_dir)
        db_workflow_dir = Path(self.db_path).parent / "workflow"
        db_workflow_dir.mkdir(exist_ok=True)

        # Copy main files
        shutil.copy2(workflow_dir / "main.nf", db_workflow_dir / "main.nf")
        shutil.copy2(workflow_dir / "nextflow.config", db_workflow_dir / "nextflow.config")

        # Copy module files
        modules_dir = workflow_dir / "modules"
        if modules_dir.exists():
            db_modules_dir = db_workflow_dir / "modules"
            db_modules_dir.mkdir(exist_ok=True)
            for module_file in modules_dir.glob("*.nf"):
                rel_path = module_file.relative_to(workflow_dir)
                target_path = db_workflow_dir / rel_path
                target_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(module_file, target_path)
                self.logger.debug(f"Copied module file: {rel_path}")

        if self.params_file:
            shutil.copy2(self.params_file, db_workflow_dir / self.params_file.name)

        # Store blueprint snapshot and workflow files
        shutil.copy2(self.info_file, self.run_dir / "blueprint.snapshot")
        workflow_files = self._store_workflow_files(self.workflow_dir, self.run_dir)

        try:
            self.logger.info("Starting annotation workflow")

            start_time = datetime.now()
            result = self._nf_workflow.run(
                input_file=self.variants_bcf,
                output_dir=self.run_dir,
                db_mode=True,
                nextflow_args=self.nextflow_args
            )
            duration = datetime.now() - start_time
            self.logger.info(f"Annotation completed in {duration.total_seconds():.2f} seconds")

            # Log workflow files
            for name, hash_value in workflow_files.items():
                self.logger.debug(f"Workflow {name} MD5: {hash_value}")

            # Create archive
            archive_path = self._create_archive(self.run_dir)
            self.logger.info(f"Created annotation archive: {archive_path}")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Workflow execution failed: {e.stderr}")
            shutil.rmtree(self.run_dir)
            sys.exit(1)

    def _create_archive(self, run_dir: Path) -> Path:
        """Create a single archive file from annotation run directory"""
        archive_name = f"{run_dir.name}.vepstash"
        archive_path = run_dir.parent / archive_name

        # Create archive with metadata
        with tarfile.open(archive_path, "w:gz") as tar:
            # Add all files from run directory
            tar.add(run_dir, arcname="")

            # Add metadata file
            metadata = {
                "created": datetime.now().isoformat(),
                "format_version": "1.0",
                "run_id": run_dir.name
            }
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as tf:
                json.dump(metadata, tf, indent=2)
                tf_path = Path(tf.name)

            tar.add(tf_path, arcname="metadata.json")
            tf_path.unlink()

        # Remove original directory after successful archive creation
        shutil.rmtree(run_dir)
        return archive_path

    def _read_archive(self, archive_path: Path, extract_path: Optional[Path] = None) -> Dict:
        """Read metadata from annotation archive without extracting"""
        with tarfile.open(archive_path, "r:gz") as tar:
            try:
                meta_file = tar.extractfile("metadata.json")
                if meta_file:
                    metadata = json.loads(meta_file.read().decode())
                    if extract_path:
                        tar.extractall(path=extract_path)
                    return metadata
            except KeyError:
                raise ValueError("Invalid annotation archive: metadata.json not found")
        return {}

    def _create_run_dir(self, db_dir: Path, workflow_hash: str) -> Path:
        """Create a unique directory for this annotation run using timestamp and workflow hash"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        annotations_dir = Path(db_dir) / "annotations"
        run_dir = annotations_dir / f"{timestamp}_{workflow_hash[:8]}"  # Use first 8 chars of hash
        run_dir.mkdir(parents=True, exist_ok=True)
        return run_dir

    def _store_workflow_files(self, workflow_dir: Path, target_dir: Path) -> Dict[str, str]:
        """Store workflow files with their hash for version control"""
        # Get main workflow files
        stored_files = {}

        # Copy main workflow files
        main_files = ['main.nf', 'nextflow.config']
        for name in main_files:
            path = workflow_dir / name
            if not path.exists():
                self.logger.error(f"Required workflow file not found: {path}")
                raise FileNotFoundError(f"Required workflow file not found: {path}")

            # Copy file
            target_path = target_dir / name
            with open(path, 'rb') as src, open(target_path, 'wb') as dst:
                dst.write(src.read())
            stored_files[name] = compute_md5(target_path)

        # Copy module files
        modules_dir = workflow_dir / "modules"
        if modules_dir.exists():
            target_modules_dir = target_dir / "modules"
            target_modules_dir.mkdir(exist_ok=True)

            for module_file in modules_dir.glob("*.nf"):
                rel_path = module_file.relative_to(workflow_dir)
                target_path = target_dir / rel_path
                target_path.parent.mkdir(parents=True, exist_ok=True)

                with open(module_file, 'rb') as src, open(target_path, 'wb') as dst:
                    dst.write(src.read())
                stored_files[str(rel_path)] = compute_md5(target_path)
                self.logger.debug(f"Stored module file: {rel_path}")

        return stored_files


class VCFAnnotator(VEPDatabase):
    """Handles annotation of user VCF files using the database"""

    def __init__(self, db_path: Path, input_vcf: Path, output_dir: Path, workflow_dir: Optional[Path] = None,
                 params_file: Optional[Path] = None, threads: int = 4):
        super().__init__(db_path)
        self.input_vcf = Path(input_vcf)
        self.output_dir = Path(output_dir)
        self.workflow_dir = workflow_dir or self.db_path / "workflow"
        self.params_file = params_file
        self.threads = max(int(threads), 1)
        self._nf_workflow = NextflowWorkflow(self.workflow_dir, params_file=self.params_file)
        self.run_dir = self._create_run_dir(workflow_hash=self._nf_workflow.workflow_hash)

        # Set up logging
        log_file = self.run_dir / "annotation.log"
        self.logger = setup_logger(log_file)
        self.logger.info(f"Initializing VCF annotation for {input_vcf}")
        self.logger.debug(f"Output directory: {output_dir}")
        self.logger.debug(f"Using {threads} threads")

    def _process_region(self, args: tuple) -> pd.DataFrame:
        """Process a single genomic region from BCF file."""
        bcf_path, region = args
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
                    if chrom[:3] != 'chr':
                        continue

                    pos = record.pos
                    ref = record.ref
                    alt = record.alts[0]  # Assuming single ALT
                    if not all([x in BASES for x in alt]):
                        continue

                    # Extract FORMAT fields
                    if len(record.samples) == 0:
                        self.logger.warning(f"No samples found at {record.chrom}:{pos}")
                        continue

                    sample = record.samples[0]
                    ad = sample.get('AD', None)
                    dp = sample.get('DP', None)

                    # Calculate AF
                    af = None
                    if ad and len(ad) >= 2:
                        ref_depth = ad[0]
                        alt_depth = ad[1]
                        af = alt_depth / (ref_depth + alt_depth) if (ref_depth + alt_depth) > 0 else None

                    # Process INFO fields
                    info = {key: None for key in INFO_FIELDS}
                    info |= {key: record.info.get(key, None) for key in INFO_FIELDS if key in record.info}
                    info |= {
                        'GT': sample.get('GT', None),
                        'AD': ad[1] if ad and len(ad) > 1 else None,
                        'DP': dp,
                        'AF': af
                    }

                    # Process clinvar and gnomad fields
                    self._process_variant_annotations(info)

                    # Process VEP annotations
                    vep_csqs = [dict(zip(vcf.header.info['CSQ'].description.split(' ')[-1].split('|'),
                                       x.split("|"))) for x in record.info["CSQ"]]
                    expanded_transcripts = parse_vep_info(vep_csqs)

                    for transcript in expanded_transcripts:
                        row = {
                            "CHROM": chrom,
                            "POS": pos,
                            "REF": ref,
                            "ALT": alt,
                            **info,
                            **transcript
                        }
                        records.append(row)

                except Exception as e:
                    excluded_count += 1
                    self.logger.error(f"Error processing variant at {record.chrom}:{record.pos}: {e}")
                    continue

            self.logger.debug(f"Processed {variant_count} variants, excluded {excluded_count} in {region}")
            return pd.DataFrame(records)

        except Exception as e:
            self.logger.error(f"Error processing region {region}: {e}")
            raise

    def _process_variant_annotations(self, info: dict) -> None:
        """Process clinvar and gnomad annotations."""
        # Process clinvar
        clinvar_clnsig = info.get("clinvar_clnsig", None)
        if clinvar_clnsig and clinvar_clnsig[0] != 'null':
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
        gnomade_ac = int(info.get("gnomade_ac", 0)) if info.get("gnomade_ac", None) else 0
        gnomadg_af = float(info["gnomadg_af"]) if info.get("gnomadg_af", None) else None
        gnomade_af = float(info["gnomade_af"]) if info.get("gnomade_af", None) else None

        info["gnomad_af"] = wavg(gnomadg_af, gnomade_af, gnomadg_ac, gnomade_ac)

        # Remove individual gnomad fields
        for field in gnomad_fields:
            info.pop(field, None)



    def annotate(self, convert_parquet: bool = True) -> Path:
        """Run annotation workflow on input VCF file.

        Args:
            convert_parquet: Whether to convert output to Parquet format

        Returns:
            Path to output file (BCF or Parquet)
        """
        start_time = time.time()
        self.logger.info("Starting VCF annotation")

        try:
            self._validate_inputs()
            annotated_bcf = self._run_workflow()
            vep_time = time.time()
            self.logger.info(f"VEP workflow completed in {vep_time - start_time:.1f}s")

            if not convert_parquet:
                return annotated_bcf

            self.logger.info("Converting to Parquet format")

            output_file = self._convert_to_parquet(annotated_bcf)
            end_time = time.time()

            self.logger.info(f"Parquet conversion completed in {end_time - vep_time:.1f}s")
            self.logger.info(f"Total execution time: {end_time - start_time:.1f}s")


            return output_file

        except Exception as e:
            self.logger.error("Annotation failed", exc_info=True)
            if self.run_dir.exists():
                shutil.rmtree(self.run_dir)
            raise


    def _validate_inputs(self) -> None:
        """Validate input VCF format and requirements"""
        if not self.input_vcf.exists():
            self.logger.error(f"Input VCF file not found: {self.input_vcf}")

            raise FileNotFoundError(f"Input VCF file not found: {self.input_vcf}")

        is_valid, error = validate_vcf_format(self.input_vcf)
        if not is_valid:
            self.logger.error(f"Invalid VCF file: {error}")

            raise ValueError(f"Invalid VCF file: {error}")

        self.ensure_indexed(self.input_vcf)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info("Input validation completed")

    def _run_workflow(self) -> Path:
        """Run the Nextflow annotation workflow"""
        sample_name = self.input_vcf.stem.split('.')[0]
        final_bcf = self.run_dir / f"{sample_name}_norm_final.bcf"

        cmd = [
            "nextflow", "run", str(self.workflow_dir / "main.nf"),
            "--input", str(self.input_vcf),
            "--output", str(self.run_dir),
            "--db_bcf", str(self.variants_bcf),
            "--db_mode", "false",
            "--vep_max_forks", str(min(self.threads, 4)),
            "--vep_max_chr_parallel", str(min(self.threads, 4)),
            "-with-trace"
        ]

        try:
            start_time = datetime.now()
            self.logger.info("Starting Nextflow workflow")
            self.logger.debug(f"Command: {' '.join(cmd)}")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            duration = datetime.now() - start_time

            self._nf_workflow.store_workflow_dag(self.run_dir, cmd)
            stats = get_bcf_stats(final_bcf)

            # Log annotation details
            self.logger.info(f"Workflow completed in {duration.total_seconds():.2f} seconds")
            self.logger.info(f"Output BCF: {final_bcf}")
            self.logger.debug("Output statistics:")
            for key, value in stats.items():
                self.logger.debug(f"{key}: {value}")

            return final_bcf

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Workflow execution failed: {e.stderr}")
            raise RuntimeError(f"Annotation workflow failed: {e}")

    def _convert_to_parquet(self, bcf_path: Path) -> Path:
        """Convert annotated BCF to optimized Parquet format"""
        vcf = pysam.VariantFile(str(bcf_path))
        regions = list(vcf.header.contigs.keys())
        args_list = [(str(bcf_path), region) for region in regions]

        self.logger.info(f"Converting BCF to Parquet: {bcf_path}")
        self.logger.debug(f"Processing {len(regions)} regions using {self.threads} threads")

        with Pool(self.threads) as pool:
            dataframes = pool.map(self._process_region, args_list)

        # Filter and combine dataframes
        dataframes = [df for df in dataframes if not df.empty]
        if not dataframes:
            self.logger.error("No valid variants found in annotated file")
            raise ValueError("No valid variants found in annotated file")

        combined_df = pd.concat(dataframes, ignore_index=True)
        output_file = self.output_dir / f"{self.input_vcf.stem}.parquet"
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
            write_statistics=True
        )
        self.logger.info("Parquet conversion completed")
        return output_file

    def _create_run_dir(self, workflow_hash: str) -> Path:
        """Create a unique directory for this annotation run.

        Args:
            workflow_hash: Hash of the workflow files

        Returns:
            Path to the created run directory
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        annotations_dir = self.output_dir / "annotations"
        run_dir = annotations_dir / f"{timestamp}_{workflow_hash[:8]}"  # Use first 8 chars of hash
        self.logger.debug(f"Creating run directory: {run_dir}")
        run_dir.mkdir(parents=True, exist_ok=True)
        self.logger.debug("Run directory created successfully")

        return run_dir
