
import shutil
import sys
import time
from vcfstash.database.base import VCFDatabase, NextflowWorkflow
from vcfstash.database.outputs import AnnotatedStashOutput, AnnotatedUserOutput
from pathlib import Path
import subprocess
from datetime import datetime
from multiprocessing import Pool
import pysam
from vcfstash.utils.logging import setup_logging

class DatabaseAnnotator(VCFDatabase):
    """Handles database annotation workflows"""

    def __init__(self, annotation_name: str, db_path: Path | str, anno_config_file: Path | str,
                 params_file: Path | str = None, config_file: Path | str = None, verbosity: int = 0,
                 force: bool = False, debug: bool = False):
        """Initialize database annotator.

        self = DatabaseAnnotator(annotation_name="testor", anno_config_file=Path('~/projects/vcfstash/tests/config/annotation.config'),
         db_path=Path('~/tmp/vcfstash/test_stash'),force=True)

        Args:
            db_path: Path to the database
            workflow_dir: Optional workflow directory (defaults to <db_path>/workflow)
            params_file: Optional parameters file
            verbosity: Logging verbosity level (0=WARNING, 1=INFO, 2=DEBUG)
        """

        super().__init__(db_path, verbosity, debug)

        self.stashed_annotations = AnnotatedStashOutput(str(self.stash_dir / annotation_name))
        self.stashed_annotations.validate_label(annotation_name)
        self.annotation_name = annotation_name
        self.logger = self.connect_loggers()
        self._setup_annotation_stash(force)
        self.output_dir = self.stashed_annotations.annotation_dir

        self.config_file = None
        if config_file:
            self.config_file = self.output_dir / 'process.config'
            shutil.copyfile(config_file.expanduser().resolve(), self.config_file)

        self.info_snapshot_file = self.output_dir / "blueprint_snapshot.info"

        anno_config_file = Path(anno_config_file).expanduser().resolve()
        self.anno_config_file = self.output_dir / 'annotation.config'
        if anno_config_file != self.anno_config_file:
            shutil.copyfile(anno_config_file.expanduser().resolve() , self.anno_config_file)


        self.params_file = self.output_dir / f'annotation.yaml'
        if params_file:
            shutil.copyfile(params_file.expanduser().resolve(), self.params_file)
        else:
            wfi = self.workflow_dir / "init.yaml"
            assert wfi.exists(), f"Workflow init params file not found: {wfi}"
            shutil.copyfile(wfi, self.params_file)
            assert self.params_file.exists(), f"Workflow params file not found: {self.params_file}"


        # Initialize NextflowWorkflow
        self.nx_workflow = NextflowWorkflow(
            input_file=self.blueprint_bcf,
            output_dir=self.output_dir,
            name=self.annotation_name,
            workflow=self.workflow_dir / "main.nf",
            config_file=self.config_file,
            anno_config_file=self.anno_config_file,
            params_file=self.params_file,
            verbosity=self.verbosity
        )

        # Log initialization parameters
        self.logger.info("Initializing database annotation")
        self.logger.debug(f"Annotation directory: {self.output_dir}")
        self.logger.debug(f"Config file: {self.config_file}")

    def _setup_annotation_stash(self, force:bool) -> None:
        # Remove destination directory if it exists to ensure clean copy
        if self.stashed_annotations.annotation_dir.exists():
            if self.stashed_output.validate_structure():  # we dont want to remove a random dir....
                if force:
                    print(f"Stash directory already exists, removing: {self.stashed_output.root_dir}")
                    shutil.rmtree(self.stashed_annotations.annotation_dir)
                else:
                    raise FileExistsError(
                        f"Output directory already exists: {self.stashed_annotations.annotation_dir}\nIf intended, use --force to overwrite.")
            else:
                if not force:
                    raise FileNotFoundError(
                        f"Output directory must not exist if --force is not set and a valid stash directory: {self.stashed_annotations.annotation_dir}")

        print(f"Creating stash structure: {self.stashed_annotations.annotation_dir}")
        self.stashed_annotations.create_structure()

    def annotate(self, extra_files:bool = True) -> None:
        """Run annotation workflow on database"""

        # Store blueprint snapshot and workflow files
        shutil.copy2(self.info_file, self.info_snapshot_file)

        try:
            self.logger.info("Starting annotation workflow")

            start_time = datetime.now()
            self.nx_workflow.run(
                db_mode='stash-annotate',
                db_bcf= self.blueprint_bcf,
                trace=extra_files,
                dag=extra_files,
                report=extra_files
            )
            if not self.debug:
                self.nx_workflow.cleanup_work_dir()

            duration = datetime.now() - start_time
            self.logger.info(f"Annotation to {self.output_dir} completed in {duration.total_seconds():.2f} seconds")


        except subprocess.CalledProcessError as e:
            self.logger.error(f"Workflow execution failed: {e.stderr}")
            self.logger.warning(f"Removing output directory: {self.output_dir}")
            shutil.rmtree(self.output_dir)
            sys.exit(1)



class VCFAnnotator(VCFDatabase):
    """Handles annotation of user VCF files using the database"""

    INFO_FIELDS = ['GT', 'DP', 'AF', "gnomadg_af", "gnomade_af", "gnomadg_ac", "gnomade_ac", 'clinvar_clnsig',
                   "deeprvat_score"]
    BASES = {"A", "C", "G", "T"}

    def __init__(self, input_vcf: Path | str, annotation_db: Path | str, output_dir: Path | str,
                 config_file: Path | str = None, params_file: Path | str = None, verbosity: int = 0,
                 force: bool = False, debug: bool = False):
        """Initialize database annotator.

        Args:
            input_vcf: Path to the input BCF/VCF file, needs to be indexed!
            annotation_db: Path to the annotation database
            output_dir: Path to the output directory
            config_file: Optional configuration file with process configurations for resoureces.
            params_file: Optional parameters file defaults to <db_path>/workflow/init.yaml
            force: Whether to overwrite existing output directory
            debug: Whether to enable debug mode
            verbosity: Logging verbosity level (0=WARNING, 1=INFO, 2=DEBUG)
        """

        self.input_vcf = Path(input_vcf).expanduser().resolve()
        self.vcf_name, fext = self._validate_and_extract_sample_name()

        if not self.input_vcf.exists():
            raise FileNotFoundError(f"Input VCF file not found: {self.input_vcf}")

        self.stashed_annotations = AnnotatedStashOutput(str(annotation_db))
        if not self.stashed_annotations.validate_structure():
            raise FileNotFoundError(f"Annotation database annotation_db not valid: {self.stashed_annotations.annotation_dir}")
        self.annotation_db_path = self.stashed_annotations.annotation_dir
        self.annotation_name = self.stashed_annotations.name
        super().__init__(self.stashed_annotations.stash_output.root_dir, verbosity, debug)


        self.output_annotations = AnnotatedUserOutput(str(output_dir))
        self.output_annotations.validate_label(self.output_annotations.name)

        self._setup_output(force=force)
        self.output_dir = self.output_annotations.root_dir
        self.logger = setup_logging(
            verbosity=self.verbosity,
            log_file=self.output_dir / "annotation.log",
        )

        self.annotation_wfl_path = self.output_annotations.workflow_dir
        self.annotation_wfl_path.mkdir(parents=True, exist_ok=True)

        self.output_vcf = Path(self.output_dir / f'{self.vcf_name}_vst{fext}')

        # now also import the mandatory annotation file, that cannot be provided by the user at this stage
        self.anno_config_file = self.annotation_db_path / 'annotation.config'
        if not self.anno_config_file.exists():
            raise FileNotFoundError(f"Annotation config file not found: {self.anno_config_file}")

        self.config_file = Path(config_file).expanduser().resolve() if config_file else self.annotation_db_path / 'process.config'
        if not self.config_file or not self.config_file.exists():
            self.config_file = None
        else:
            shutil.copyfile(self.config_file, self.annotation_wfl_path / self.config_file.name)
            self.config_file = self.annotation_wfl_path / self.config_file

        self.params_file = self.annotation_wfl_path / f'annotation.yaml'
        if params_file:
            shutil.copyfile(params_file.expanduser().resolve(), self.params_file)
        else:
            wfi = self.annotation_db_path / "annotation.yaml"
            assert wfi.exists(), f"Workflow annotation params file not found: {wfi}, required if no yaml provided!"
            shutil.copyfile(wfi, self.params_file)
        assert self.params_file.exists(), f"Workflow params file not found: {self.params_file}"


        self.stash_file = self.annotation_db_path / "vcfstash_annotated.bcf"
        if not self.stash_file.exists():
            raise FileNotFoundError(f"Stash file not found: {self.stash_file}")

        # Initialize NextflowWorkflow
        self.nx_workflow = NextflowWorkflow(
            input_file=self.input_vcf,
            output_dir=self.output_dir,
            name=self.annotation_name,
            workflow=self.workflow_dir / "main.nf",
            config_file=self.config_file,
            anno_config_file=self.anno_config_file,
            params_file=self.params_file,
            verbosity=self.verbosity
        )

        # Log initialization parameters
        self.logger.info(f"Initializing annotation of {self.input_vcf.name}")
        self.logger.debug(f"Cache file: {self.stash_file}")
        self.logger.debug(f"Config file: {self.config_file}")


    def _setup_output(self, force:bool) -> None:
        # Remove destination directory if it exists to ensure clean copy
        if self.output_annotations.root_dir.exists():
            if self.output_annotations.validate_structure():  # we dont want to remove a random dir....
                if force:
                    print(f"Output directory already exists, removing: {self.output_annotations.root_dir}")
                    shutil.rmtree(self.output_annotations.root_dir)
                else:
                    raise FileExistsError(
                        f"Output directory already exists: {self.output_annotations.root_dir}\nIf intended, use --force to overwrite.")
            else:
                raise FileNotFoundError(
                    f"Output directory must not exist if --force is not set and a valid output directory: {self.output_annotations.root_dir}")

        print(f"Creating output structure: {self.output_annotations.root_dir}")
        self.output_annotations.create_structure()


    def _validate_and_extract_sample_name(self) -> tuple[str, str]:
        """
        Validates the input VCF file has an acceptable extension
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
            raise ValueError(f"Input VCF file '{input_vcf_path}' lacks a file extension.")

        # Check for valid extensions, considering multi-part extensions
        if input_vcf_path.name.endswith('.vcf.gz'):
            extension = '.vcf.gz'
            sample_name = input_vcf_path.name[:-7]  # Removes '.vcf.gz'
            index_file = input_vcf_path.with_suffix(input_vcf_path.suffix + '.tbi')
        elif input_vcf_path.name.endswith('.bcf'):
            extension = '.bcf'
            sample_name = input_vcf_path.name[:-4]  # Removes '.bcf'
            index_file = input_vcf_path.with_suffix('.bcf.csi')
        elif input_vcf_path.name.endswith('.vcf'):
            extension = '.vcf'
            sample_name = input_vcf_path.name[:-4]  # Removes '.vcf'
            index_file = input_vcf_path.with_suffix('.vcf.tbi')
        else:
            raise ValueError(
                f"Input VCF file '{input_vcf_path}' must end with one of {self.VALID_VCF_EXTENSIONS}"
            )

        if not index_file.exists():
            raise FileNotFoundError(f"Index file for '{input_vcf_path}' not found: '{index_file}'")

        return sample_name, extension


    def _process_region(self, args: tuple) -> "pd.DataFrame":
        """Process a single genomic region from BCF file."""
        try:
            import pandas as pd

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
                        if not all([x in self.BASES for x in alt]):
                            continue

                        # Process INFO fields
                        info = {key: None for key in self.INFO_FIELDS}
                        info |= {key: record.info.get(key, None) for key in self.INFO_FIELDS if key in record.info}

                        # Extract FORMAT fields
                        if len(record.samples):

                            sample = record.samples[0]
                            ad = sample.get('AD', None)
                            dp = sample.get('DP', None)

                            # Calculate AF
                            af = None
                            if ad and len(ad) >= 2:
                                ref_depth = ad[0]
                                alt_depth = ad[1]
                                af = alt_depth / (ref_depth + alt_depth) if (ref_depth + alt_depth) > 0 else None

                            info |= {
                                'GT': sample.get('GT', None),
                                'AD': ad[1] if ad and len(ad) > 1 else None,
                                'DP': dp,
                                'AF': af
                            }

                        # Process clinvar and gnomad fields
                        self._process_variant_annotations(info)

                        # Process VCF annotations
                        # TODO: currently does not work for any Tag apart from CSQ, need to pull from .config via ${params.must_contain_info_tag}
                        vcf_csqs = [dict(zip(vcf.header.info['CSQ'].description.split(' ')[-1].split('|'),
                                           x.split("|"))) for x in record.info["CSQ"]]
                        # TODO: This currently only works with VEP annotations, need to be more flexible
                        expanded_transcripts = self.parse_vcf_info(vcf_csqs)

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
        except ImportError:
            raise ImportError(
                "This feature requires pandas and pyarrow. "
                "Install with 'pip install pandas pyarrow'"
            )

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

        info["gnomad_af"] = self.wavg(gnomadg_af, gnomade_af, gnomadg_ac, gnomade_ac)

        # Remove individual gnomad fields
        for field in gnomad_fields:
            info.pop(field, None)



    def annotate(self, uncached: bool = False, convert_parquet: bool = False) -> None:
        """Run annotation workflow on input VCF file.

        Args:
            convert_parquet: Whether to convert output to Parquet format
            unchached: Whether to run the workflow in uncached mode

        Returns:
            Path to output file (BCF or Parquet)
            self = VCFAnnotator(input_vcf="~/projects/vcfstash/tests/data/nodata/sample4.bcf",
             annotation_db = "~/tmp/test/test_out/stash/testor", output_dir="~/tmp/test/aout" ,force=True)

        """
        start_time = time.time()
        self.logger.info("Starting VCF annotation")

        try:


            # Run the workflow in database mode
            self.nx_workflow.run(
                db_mode='annotate' if not uncached else 'annotate-nocache',
                db_bcf=self.stash_file,
                trace=True,
                dag=True,
                report=True
            )
            duration = time.time() - start_time
            self.logger.info(f"VCF workflow completed in {duration:.1f}s")


        except Exception as e:
            self.logger.error("Annotation failed", exc_info=True)
            raise

        if convert_parquet:
            # This should be implemented in nextflow in the long run to leverage resource management
            #threads = self.nx_workflow.nf_config_content['params'].get('vep_max_forks',1) * self.nx_workflow.nf_config_content['params'].get('vep_max_chr_parallel', 1)
            self._convert_to_parquet(self.output_vcf)# , threads=threads)

        if not self.debug:
            self.nx_workflow.cleanup_work_dir()


    def _convert_to_parquet(self, bcf_path: Path, threads:int = 2) -> Path:
        """Convert annotated BCF to optimized Parquet format"""
        try:
            import pandas as pd
            import pyarrow as pa
            import pyarrow.parquet as pq

            vcf = pysam.VariantFile(str(bcf_path))
            regions = list(vcf.header.contigs.keys())
            args_list = [(str(bcf_path), region) for region in regions]

            self.logger.info(f"Converting BCF to Parquet: {bcf_path}")
            self.logger.debug(f"Processing {len(regions)} regions using {threads} threads")

            with Pool(threads) as pool:
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

        except ImportError:
            raise ImportError(
                "Converting to Parquet requires additional dependencies. "
                "Please install them with: pip install pandas pyarrow"
            )

    @staticmethod
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
