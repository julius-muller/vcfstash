#!/bin/bash
set -euo pipefail

# Activate the virtual environment
source /home/appuser/projects/vcfstash/.venv/bin/activate

# Set a variable for the verbose flag
VERBOSE="-v"

# Define allele frequencies ordered from highest to lowest: '1e-1', '5e-2', '1e-2', '1e-3', '0'
AF_VALUES=('1e-1' '5e-2' '1e-2' '1e-3' '0')
# Define the cache versions for the cached annotation runs
CV_VALUES=('ex' 'genex' 'genexd')

for AF in "${AF_VALUES[@]}"; do
    echo "Processing allele frequency: ${AF}"

    # Define the base stash directory for this allele frequency
    STASH_DIR="/mnt/data/samples/bench/vcfcache/gnomad_${AF}"

    # --- Step 1: Stash Initialization ---
    if [ ! -d "${STASH_DIR}" ]; then
      echo "Running stash-init for allele frequency ${AF}"
      vcfstash stash-init $VERBOSE \
        --vcf /mnt/data/resources/gnomad/vcf_gnomad_v4_hg19_exomes/gnomad.exomes.v4.1.sites.grch37.trimmed_liftover_norm_${AF}.bcf \
        --output ${STASH_DIR} \
        -y /home/appuser/projects/vcfstash/tests/config/user_params.yaml
    else
      echo "Skipping stash-init for ${AF} as ${STASH_DIR} exists."
    fi

    # --- Step 2: First Annotation (gnomad_ex) ---
    if [ ! -d "${STASH_DIR}/stash/gnomad_ex_${AF}" ]; then
      echo "Running stash-annotate (gnomad_ex_${AF})"
      vcfstash stash-annotate $VERBOSE \
        --name gnomad_ex_${AF} \
        -a /home/appuser/projects/vcfstash/tests/config/annotation.config \
        --db ${STASH_DIR}
    else
      echo "Skipping stash-annotate (gnomad_ex_${AF}) as ${STASH_DIR}/stash/gnomad_ex_${AF} exists."
    fi

    # --- Step 3: Add Genomes Dataset ---
    if [ ! -f "${STASH_DIR}/genomes_added.flag" ]; then
      echo "Running stash-add (genomes dataset) for ${AF}"
      vcfstash stash-add $VERBOSE \
        --db ${STASH_DIR} \
        -i /mnt/data/resources/gnomad/vcf_gnomad_v4_hg19_genomes/gnomad.genomes.v4.1.sites.grch37.trimmed_liftover_norm_${AF}.bcf
      touch "${STASH_DIR}/genomes_added.flag"
    else
      echo "Skipping stash-add (genomes) for ${AF} as marker file exists."
    fi

    # --- Step 4: Second Annotation (gnomad_genex) ---
    if [ ! -d "${STASH_DIR}/stash/gnomad_genex_${AF}" ]; then
      echo "Running stash-annotate (gnomad_genex_${AF})"
      vcfstash stash-annotate $VERBOSE \
        --name gnomad_genex_${AF} \
        -a /home/appuser/projects/vcfstash/tests/config/annotation.config \
        --db ${STASH_DIR}
    else
      echo "Skipping stash-annotate (gnomad_genex_${AF}) as ${STASH_DIR}/stash/gnomad_genex_${AF} exists."
    fi

    # --- Step 5: Add dbSNP Dataset ---
    if [ ! -f "${STASH_DIR}/dbsnp_added.flag" ]; then
      echo "Running stash-add (dbSNP) for ${AF}"
      vcfstash stash-add $VERBOSE \
        --db ${STASH_DIR} \
        -i /mnt/data/resources/dbsnp/b151_GRCh37p13/common_all_20180423_norms.bcf
      touch "${STASH_DIR}/dbsnp_added.flag"
    else
      echo "Skipping stash-add (dbSNP) for ${AF} as marker file exists."
    fi

    # --- Step 6: Third Annotation (gnomad_genexd) ---
    if [ ! -d "${STASH_DIR}/stash/gnomad_genexd_${AF}" ]; then
      echo "Running stash-annotate (gnomad_genexd_${AF})"
      vcfstash stash-annotate $VERBOSE \
        --name gnomad_genexd_${AF} \
        -a /home/appuser/projects/vcfstash/tests/config/annotation.config \
        --db ${STASH_DIR}
    else
      echo "Skipping stash-annotate (gnomad_genexd_${AF}) as ${STASH_DIR}/stash/gnomad_genexd_${AF} exists."
    fi

    # --- Cached Annotation Runs for Each Cache Version ---
    for CV in "${CV_VALUES[@]}"; do
        echo "Annotating cache version: ${CV} for allele frequency ${AF}"

        # Define the output directory and the annotation directory for this cache version
        OUTPUT_WITH="/mnt/data/samples/bench/vst_results/mgm32_gnomad_${CV}_${AF}"
        ANNOTATION_DIR="${STASH_DIR}/stash/gnomad_${CV}_${AF}"

        if [ ! -d "${OUTPUT_WITH}" ]; then
          echo "Running annotation with cache for ${CV}_${AF}"
          vcfstash annotate $VERBOSE \
            -a ${ANNOTATION_DIR} \
            --vcf /mnt/data/samples/test_mgm/mgm_WGS_32.gatkWGS_norm.bcf \
            --output ${OUTPUT_WITH}
        else
          echo "Skipping annotation with cache for ${CV}_${AF} as ${OUTPUT_WITH} exists."
        fi
    done

done

# --- Single Uncached Annotation Run (Only Once) ---
# Use a fixed output directory (no allele frequency) and a representative annotation directory.
# Here we use the annotation directory from allele frequency '1e-1' and the 'gnomad_ex' step.
UNCACHED_OUTPUT="/mnt/data/samples/bench/vst_results/mgm32_gnomad_uc"
ANNOTATION_DIR_UC="/mnt/data/samples/bench/vst_caches/gnomad_genex_1e-1/stash/gnomad_ex_1e-1"

if [ ! -d "${UNCACHED_OUTPUT}" ]; then
  echo "Running annotation without cache (uncached) once"
  vcfstash annotate $VERBOSE \
    -a ${ANNOTATION_DIR_UC} \
    --vcf /mnt/data/samples/test_mgm/mgm_WGS_32.gatkWGS_norm.bcf \
    --output ${UNCACHED_OUTPUT} \
    --uncached
else
  echo "Skipping uncached annotation as ${UNCACHED_OUTPUT} exists."
fi
