process CaptureToolVersions {
    publishDir "${output_dir}", mode: 'copy'

    input:
    val sample_name
    path output_dir

    output:
    path "${sample_name}_final.tool_version.log", emit: version_log

    script:
    """
	echo "${params.bcftools_cmd} version:" > ${sample_name}_final.tool_version.log
    ${params.bcftools_cmd} --version >> ${sample_name}_final.tool_version.log
    """

    stub:
    """
    touch ${sample_name}_final.tool_version.log
    """
}

process ValidateInputs {
      stageInMode 'symlink'

      input:
      path vcf
      path vcf_index  // Explicitly passing the index file
      path reference
      path reference_index
      path chr_add

      output:
      val true, emit: validated

      script:
      """
      # Check file existence
      test -e "${vcf}" || { echo "VCF file not found: ${vcf}"; exit 1; }
      test -e "${vcf_index}" || { echo "VCF index file not found: ${vcf_index}"; exit 1; }
      test -e "${reference}" || { echo "Reference file not found: ${reference}"; exit 1; }
      test -e "${reference_index}" || { echo "Reference index file not found: ${reference_index}"; exit 1; }
      test -e "${chr_add}" || { echo "Chr add file not found: ${chr_add}"; exit 1; }

      # Basic header check
      ${params.bcftools_cmd} view -h "${vcf}" | head -n 1 >/dev/null || { echo "Invalid VCF/BCF format"; exit 1; }

      # ===========================================
      # Reference genome compatibility check
      # ===========================================
      echo "Validating reference genome compatibility with VCF..."

      # Get the list of available chromosomes in the reference
      echo "Available chromosomes in reference:"
      cut -f1 ${reference_index}

      # Get all unique chromosomes from the VCF file
      echo "Extracting chromosomes from VCF..."
      ${params.bcftools_cmd} view -H ${vcf} | cut -f1 | sort | uniq > vcf_chroms.txt
      
      # Get all chromosomes from the reference index
      cut -f1 ${reference_index} > ref_chroms.txt
      
      # Show first few variants for debugging
      echo "First few variants in VCF:"
      ${params.bcftools_cmd} view -H ${vcf} | head -n 3
      
      # Check if we need to adjust chromosome names (chr prefix handling)
      NEEDS_CHR_PREFIX=0
      NEEDS_CHR_REMOVAL=0
      
      # Sample the first chromosome from VCF
      FIRST_CHR=\$(head -n 1 vcf_chroms.txt)
      echo "First chromosome in VCF: \$FIRST_CHR"
      
      # Check if VCF has chr prefix but reference doesn't
      if [[ "\$FIRST_CHR" == chr* ]] && ! grep -q "^chr" ref_chroms.txt; then
          echo "VCF uses 'chr' prefix but reference doesn't. Will adjust for comparison."
          NEEDS_CHR_REMOVAL=1
      fi
      
      # Check if reference has chr prefix but VCF doesn't
      if [[ "\$FIRST_CHR" != chr* ]] && grep -q "^chr" ref_chroms.txt; then
          echo "Reference uses 'chr' prefix but VCF doesn't. Will adjust for comparison."
          NEEDS_CHR_PREFIX=1  
      fi
      
      # Validate all chromosomes in the VCF exist in the reference
      echo "Validating all chromosomes..."
      MISSING_CHROMS=0
      
      while read -r chrom; do
          # Skip empty lines
          [ -z "\$chrom" ] && continue
          
          REF_CHROM="\$chrom"
          
          # Adjust chromosome name if needed
          if [ \$NEEDS_CHR_REMOVAL -eq 1 ]; then
              REF_CHROM=\${chrom#"chr"}
          elif [ \$NEEDS_CHR_PREFIX -eq 1 ]; then
              REF_CHROM="chr\$chrom"
          fi
          
          # Check if this chromosome exists in reference
          if ! grep -q "^\$REF_CHROM\$" ref_chroms.txt && ! grep -q "^\$REF_CHROM\t" ref_chroms.txt; then
              echo "ERROR: Chromosome '\$chrom' (adjusted to '\$REF_CHROM' for reference) not found in reference genome!"
              MISSING_CHROMS=1
          fi
      done < vcf_chroms.txt
      
      # Fail if any chromosomes are missing
      if [ \$MISSING_CHROMS -eq 1 ]; then
          echo "ERROR: One or more chromosomes from the VCF are missing in the reference genome."
          echo "This suggests the reference genome and VCF use different chromosome naming conventions,"
          echo "or the reference is incomplete for this dataset."
          echo "Available chromosomes in reference:"
          cat ref_chroms.txt
          echo "Chromosomes in VCF:"
          cat vcf_chroms.txt
          exit 1
      fi

      echo "All chromosomes in VCF are present in the reference genome."
      
      # ===========================================
      # Reference allele compatibility check
      # ===========================================
      echo "Validating reference allele compatibility with VCF..."
      echo "This checks that the REF allele in the VCF matches the reference genome"
      
      # Extract a few variants for testing (first 5 variants)
      ${params.bcftools_cmd} view -H ${vcf} | head -n 5 > sample_variants.txt
      
      # Initialize counter for mismatches
      REF_MISMATCHES=0
      TOTAL_CHECKED=0
      
      while read -r variant_line; do
          # Parse VCF fields
          CHROM=$(echo "$variant_line" | cut -f1)
          POS=$(echo "$variant_line" | cut -f2)
          REF_ALLELE=$(echo "$variant_line" | cut -f4)
          
          # Skip if the REF allele is not a simple nucleotide sequence
          if [[ ! $REF_ALLELE =~ ^[ACGTN]+$ ]]; then
              echo "Skipping variant with non-standard REF allele: $REF_ALLELE"
              continue
          fi
          
          # Adjust chromosome name if needed
          QUERY_CHROM="$CHROM"
          if [ \$NEEDS_CHR_REMOVAL -eq 1 ]; then
              QUERY_CHROM=\${CHROM#"chr"}
          elif [ \$NEEDS_CHR_PREFIX -eq 1 ]; then
              QUERY_CHROM="chr\$CHROM"
          fi
          
          # Use samtools to extract the reference sequence at this position
          REF_BASE=$(samtools faidx ${reference} $QUERY_CHROM:$POS-$(($POS + ${#REF_ALLELE} - 1)) | grep -v ">" | tr -d '\\n')
          
          echo "Checking variant $CHROM:$POS (querying as $QUERY_CHROM)"
          echo "  VCF REF allele: $REF_ALLELE"
          echo "  Genome sequence: $REF_BASE"
          
          # Compare the sequences
          if [[ "$REF_ALLELE" != "$REF_BASE" ]]; then
              echo "WARNING: Reference mismatch at $CHROM:$POS (querying as $QUERY_CHROM:$POS)"
              echo "  VCF REF allele: $REF_ALLELE"
              echo "  Reference genome: $REF_BASE"
              REF_MISMATCHES=$((REF_MISMATCHES + 1))
          fi
          
          TOTAL_CHECKED=$((TOTAL_CHECKED + 1))
      done < sample_variants.txt
      
      # Report the results
      if [ $REF_MISMATCHES -gt 0 ]; then
          MISMATCH_PERCENT=$((REF_MISMATCHES * 100 / TOTAL_CHECKED))
          echo "WARNING: Found $REF_MISMATCHES reference allele mismatches out of $TOTAL_CHECKED variants checked ($MISMATCH_PERCENT%)."
          echo "This suggests the VCF may have been created with a different reference genome version."
          
          if [ $REF_MISMATCHES -gt 2 ]; then
              echo "ERROR: Too many reference allele mismatches detected ($REF_MISMATCHES out of $TOTAL_CHECKED)."
              echo "The VCF file appears to use a different reference genome than the one provided."
              echo "Normalization is likely to fail. Please check your reference genome."
              exit 1
          else
              echo "Continuing with caution - a small number of mismatches may be acceptable."
          fi
      else
          echo "Reference allele check passed: All sampled variants match the reference genome."
      fi
      
      echo "Reference validation passed."
      """
  }


workflow UTILS {
    take:
    sample_name
    output_dir
    vcf
    chr_add
    reference

    main:
    // Make sure reference file index exists
    reference_index = file("${reference}.fai")

    // Explicitly look for VCF index file - try both csi and tbi formats
    vcf_index = file("${vcf}.{csi,tbi}", checkIfExists: true)[0]

    validateInputResult = ValidateInputs(
        vcf,
        vcf_index,
        reference,
        reference_index,
        chr_add,
    )

    emit:
    validate = validateInputResult.validated
}