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

      # Basic header check (quick check only)
      ${params.bcftools_cmd} view -h "${vcf}" | head -n 1 >/dev/null || { echo "Invalid VCF/BCF format"; exit 1; }

      # ===========================================
      # Reference genome compatibility check - OPTIMIZED FOR LARGE FILES
      # ===========================================
      echo "Validating reference genome compatibility with VCF..."
      echo "Loading reference chromosome names..."
      cut -f1 ${reference_index} > ref_chroms.txt

      # Get a sample of chromosomes from the VCF - don't process the entire file
      echo "Sampling chromosomes from VCF (using index for efficiency)..."
      ${params.bcftools_cmd} index -s ${vcf} 2>/dev/null | cut -f1 > vcf_chroms.txt
      
      # If the above fails (some BCF versions don't support -s), try a limited sample
      if [ ! -s vcf_chroms.txt ]; then
        echo "Fallback: extracting sample chromosomes from VCF directly..."
        ${params.bcftools_cmd} view -H ${vcf} | head -n 1000 | cut -f1 | sort | uniq > vcf_chroms.txt
      fi
      
      # Show first few variants for debugging
      echo "First few variants in VCF:"
      ${params.bcftools_cmd} view -H ${vcf} | head -n 3
      
      # Create chromosome mapping lookup tables for efficiency
      echo "Loading chromosome mapping from ${chr_add}..."
      awk '{print \$1"\\t"\$2}' ${chr_add} > chr_mapping.txt
      
      # Check if we need to adjust chromosome names based on the first chromosome
      FIRST_CHR=\$(head -n 1 vcf_chroms.txt)
      echo "First chromosome in VCF: \$FIRST_CHR"
      
      # Determine if we need prefix adjustments
      NEEDS_CHR_PREFIX=0
      NEEDS_CHR_REMOVAL=0
      
      if [[ "\$FIRST_CHR" == chr* ]] && ! grep -q "^chr" ref_chroms.txt; then
          echo "VCF uses 'chr' prefix but reference doesn't. Will adjust for comparison."
          NEEDS_CHR_REMOVAL=1
      elif [[ "\$FIRST_CHR" != chr* ]] && grep -q "^chr" ref_chroms.txt; then
          echo "Reference uses 'chr' prefix but VCF doesn't. Will adjust for comparison."
          NEEDS_CHR_PREFIX=1  
      fi
      
      # Process chromosomes in bulk using AWK for efficiency
      echo "Validating chromosomes against reference..."
      awk -v needs_prefix=\$NEEDS_CHR_PREFIX -v needs_removal=\$NEEDS_CHR_REMOVAL '
      # First load the reference chromosomes into a lookup hash
      FILENAME == "ref_chroms.txt" {
          ref_chroms[\$1] = 1;
          next;
      }
      
      # Then load chromosome mappings into a lookup hash
      FILENAME == "chr_mapping.txt" {
          chr_map[\$1] = \$2;
          rev_map[\$2] = \$1;
          next;
      }
      
      # Finally, check each VCF chromosome against reference
      FILENAME == "vcf_chroms.txt" {
          chrom = \$1;
          if (chrom == "") next;
          
          # Try different mappings in this order:
          # 1. Direct mapping from chr_add.txt
          if (chrom in chr_map) {
              ref_chrom = chr_map[chrom];
          }
          # 2. Standard prefix adjustments
          else if (needs_prefix == 1) {
              ref_chrom = "chr" chrom;
          }
          else if (needs_removal == 1 && chrom ~ /^chr/) {
              ref_chrom = substr(chrom, 4);  # Remove chr prefix
          }
          else {
              ref_chrom = chrom;  # No adjustment needed
          }
          
          # Check if this chromosome exists in reference
          if (ref_chrom in ref_chroms) {
              # Found direct match
              next;
          }
          else if (ref_chrom in rev_map && rev_map[ref_chrom] in ref_chroms) {
              # Found via reverse mapping
              next;
          }
          else {
              # Not found - report error
              print "ERROR: Chromosome " chrom " (adjusted to " ref_chrom ") not found in reference genome!";
              error_count++;
          }
      }
      
      END {
          exit error_count > 0 ? 1 : 0;
      }
      ' ref_chroms.txt chr_mapping.txt vcf_chroms.txt
      
      # Check if any chromosomes are missing
      if [ \$? -ne 0 ]; then
          echo "ERROR: One or more chromosomes from the VCF are missing in the reference genome."
          echo "This suggests the reference genome and VCF use different chromosome naming conventions,"
          echo "or the reference is incomplete for this dataset."
          echo "Available chromosomes in reference:"
          head -n 10 ref_chroms.txt
          echo "First 10 chromosomes in VCF:"
          head -n 10 vcf_chroms.txt
          echo "Consider updating chr_add.txt to include additional mappings."
          exit 1
      fi

      echo "All sampled chromosomes in VCF are present in the reference genome or have valid mappings."
      
      # ===========================================
      # Reference allele compatibility check - LIMITED SAMPLING FOR LARGE FILES
      # ===========================================
      echo "Validating reference allele compatibility with small sample of variants..."
      
      # Only check a very small number of variants for large files
      SAMPLE_SIZE=3
      
      # Extract a random sampling of variants for testing
      # Use head/tail combination to get variants from different parts of the file
      ${params.bcftools_cmd} view -H ${vcf} | head -n 100 | tail -n \$SAMPLE_SIZE > sample_variants.txt
      
      # Initialize counter for mismatches
      REF_MISMATCHES=0
      TOTAL_CHECKED=0
      
      while read -r variant_line; do
          # Parse VCF fields
          CHROM=\$(echo "\$variant_line" | cut -f1)
          POS=\$(echo "\$variant_line" | cut -f2)
          REF_ALLELE=\$(echo "\$variant_line" | cut -f4)
          
          # Skip if the REF allele is not a simple nucleotide sequence
          if [[ ! \$REF_ALLELE =~ ^[ACGTN]+\$ ]]; then
              echo "Skipping variant with non-standard REF allele: \$REF_ALLELE"
              continue
          fi
          
          # Get the correct reference chromosome name using the same mapping logic
          QUERY_CHROM="\$CHROM"
          
          # Apply chromosome mapping
          MAPPED_CHROM=\$(grep -m 1 "^\$CHROM\t" chr_mapping.txt | cut -f2 || echo "")
          if [ ! -z "\$MAPPED_CHROM" ]; then
              QUERY_CHROM="\$MAPPED_CHROM"
          elif [ \$NEEDS_CHR_PREFIX -eq 1 ]; then
              QUERY_CHROM="chr\$CHROM"
          elif [ \$NEEDS_CHR_REMOVAL -eq 1 ]; then
              QUERY_CHROM=\${CHROM#"chr"}
          fi
          
          # Use samtools to extract the reference sequence at this position
          REF_BASE=\$(samtools faidx ${reference} \$QUERY_CHROM:\$POS-\$((\$POS + \${#REF_ALLELE} - 1)) | grep -v ">" | tr -d '\\n')
          
          echo "Checking variant \$CHROM:\$POS (querying as \$QUERY_CHROM)"
          echo "  VCF REF allele: \$REF_ALLELE"
          echo "  Genome sequence: \$REF_BASE"
          
          TOTAL_CHECKED=\$((TOTAL_CHECKED + 1))
          
          # Compare the sequences
          if [[ "\$REF_ALLELE" != "\$REF_BASE" ]]; then
              echo "WARNING: Reference mismatch at \$CHROM:\$POS (querying as \$QUERY_CHROM:\$POS)"
              echo "  VCF REF allele: \$REF_ALLELE"
              echo "  Reference genome: \$REF_BASE"
              REF_MISMATCHES=\$((REF_MISMATCHES + 1))
          fi
      done < sample_variants.txt
      
      # Report the results of allele checking
      if [ \$REF_MISMATCHES -gt 0 ]; then
          MISMATCH_PERCENT=\$((REF_MISMATCHES * 100 / TOTAL_CHECKED))
          echo "WARNING: Found \$REF_MISMATCHES reference allele mismatches out of \$TOTAL_CHECKED variants checked (\$MISMATCH_PERCENT%)."
          echo "This suggests the VCF may have been created with a different reference genome version."
          
          if [ \$REF_MISMATCHES -gt 0 ]; then
              echo "ERROR: Reference allele mismatches detected (\$REF_MISMATCHES out of \$TOTAL_CHECKED)."
              echo "The VCF file appears to use a different reference genome than the one provided."
              echo "Normalization is likely to fail. Please check your reference genome."
              exit 1
          fi
      else
          echo "Reference allele check passed: All sampled variants match the reference genome."
      fi
      
      echo "Reference validation completed successfully."
      """
  }
e

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