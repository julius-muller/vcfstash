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
      cut -f1 ${reference_index} > ref_chroms.txt

      # Get all unique chromosomes from the VCF file
      echo "Extracting chromosomes from VCF..."
      ${params.bcftools_cmd} view -H ${vcf} | cut -f1 | sort | uniq > vcf_chroms.txt
      
      # Show first few variants for debugging
      echo "First few variants in VCF:"
      ${params.bcftools_cmd} view -H ${vcf} | head -n 3
      
      # Load the chromosome mapping from chr_add.txt file
      echo "Loading chromosome mapping from ${chr_add}..."
      
      # Create more efficient lookup tables (using associative arrays)
      # Direct lookup for chromosome mapping and reverse mapping
      awk '{print "CHR_MAP["$1"]="$2}' ${chr_add} > chr_mapping.awk
      awk '{print "CHR_REV_MAP["$2"]="$1}' ${chr_add} >> chr_mapping.awk
      
      # Sample the first chromosome from VCF
      FIRST_CHR=\$(head -n 1 vcf_chroms.txt)
      echo "First chromosome in VCF: \$FIRST_CHR"
      
      # Check if we need to adjust chromosome names (chr prefix handling)
      NEEDS_CHR_PREFIX=0
      NEEDS_CHR_REMOVAL=0
      
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
      
      # Add prefix handling logic to the AWK script
      cat >> chr_mapping.awk << EOF
      BEGIN {
          needs_chr_prefix = ${NEEDS_CHR_PREFIX};
          needs_chr_removal = ${NEEDS_CHR_REMOVAL};
      }
      
      function getRefChrom(chrom) {
          # First try direct mapping
          if (chrom in CHR_MAP) {
              return CHR_MAP[chrom];
          }
          # Then try prefix logic
          else if (needs_chr_prefix) {
              return "chr" chrom;
          }
          else if (needs_chr_removal && chrom ~ /^chr/) {
              return substr(chrom, 4);  # Remove first 3 chars ("chr")
          }
          else {
              return chrom;  # No change needed
          }
      }
      
      function checkRefChrom(chrom, ref_chrom) {
          # Check if ref_chrom exists in reference
          if (ref_chrom in REF_EXISTS) {
              return 1;  # Found
          }
          
          # Try reverse mapping for special cases (MT->chrM->M)
          if (ref_chrom in CHR_REV_MAP) {
              rev_chrom = CHR_REV_MAP[ref_chrom];
              if (rev_chrom in REF_EXISTS) {
                  return 1;  # Found via reverse mapping
              }
          }
          
          return 0;  # Not found
      }
EOF
      
      # Load reference chromosomes into AWK array
      while read -r ref; do
          echo "REF_EXISTS[\\"$ref\\"] = 1;" >> chr_mapping.awk
      done < ref_chroms.txt
      
      # Validate all chromosomes in the VCF exist in the reference
      echo "Validating all chromosomes..."
      
      # Use AWK for efficient validation of all chromosomes at once
      awk -f chr_mapping.awk '
      {
          chrom = $0;
          if (chrom == "") next;
          
          ref_chrom = getRefChrom(chrom);
          
          # If ref_chrom is found in reference or via mapping
          found = checkRefChrom(chrom, ref_chrom);
          
          if (!found) {
              print "ERROR: Chromosome " chrom " (adjusted to " ref_chrom ") not found in reference genome!";
              missing_count++;
          }
      }
      END {
          exit missing_count > 0 ? 1 : 0;
      }' vcf_chroms.txt
      
      # Check exit status
      MISSING_CHROMS=$?
      
      # Fail if any chromosomes are missing
      if [ \$MISSING_CHROMS -ne 0 ]; then
          echo "ERROR: One or more chromosomes from the VCF are missing in the reference genome."
          echo "This suggests the reference genome and VCF use different chromosome naming conventions,"
          echo "or the reference is incomplete for this dataset."
          echo "Available chromosomes in reference:"
          cat ref_chroms.txt
          echo "Chromosomes in VCF:"
          cat vcf_chroms.txt
          echo "Consider updating chr_add.txt to include additional mappings."
          exit 1
      fi

      echo "All chromosomes in VCF are present in the reference genome or have valid mappings."
      
      # ===========================================
      # Reference allele compatibility check - OPTIMIZED for large VCFs
      # ===========================================
      echo "Validating reference allele compatibility with VCF..."
      
      # Only sample up to 5 variants to validate reference compatibility
      # This is a reasonable tradeoff for large VCFs
      SAMPLE_SIZE=5
      
      # Extract a limited number of variants for faster testing
      ${params.bcftools_cmd} view -H ${vcf} | head -n ${SAMPLE_SIZE} > sample_variants.txt
      
      # Create an AWK script for efficient reference checking
      cat > check_ref.awk << 'EOF'
      BEGIN {
          FS="\t";
          mismatch_count = 0;
          total_checked = 0;
      }
      {
          chrom = $1;
          pos = $2;
          ref_allele = $4;
          
          # Skip if allele is not simple nucleotides
          if (ref_allele !~ /^[ACGTN]+$/) {
              print "Skipping variant with non-standard REF allele: " ref_allele;
              next;
          }
          
          # Get correct chromosome name for reference lookup
          query_chrom = getRefChrom(chrom);
          
          # Execute samtools command to get reference base
          cmd = "samtools faidx ${reference} " query_chrom ":" pos "-" (pos + length(ref_allele) - 1) " | grep -v '>' | tr -d '\\n'";
          cmd | getline ref_base;
          close(cmd);
          
          print "Checking variant " chrom ":" pos " (querying as " query_chrom ")";
          print "  VCF REF allele: " ref_allele;
          print "  Genome sequence: " ref_base;
          
          total_checked++;
          
          if (ref_allele != ref_base) {
              print "WARNING: Reference mismatch at " chrom ":" pos " (querying as " query_chrom ":" pos ")";
              print "  VCF REF allele: " ref_allele;
              print "  Reference genome: " ref_base;
              mismatch_count++;
          }
      }
      END {
          print "REF_MISMATCHES=" mismatch_count;
          print "TOTAL_CHECKED=" total_checked;
          
          # Exit with error code if too many mismatches
          if (mismatch_count > 0) {
              exit 1;
          } else {
              exit 0;
          }
      }
EOF
      
      # Run the reference check with the previous chromosome mapping logic
      awk -f chr_mapping.awk -f check_ref.awk sample_variants.txt > ref_check_results.txt
      REF_CHECK_STATUS=$?
      
      # Extract the mismatch counts
      REF_MISMATCHES=\$(grep "REF_MISMATCHES=" ref_check_results.txt | cut -d= -f2)
      TOTAL_CHECKED=\$(grep "TOTAL_CHECKED=" ref_check_results.txt | cut -d= -f2)
      
      # Display the results
      cat ref_check_results.txt | grep -v "REF_MISMATCHES=" | grep -v "TOTAL_CHECKED="
      
      # Report the results
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