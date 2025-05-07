
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
	stageInMode 'symlink'  // Use symbolic links instead of hard links

	input:
	path vcf
	path reference
	path reference_index
	path chr_add

	output:
	val true, emit: validated

	script:
	"""
	# Check file existence - we can use test -e for symlinks
	test -e "${vcf}" || { echo "VCF file not found or not properly linked: ${vcf}"; exit 1; }
	test -e "${reference}" || { echo "Reference file not found or not properly linked: ${reference}"; exit 1; }
	test -e "${reference_index}" || { echo "Reference index not found or not properly linked: ${reference_index}"; exit 1; }
	test -e "${chr_add}" || { echo "Chr add file not found or not properly linked: ${chr_add}"; exit 1; }

	# Test for VCF index files
	test -f "${vcf}.csi" || test -f "${vcf}.tbi" || { echo "No valid index found for VCF file: ${vcf} (need .csi or .tbi)"; exit 1; }

	# Basic header check - just validate without processing entire file
	${params.bcftools_cmd} view -h "${vcf}" | head -n 1 >/dev/null || { echo "Invalid VCF/BCF format"; exit 1; }

	# ===========================================
	# Reference genome compatibility check
	# ===========================================
	echo "Validating reference genome compatibility with VCF..."

	# Get a sample variant from input VCF for testing
	FIRST_VAR=\$(${params.bcftools_cmd} view -H ${vcf} | head -n 1)
	if [ -z "\$FIRST_VAR" ]; then
	  echo "Error: Input VCF appears to be empty"
	  exit 1
	fi

	# Extract chromosome, position and reference allele
	CHROM=\$(echo "\$FIRST_VAR" | cut -f 1)
	POS=\$(echo "\$FIRST_VAR" | cut -f 2)
	REF=\$(echo "\$FIRST_VAR" | cut -f 4)

	echo "Testing variant: \$CHROM:\$POS \$REF"

	# Check if reference matches at this position
	REF_BASE=\$(${params.bcftools_cmd} consensus -f ${reference} -r "\$CHROM:\$POS-\$POS" ${vcf} 2>/dev/null | grep -v ">" | tr -d '\n')

	if [ "\$REF_BASE" != "\$REF" ]; then
	  echo "ERROR: Reference genome mismatch detected!"
	  echo "VCF reference allele at \$CHROM:\$POS is \$REF, but reference genome has \$REF_BASE"
	  echo "This indicates that the provided reference genome does not match the VCF file"
	  echo "This would cause a segmentation fault in bcftools norm"
	  exit 1
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

    validateInputResult = ValidateInputs(
        vcf,
        reference,
        reference_index,
        chr_add,
    )

    emit:
    validate = validateInputResult.validated
}