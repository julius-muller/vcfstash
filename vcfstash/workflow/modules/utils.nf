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
	path chr_add

	output:
	val true, emit: validated

	script:
	"""
	# Check file existence
	test -e "${vcf}" || { echo "VCF file not found: ${vcf}"; exit 1; }
	test -e "${vcf_index}" || { echo "VCF index file not found: ${vcf_index}"; exit 1; }
	test -e "${chr_add}" || { echo "Chr add file not found: ${chr_add}"; exit 1; }

      """
  }


workflow UTILS {
    take:
    sample_name
    output_dir
    vcf
    chr_add

    main:
    // Explicitly look for VCF index file - try both csi and tbi formats
    vcf_index = file("${vcf}.{csi,tbi}", checkIfExists: true)[0]

    validateInputResult = ValidateInputs(
        vcf,
        vcf_index,
        chr_add,
    )

    emit:
    validate = validateInputResult.validated
}
