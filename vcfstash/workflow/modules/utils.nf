
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
    input:
    path vcf
    path vcf_index
    path reference
    path reference_index
    path chr_add

    output:
    val true, emit: validated

    script:
    """
    # Check if files exist
    [ -f "${vcf}" ] || { echo "VCF file not found: ${vcf}"; exit 1; }
	[ -f "${vcf}.csi" -o -f "${vcf}.tbi" ] || { echo "No valid index found for VCF file: ${vcf} (need .csi or .tbi)"; exit 1; }
    [ -f "${reference}" ] || { echo "Reference file not found: ${reference}"; exit 1; }
    [ -f "${reference_index}" ] || { echo "Reference index not found: ${reference_index}"; exit 1; }
    [ -f "${chr_add}" ] || { echo "Chr add file not found: ${chr_add}"; exit 1; }

    # Basic file format checks
    ${params.bcftools_cmd} view -h "${vcf}" >/dev/null || { echo "Invalid VCF/BCF format"; exit 1; }
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
        file("${vcf}.csi"),
        reference,
        reference_index,
        chr_add,
    )

    emit:
    validate = validateInputResult.validated
}