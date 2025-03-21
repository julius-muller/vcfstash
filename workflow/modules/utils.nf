
process CaptureToolVersions {
    publishDir "${output_dir}", mode: 'copy'

    input:
    val sample_name
    path output_dir

    output:
    path "${sample_name}_final.tool_version.log", emit: version_log

    script:
    """
    echo "bcftools version:" > ${sample_name}_final.tool_version.log
    bcftools --version >> ${sample_name}_final.tool_version.log
    echo "vep version:" >> ${sample_name}_final.tool_version.log
    ${params.container_engine} exec ${params.container_bind} ${params.vep_container} vep 2>&1 | grep "ensembl-.*" >> ${sample_name}_final.tool_version.log
    echo "echtvar version:" >> ${sample_name}_final.tool_version.log
    ${params.echtvar_path} --version >> ${sample_name}_final.tool_version.log
    """

    stub:
    """
    touch ${sample_name}_final.tool_version.log
    """
}

process ValidateInputs {
    input:
    path input_file
    path chr_add
    path reference
    path vep_cache

    output:
    val true, emit: validated

    script:
    """
    # Check if all required files exist
    [ -f "${input_file}" ] || (echo "Input file not found: ${input_file}" && exit 1)
    [ -f "${chr_add}" ] || (echo "Chr add file not found: ${chr_add}" && exit 1)
    [ -f "${reference}" ] || (echo "Reference genome not found: ${reference}" && exit 1)
    [ -d "${vep_cache}" ] || (echo "VEP cache directory not found: ${vep_cache}" && exit 1)

    # Check if input file is a valid VCF/BCF
    bcftools view -h ${input_file} > /dev/null || (echo "Invalid VCF/BCF file: ${input_file}" && exit 1)
    """
}

workflow UTILS {
    take:
    sample_name
    output_dir
    input_file
    chr_add
    reference
    vep_cache

    main:
    ValidateInputs(input_file, chr_add, reference, vep_cache)
    CaptureToolVersions(sample_name, output_dir)

    emit:
    version_log = CaptureToolVersions.out.version_log
}