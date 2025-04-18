process Annotate {
    debug true

    input:
    path input_bcf
    path input_bcf_index

    output:
    path "vcfstash_annotated.bcf", emit: anno_bcf
    path "vcfstash_annotated.bcf.csi", emit: anno_bcf_index

    script:
    """
    export INPUT_BCF="${input_bcf}"
    export OUTPUT_BCF="vcfstash_annotated.bcf"


    echo "Checking tool version..."
    TOOL_VERSION_OUTPUT=\$(${params.annotation_tool_cmd} ${params.tool_version_command})
	TOOL_VERSION=\$(echo "\$TOOL_VERSION_OUTPUT" | ${params.tool_version_regex})
    if [ -z "\$TOOL_VERSION" ]; then
        echo "Error: Unable to determine tool version. Output: \$TOOL_VERSION_OUTPUT" >&2
        exit 1
    fi
    echo "Tool version detected: \$TOOL_VERSION"

	if [ "\$TOOL_VERSION" != "${params.required_tool_version}" ]; then
		echo "Error: Tool version mismatch. Found \$TOOL_VERSION but required ${params.required_tool_version}" >&2
		exit 1
	fi

    echo "Executing command: ${params.annotation_cmd}"
    ${params.annotation_cmd}

    # # Check output files and report variant counts
    if [ ! -f "\$OUTPUT_BCF" ]; then
        echo "Error: Output file vcfstash_annotated.bcf was not created" >&2
        exit 1
    fi

    input_variants=\$(bcftools index -n "\$INPUT_BCF".csi)
    output_variants=\$(bcftools index -n "\$OUTPUT_BCF".csi)
    echo "Annotation complete: Input variants: \${input_variants}, Output variants: \${output_variants}"

    """
}

workflow ANNOTATE {
    take:
    input_bcf
    input_bcf_index

    main:
    annotation_result = Annotate(
        input_bcf,
        input_bcf_index
    )

    emit:
    annotated_bcf = annotation_result.anno_bcf
    annotated_bcf_index = annotation_result.anno_bcf_index
}