process Annotate {
    debug true
    cpus params.cpus
    memory params.memory

    input:
    path input_bcf
    path input_bcf_index

    output:
    path "annotated.bcf", emit: anno_bcf
    path "annotated.bcf.csi", emit: anno_bcf_index

    script:
    """
    export INPUT_BCF="${input_bcf}"
    export OUTPUT_BCF="annotated.bcf"

    ${params.annotation_cmd}
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