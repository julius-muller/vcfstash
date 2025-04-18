// workflow/modules/annotate_cache.nf

process AnnotateFromDB {
    scratch true

    input:
    path input_bcf
    path input_bcf_index
    path db_bcf
    path db_bcf_index
    val sample_name

    output:
    path "${sample_name}_isecvst.bcf", emit: annotated_bcf
    path "${sample_name}_isecvst.bcf.csi", emit: annotated_bcf_index

    script:
    """
    ${params.bcftools_cmd} annotate -a ${db_bcf} ${input_bcf} -c INFO -o ${sample_name}_isecvst.bcf -Ob --threads ${task.cpus} -W
    """
}

process FilterMissingAnnotations {
    scratch true

    input:
    path annotated_bcf
    path annotated_bcf_index
    val sample_name

    output:
    path "${sample_name}_isecvst_miss.bcf", emit: filtered_bcf
    path "${sample_name}_isecvst_miss.bcf.csi", emit: filtered_bcf_index

    script:
    """
    ${params.bcftools_cmd} filter -i 'INFO/CSQ==""' -Ob -o ${sample_name}_isecvst_miss.bcf ${annotated_bcf} --threads ${task.cpus} -W
    """
}

process MergeAnnotations {
    scratch true

    input:
    path original_annotated_bcf
    path original_annotated_bcf_index
    path newly_annotated_bcf
    path newly_annotated_bcf_index
    val sample_name

    output:
    path "${sample_name}_vst.bcf", emit: merged_bcf
    path "${sample_name}_vst.bcf.csi", emit: merged_bcf_index

    script:
    """
    ${params.bcftools_cmd} annotate -a ${newly_annotated_bcf} ${original_annotated_bcf} -c INFO -o ${sample_name}_vst.bcf -Ob --threads ${task.cpus} -W
    """
}