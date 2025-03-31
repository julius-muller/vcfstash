nextflow.enable.dsl = 2

include { ANNOTATE } from './annotate'

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
    bcftools annotate -a ${db_bcf} ${input_bcf} -c INFO -o ${sample_name}_isecvst.bcf -Ob --threads ${task.cpus} -W
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
    bcftools filter -i 'INFO/CSQ==""' -Ob -o ${sample_name}_isecvst_miss.bcf ${annotated_bcf} --threads ${task.cpus} -W
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
    bcftools annotate -a ${newly_annotated_bcf} ${original_annotated_bcf} -c INFO -o ${sample_name}_vst.bcf -Ob --threads ${task.cpus} -W
    """
}

workflow ANNOTATE_INFO {
    take:
    input_bcf
    input_bcf_index
    db_bcf
    db_bcf_index
    vep_cache
    reference
    sample_name

    main:
    // Step 1: Annotate input with database INFO fields
    annotateFromDBResult = AnnotateFromDB(
        input_bcf,
        input_bcf_index,
        db_bcf,
        db_bcf_index,
        sample_name
    )

    // Step 2: Filter to get variants with missing CSQ annotations
    filterMissingResult = FilterMissingAnnotations(
        annotateFromDBResult.annotated_bcf,
        annotateFromDBResult.annotated_bcf_index,
        sample_name
    )

    // Step 3: Run standard ANNOTATE workflow on missing annotations
    annotateResult = ANNOTATE(
        filterMissingResult.filtered_bcf,
        filterMissingResult.filtered_bcf_index,
        vep_cache,
        reference,
        "${sample_name}_isecvst_miss"
    )

    // Step 4: Merge newly annotated variants back into original file
    mergeResult = MergeAnnotations(
        annotateFromDBResult.annotated_bcf,
        annotateFromDBResult.annotated_bcf_index,
        annotateResult.vep_bcf,
        annotateResult.vep_bcf_index,
        sample_name
    )

    emit:
    merged_bcf = mergeResult.merged_bcf
    merged_bcf_index = mergeResult.merged_bcf_index
}