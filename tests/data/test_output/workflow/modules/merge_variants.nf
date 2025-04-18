process MergeVariants {
    input:
    path existing_bcf
    path existing_bcf_index
    path new_bcf
    path new_bcf_index

    output:
    path "merged.bcf", emit: merged_bcf
    path "merged.bcf.csi", emit: merged_bcf_index

    script:
    """
    # Merge variants
    ${params.bcftools_cmd} concat \
        --allow-overlaps \
        --rm-dup all \
        -Ob \
        --write-index \
        -o merged.bcf \
        -G \
        --threads ${task.cpus} \
        ${existing_bcf} \
        ${new_bcf}
    """
}

workflow MERGE_VARIANTS {
    take:
    existing_bcf
    existing_bcf_index
    new_bcf
    new_bcf_index

    main:
    MergeVariants(existing_bcf, existing_bcf_index, new_bcf, new_bcf_index)

    emit:
    merged_bcf = MergeVariants.out.merged_bcf
    merged_bcf_index = MergeVariants.out.merged_bcf_index
}