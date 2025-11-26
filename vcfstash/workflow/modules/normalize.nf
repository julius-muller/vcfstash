
process RemoveGTAndInfo {
    scratch = false
    debug = true

    input:
    path vcf_file
    path vcf_index  // Explicitly require the index file as input
    val sample_name
    val remove_gt
    val remove_info

    output:
    path "${sample_name}_filtered.bcf", emit: filtered_bcf
    path "${sample_name}_filtered.bcf.csi", emit: filtered_bcf_index
    path "${sample_name}_filtered.log", emit: filtered_bcf_log

    script:
    def gt_option = remove_gt ? "-G" : ""
    // Only remove specific INFO tag during annotation, otherwise remove all INFO
    def info_option = remove_info ? "-x INFO" : (
        // For non-stash modes, must_contain_info_tag is required
        !params.must_contain_info_tag ?
        { error "must_contain_info_tag parameter is required and cannot be empty for annotation mode" } :
        "-x INFO/${params.must_contain_info_tag}"
    )

    """
    echo "Removing GT and INFO fields..." | tee ${sample_name}_filtered.log

    # Remove GT and INFO fields
    ${params.bcftools_cmd} view ${gt_option} -Ou ${vcf_file} --threads ${(task.cpus).intdiv(3)} |
    ${params.bcftools_cmd} annotate ${info_option} -o ${sample_name}_filtered.bcf -Ob --threads ${(task.cpus).intdiv(3)} --write-index 2>&1 | tee -a ${sample_name}_filtered.log

    # Check if filtering was successful
    if [ ! -f "${sample_name}_filtered.bcf" ] || [ ! -f "${sample_name}_filtered.bcf.csi" ]; then
        echo "ERROR: Filtering failed! Output file not created." | tee -a ${sample_name}_filtered.log
        exit 1
    fi

    input_variants=\$(${params.bcftools_cmd} index -n "${vcf_file}".csi)
    output_variants=\$(${params.bcftools_cmd} index -n "${sample_name}_filtered.bcf".csi)

    echo "Filtering complete: Input variants: \${input_variants}, Output variants: \${output_variants}" | tee -a ${sample_name}_filtered.log
    """
}

process RenameAndNormalizeVCF {
    scratch = false
    debug = true

    input:
    path chr_add_file
    path vcf_file
    path vcf_index  // Explicitly require the index file as input
    val sample_name
    val remove_gt
    val remove_info

    output:
    path "${sample_name}_norm.bcf", emit: norm_bcf
    path "${sample_name}_norm.bcf.csi", emit: norm_bcf_index
    path "${sample_name}_norm.log", emit: norm_bcf_log

    script:
    def gt_option = remove_gt ? "-G" : ""
    // Only remove specific INFO tag during annotation, otherwise remove all INFO
    def info_option = remove_info ? "-x INFO" : (
        // For non-stash modes, must_contain_info_tag is required
        !params.must_contain_info_tag ?
        { error "must_contain_info_tag parameter is required and cannot be empty for annotation mode" } :
        "-x INFO/${params.must_contain_info_tag}"
    )

    """
    echo "Proceeding with normalization..." | tee ${sample_name}_norm.log

    # Now perform the normalization with more confidence
    ${params.bcftools_cmd} view ${gt_option} -Ou ${vcf_file} --threads ${(task.cpus).intdiv(3)} |
    ${params.bcftools_cmd} annotate ${info_option} --rename-chrs ${chr_add_file} --threads ${(task.cpus).intdiv(3)} -Ou |
    ${params.bcftools_cmd} filter -i 'CHROM ~ "^chr[1-9,X,Y,M]" && CHROM ~ "[0-9,X,Y,M]\$"' --threads ${(task.cpus).intdiv(3)} -Ou |
    ${params.bcftools_cmd} norm -m- -o ${sample_name}_norm.bcf -Ob --threads ${(task.cpus).intdiv(3)} --write-index 2>&1 | tee -a ${sample_name}_norm.log

    # Check if normalization was successful
    if [ ! -f "${sample_name}_norm.bcf" ] || [ ! -f "${sample_name}_norm.bcf.csi" ]; then
        echo "ERROR: Normalization failed! Output file not created." | tee -a ${sample_name}_norm.log
        exit 1
    fi

    input_variants=\$(${params.bcftools_cmd} index -n "${vcf_file}".csi)
    output_variants=\$(${params.bcftools_cmd} index -n "${sample_name}_norm.bcf".csi)

    echo "Normalization complete: Input variants: \${input_variants}, Output variants: \${output_variants}" | tee -a ${sample_name}_norm.log
    """
}

workflow FILTER_ONLY {
    take:
    vcf
    vcf_index  // Add the index as a parameter to the workflow
    sample_name
    remove_gt
    remove_info

    main:
    RemoveGTAndInfo(vcf, vcf_index, sample_name, remove_gt, remove_info)

    emit:
    filtered_bcf = RemoveGTAndInfo.out.filtered_bcf
    filtered_bcf_index = RemoveGTAndInfo.out.filtered_bcf_index
    filtered_bcf_log = RemoveGTAndInfo.out.filtered_bcf_log
}

workflow NORMALIZE {
    take:
    chr_add
    vcf
    vcf_index  // Add the index as a parameter to the workflow
    sample_name
    remove_gt
    remove_info

    main:
    RenameAndNormalizeVCF(chr_add, vcf, vcf_index, sample_name, remove_gt, remove_info)

    emit:
    norm_bcf = RenameAndNormalizeVCF.out.norm_bcf
    norm_bcf_index = RenameAndNormalizeVCF.out.norm_bcf_index
    norm_bcf_log = RenameAndNormalizeVCF.out.norm_bcf_log
}
