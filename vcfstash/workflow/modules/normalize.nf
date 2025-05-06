process RenameAndNormalizeVCF {
    scratch = true
    debug = true

    input:
    path chr_add_file
    path vcf_file
    path reference_genome
    val sample_name
    val do_sort
    val remove_gt
    val remove_info

    output:
    path "${sample_name}_norm.bcf", emit: norm_bcf
    path "${sample_name}_norm.bcf.csi", emit: norm_bcf_index

    script:
    // -m ${params.bcftools_sort_memory}
    def sort_command = do_sort ? "${params.bcftools_cmd} sort -T ${task.scratch} |" : ""
    def gt_option = remove_gt ? "-G" : ""
    // Only remove specific INFO tag during annotation, otherwise remove all INFO
	def info_option = remove_info ? "-x INFO" : (
		// For non-stash modes, must_contain_info_tag is required
		!params.must_contain_info_tag ?
		{ error "must_contain_info_tag parameter is required and cannot be empty for annotation mode" } :
		"-x INFO/${params.must_contain_info_tag}"
	)

    """

    # First, convert to BCF and index to ensure proper format
    ${params.bcftools_cmd} view ${gt_option} -Ou ${vcf_file} --threads ${(task.cpus).intdiv(3)} |
    ${params.bcftools_cmd} annotate ${info_option} --rename-chrs ${chr_add_file} --threads ${(task.cpus).intdiv(3)} -Ou |
    ${params.bcftools_cmd} filter -i 'CHROM ~ "^chr[1-9,X,Y,M]" && CHROM ~ "[0-9,X,Y,M]\$"' --threads ${(task.cpus).intdiv(3)} -Ou | \\
    ${sort_command} \\
    ${params.bcftools_cmd} norm -m- -c x -f ${reference_genome} -o ${sample_name}_norm.bcf -Ob --threads ${(task.cpus).intdiv(3)} --write-index
    """
}

workflow NORMALIZE {
    take:
    chr_add
    vcf
    reference
    sample_name
    do_sort
    remove_gt
    remove_info

    main:
    RenameAndNormalizeVCF(chr_add, vcf, reference, sample_name, do_sort, remove_gt, remove_info)

    emit:
    norm_bcf = RenameAndNormalizeVCF.out.norm_bcf
    norm_bcf_index = RenameAndNormalizeVCF.out.norm_bcf_index
}