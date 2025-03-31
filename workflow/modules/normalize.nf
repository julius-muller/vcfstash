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
    def sort_command = do_sort ? "bcftools sort -m ${params.bcftools_sort_memory} -T ${task.scratch} |" : ""
    def gt_option = remove_gt ? "-G" : ""
    def info_option = remove_info ? "-x INFO" : ""  // Create annotation option based on parameter

    """
    # First, convert to BCF and index to ensure proper format
    bcftools view ${gt_option} -Ou ${vcf_file} --threads ${(task.cpus).intdiv(3)} |
    bcftools annotate ${info_option} --rename-chrs ${chr_add_file} --threads ${(task.cpus).intdiv(3)} -Ou |
    bcftools filter -i 'CHROM ~ "^chr[1-9,X,Y,M]" && CHROM ~ "[0-9,X,Y,M]\$"' --threads ${(task.cpus).intdiv(3)} -Ou | \\
    ${sort_command} \\
    bcftools norm -m- -c x -f ${reference_genome} -o ${sample_name}_norm.bcf -Ob --threads ${(task.cpus).intdiv(3)} --write-index
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