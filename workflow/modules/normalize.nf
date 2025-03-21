
process RenameAndNormalizeVCF {
    scratch = true

    input:
    path chr_add_file
    path vcf_file
    path reference_genome
    val sample_name

    output:
    path "${sample_name}_norm.bcf", emit: norm_bcf
    path "${sample_name}_norm.bcf.csi", emit: norm_bcf_index

    script:
    """
    # First, convert to BCF and index to ensure proper format
    bcftools annotate -x INFO --rename-chrs ${chr_add_file} ${vcf_file} --threads ${(task.cpus).intdiv(2)} -Ob --write-index -o tmp1.bcf
    bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY tmp1.bcf --threads ${(task.cpus).intdiv(2)} | \
    bcftools sort -m 20G -T ${task.scratch} | \
    bcftools norm -m- -c x -f ${reference_genome} -o ${sample_name}_norm.bcf -O b --threads ${(task.cpus).intdiv(2)} --write-index

    """
}

workflow NORMALIZE {
    take:
    chr_add
    vcf
    reference
    sample_name

    main:
    RenameAndNormalizeVCF(chr_add, vcf, reference, sample_name)

    emit:
    norm_bcf = RenameAndNormalizeVCF.out.norm_bcf
    norm_bcf_index = RenameAndNormalizeVCF.out.norm_bcf_index
}