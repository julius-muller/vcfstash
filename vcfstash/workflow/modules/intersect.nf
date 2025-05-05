
process PerformIntersection {
    scratch = true

    input:
    path norm_bcf
    path norm_bcf_index
    path db_bcf
    path db_bcf_index

    output:
    path "out", emit: isec_dir
    //  shared could be: -w 1 -n =2 and then followed by -w 1 -n =1
    script:
    """
    ${params.bcftools_cmd} isec $norm_bcf $db_bcf --threads ${task.cpus} -p out -O b

    """
}

process AnnotateIntersection {
    scratch = true

    input:
    path isec_dir
    val sample_name

    output:
    path "${sample_name}_norm_vcf.bcf", emit: annotated_bcf
    path "${sample_name}_norm_vcf.bcf.csi", emit: annotated_bcf_index

    script:
    """
    ${params.bcftools_cmd} annotate -a ${isec_dir}/0003.bcf \
      -c INFO/${params.must_contain_info_tag},INFO/gnomadg_ac,INFO/gnomadg_an,INFO/gnomadg_af,INFO/gnomadg_filter,INFO/gnomade_ac,INFO/gnomade_an,INFO/gnomade_af,INFO/gnomade_filter,INFO/clinvar_clnsig,INFO/clinvar_clnrevstat,INFO/clinvar_clndn \
      -o ${sample_name}_norm_vcf.bcf ${isec_dir}/0000.bcf -O b --threads ${task.cpus} --write-index
    """
}

process SubsetIntersection {
    scratch = true

    input:
    path norm_bcf
    path norm_bcf_index
    path db_bcf
    path db_bcf_index
    val sample_name

    output:
    path "${sample_name}_norm_sub.bcf", emit: subset_bcf
    path "${sample_name}_norm_sub.bcf.csi", emit: subset_bcf_index

    script:
    """
    ${params.bcftools_cmd} isec $norm_bcf $db_bcf -n=1 -w1 --threads ${task.cpus} -o ${sample_name}_norm_sub.bcf --write-index -O b
    """
}

workflow INTERSECT {
    take:
    norm_bcf
    norm_bcf_index
    db_bcf
    db_bcf_index
    sample_name

    main:
    PerformIntersection(norm_bcf, norm_bcf_index, db_bcf, db_bcf_index)
    AnnotateIntersection(PerformIntersection.out.isec_dir, sample_name)
    SubsetIntersection(norm_bcf, norm_bcf_index, db_bcf, db_bcf_index, sample_name)

    emit:
    isec_dir = PerformIntersection.out.isec_dir
    annotated_bcf = AnnotateIntersection.out.annotated_bcf
    annotated_bcf_index = AnnotateIntersection.out.annotated_bcf_index
    subset_bcf = SubsetIntersection.out.subset_bcf
    subset_bcf_index = SubsetIntersection.out.subset_bcf_index
}