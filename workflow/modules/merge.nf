
process MergeAndSort {
    scratch = true

    input:
    path annotated_bcf
    path annotated_bcf_index
    path vep_bcf
    path vep_bcf_index
    val output_dir
    val sample_name

    publishDir "${output_dir}", mode: 'copy'

    output:
    tuple path("${sample_name}_norm_final.bcf"), path("${sample_name}_norm_final.bcf.csi"), emit: final_files

    script:
    """
    # Merge and ensure single sample
    bcftools merge ${annotated_bcf} ${vep_bcf} --threads ${task.cpus} | \
    bcftools view -s mgm_WGS_32 -Ob -o "${sample_name}_norm_final.bcf" --threads ${task.cpus}
    bcftools index "${sample_name}_norm_final.bcf" --threads ${task.cpus}
    """
}

process CreateParquet {
    scratch = true

    input:
    path final_bcf
    path final_bcf_index
    val sample_name
    val output_dir
    val project_path

    publishDir "${output_dir}", mode: 'copy'

    output:
    path "${sample_name}_final.parquet", emit: parquet_file

    script:
    """
    ${project_path}/.venv/bin/python ${project_path}/scripts/pq_pl.py -i $final_bcf -o ${sample_name}_final.parquet -t 10
    """
}

workflow MERGE {
    take:
    annotated_bcf
    annotated_bcf_index
    vep_bcf
    vep_bcf_index
    output_dir
    sample_name
    project_path

    main:
    MergeAndSort(annotated_bcf, annotated_bcf_index, vep_bcf, vep_bcf_index, output_dir, sample_name)
    CreateParquet(MergeAndSort.out.final_files[0], MergeAndSort.out.final_files[1], sample_name, output_dir, project_path)

    emit:
    final_bcf = MergeAndSort.out.final_files[0]
    final_bcf_index = MergeAndSort.out.final_files[1]
    parquet_file = CreateParquet.out.parquet_file
}