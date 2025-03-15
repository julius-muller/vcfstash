nextflow.enable.dsl = 2

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

process PerformIntersection {
    scratch = true

    input:
    path norm_bcf
    path norm_bcf_index
    path db_bcf
    path db_bcf_index

    output:
    path "out", emit: isec_dir

    script:
    """
    bcftools isec $norm_bcf $db_bcf --threads ${task.cpus} -p out -O b
    """
}

process AnnotateIntersection {
    scratch = true

    input:
    path isec_dir
    val sample_name

    output:
    path "${sample_name}_norm_vep.bcf", emit: annotated_bcf
    path "${sample_name}_norm_vep.bcf.csi", emit: annotated_bcf_index

    script:
    """
    bcftools annotate -a ${isec_dir}/0003.bcf \
      -c INFO/CSQ,INFO/gnomadg_ac,INFO/gnomadg_an,INFO/gnomadg_af,INFO/gnomadg_filter,INFO/gnomade_ac,INFO/gnomade_an,INFO/gnomade_af,INFO/gnomade_filter,INFO/clinvar_clnsig,INFO/clinvar_clnrevstat,INFO/clinvar_clndn \
      -o ${sample_name}_norm_vep.bcf ${isec_dir}/0000.bcf -O b --threads ${task.cpus} --write-index
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
    bcftools isec $norm_bcf $db_bcf -n=1 -w1 --threads ${task.cpus} -o ${sample_name}_norm_sub.bcf --write-index -O b
    """
}

process EchtvarAnnotate {
    scratch = true

    input:
    path subset_bcf
    val sample_name

    output:
    path "${sample_name}_norm_sub_ev.bcf", emit: ev_bcf
    path "${sample_name}_norm_sub_ev.bcf.csi", emit: ev_bcf_index

    script:
    """
    /mnt/data/apps/echtvar/echtvar anno -e /mnt/data/resources/gnomad/vep_gnomad_v4_hg19_genomes/echtvar_vep_gnomad4.1_genome_hg19lo.zip \
                 -e /mnt/data/resources/gnomad/vep_gnomad_v4_hg19_exomes/echtvar_vep_gnomad4.1_exome_hg19lo.zip \
                 -e /mnt/data/resources/clinvar/vcf_GRCh37/clinvar_20241111_hg19.zip \
                 $subset_bcf ${sample_name}_norm_sub_ev.bcf
    bcftools index ${sample_name}_norm_sub_ev.bcf --threads ${task.cpus}
    """
}

process EchtvarAnnotateChr {
    scratch = true
    memory '60 GB'

    input:
    path subset_bcf
    path subset_bcf_index
    val sample_name

    output:
    path "${sample_name}_norm_sub_ev.bcf", emit: ev_bcf
    path "${sample_name}_norm_sub_ev.bcf.csi", emit: ev_bcf_index

    script:
    """

        # Debug information
    echo "Current directory contents:"
    ls -la
    echo "BCF file exists: \$(test -f ${subset_bcf} && echo 'yes' || echo 'no')"
    echo "BCF index exists: \$(test -f ${subset_bcf_index} && echo 'yes' || echo 'no')"


    # Split by chromosome first
    mkdir -p split_chr
    for chr in \$(bcftools index -s ${subset_bcf} | cut -f1); do
        echo "Processing chromosome \$chr"
        bcftools view ${subset_bcf} \$chr -Ob -o split_chr/\${chr}.bcf --threads ${task.cpus} --write-index
    done

    # List files for debugging
    echo "Files in split_chr directory:"
    ls -l split_chr/

    # Process each chromosome with explicit file handling
    for chr_file in split_chr/*.bcf; do
        if [ -f "\$chr_file" ]; then
            echo "Annotating \$chr_file"
            /mnt/data/apps/echtvar/echtvar anno \
                -e /mnt/data/resources/gnomad/vep_gnomad_v4_hg19_genomes/echtvar_vep_gnomad4.1_genome_hg19lo.zip \
                -e /mnt/data/resources/gnomad/vep_gnomad_v4_hg19_exomes/echtvar_vep_gnomad4.1_exome_hg19lo.zip \
                -e /mnt/data/resources/clinvar/vcf_GRCh37/clinvar_20241111_hg19.zip \
                "\$chr_file" "\${chr_file%.bcf}_ev.bcf"
        fi
    done

    # Merge results
    bcftools concat split_chr/*_ev.bcf -Ob -o ${sample_name}_norm_sub_ev.bcf --write-index --threads ${task.cpus}
    """
}

process VEPAnnotate {
    scratch = true
    debug true

    input:
    path input_bcf
    path input_bcf_index
    path vep_cache_dir
    path reference_genome
    val sample_name

    output:
    path "${sample_name}_vep.bcf", emit: vep_bcf
    path "${sample_name}_vep.bcf.csi", emit: vep_bcf_index

    script:
    """
    # Get list of chromosomes
    mkdir -p chr_tmp
    chromosomes=(\$(bcftools index -s ${input_bcf} | cut -f1))

    # Process each chromosome in parallel, limited by max_chr_parallel
    echo "Processing \${#chromosomes[@]} chromosomes..."
    for chr in "\${chromosomes[@]}"; do
        while [ \$(jobs -p | wc -l) -ge ${params.vep_max_chr_parallel} ]; do
            sleep 5
        done

        (
            echo "Starting chromosome \$chr"
            bcftools view -r \$chr ${input_bcf} | \
            apptainer exec -B /mnt/data:/mnt/data /mnt/data/apps/ensembl-vep/113/vep.sif vep \
                -a GRCh37 \
                --cache \
                --offline \
                --force_overwrite \
                --dir_cache ${vep_cache_dir} \
                --fa ${reference_genome} \
                --format vcf \
                --transcript_version \
                --total_length \
                --flag_pick \
                --exclude_predicted \
                --no_stats \
                --buffer_size 500000 \
                --hgvs \
                --hgvsg \
                --spdi \
                --variant_class \
                --uniprot \
                --transcript_version \
                --gene_version \
                --protein \
                --symbol \
                --canonical \
                --appris \
                --mane \
                --biotype \
                --domains \
                --gencode_basic \
                --fork ${params.vep_max_forks} \
                --vcf \
                -i stdin \
                -o stdout | \
            bcftools view -Ob -o chr_tmp/\${chr}_vep.bcf --write-index
            echo "Finished chromosome \$chr"
        ) &
    done

    # Wait for all chromosomes to finish
    wait
    echo "All chromosomes processed"

    # Merge chromosome results
    echo "Merging chromosomes..."
    bcftools concat chr_tmp/*_vep.bcf -Ob -o ${sample_name}_vep.bcf --write-index --threads ${task.cpus}
    """
}

// In MergeAndSort process
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

// Capture versions of major tools
process CaptureToolVersions {
    publishDir "${output_dir}", mode: 'copy'

    input:
    val sample_name
    path output_dir

    output:
    path "${sample_name}_final.tool_version.log"

    script:
    """
    echo "bcftools version:" > ${sample_name}_final.tool_version.log
    bcftools --version >> ${sample_name}_final.tool_version.log
    echo "vep version:" >> ${sample_name}_final.tool_version.log
    /usr/bin/apptainer exec -B /mnt/data:/mnt/data /mnt/data/apps/ensembl-vep/113/vep.sif vep 2>&1 | grep "ensembl-.*" >> ${sample_name}_final.tool_version.log
    echo "echtvar version:" >> ${sample_name}_final.tool_version.log
    /mnt/data/apps/echtvar/echtvar --version >> ${sample_name}_final.tool_version.log
    """

    stub:
    """
    touch ${sample_name}_final.tool_version.log
    """
}

workflow {
    // Extract sample name from the input file
    if (!params.input) {
        error "Input file parameter is required: --input"
    }
    def sampleName = params.input.tokenize('/')[-1].replace('.vcf', '').replace('.bcf', '')

    // Define output directory
    if (!params.output) {
        error "Output directory parameter is required: --output"
    }
    def outputDir = file(params.output)
    def projectPath = new File(workflow.projectDir.toString(), "../..").getCanonicalPath()

    // Input files from config with validation
    if (!params.chr_add) {
        error "Chr add file parameter is required: --chr_add"
    }
    if (!params.reference) {
        error "Reference genome parameter is required: --reference"
    }
    if (!params.vep_cache) {
        error "VEP cache directory parameter is required: --vep_cache"
    }

    chr_add = file(params.chr_add)
    vcf = file(params.input)
    vcf_index = file("${params.input}.csi")
    vep_cache = file(params.vep_cache)
    reference = file(params.reference)

    if (params.db_mode) {
        // Database annotation mode - skip normalization and go directly to annotation
        EchtvarAnnotateChr(vcf, vcf_index, sampleName)
        VEPAnnotate(EchtvarAnnotateChr.out.ev_bcf, EchtvarAnnotateChr.out.ev_bcf_index, vep_cache, reference, sampleName)

        // Publish final database BCF
        VEPAnnotate.out.vep_bcf.subscribe { bcf ->
            file(bcf).copyTo("${outputDir}/${sampleName}.bcf")
        }
        VEPAnnotate.out.vep_bcf_index.subscribe { idx ->
            file(idx).copyTo("${outputDir}/${sampleName}.bcf.csi")
        }
    } else if (params.db_bcf == false) {
        // Direct VEP annotation without database
        RenameAndNormalizeVCF(chr_add, vcf, reference, sampleName)
        EchtvarAnnotate(RenameAndNormalizeVCF.out.norm_bcf, sampleName)
        VEPAnnotate(EchtvarAnnotate.out.ev_bcf, EchtvarAnnotate.out.ev_bcf_index, vep_cache, reference, sampleName)

        // Modified merge step to avoid duplicate inputs
        MergeAndSort(
            VEPAnnotate.out.vep_bcf,
            VEPAnnotate.out.vep_bcf_index,
            EchtvarAnnotate.out.ev_bcf,
            EchtvarAnnotate.out.ev_bcf_index,
            outputDir,
            sampleName
        ) | \
        publishDir(params.output, mode: 'copy', saveAs: { filename ->
            filename.replace('_vep.bcf', '_norm_final.bcf')
        })
        // CreateParquet(MergeAndSort.out.final_bcf, MergeAndSort.out.final_bcf_index, sampleName, outputDir, projectPath)
    } else {
        // Sample annotation mode - full workflow with normalization
        RenameAndNormalizeVCF(chr_add, vcf, reference, sampleName)
        db_bcf = file(params.db_bcf)
        db_bcf_index = file("${params.db_bcf}.csi")

        PerformIntersection(RenameAndNormalizeVCF.out.norm_bcf, RenameAndNormalizeVCF.out.norm_bcf_index, db_bcf, db_bcf_index)
        AnnotateIntersection(PerformIntersection.out.isec_dir, sampleName)
        SubsetIntersection(RenameAndNormalizeVCF.out.norm_bcf, RenameAndNormalizeVCF.out.norm_bcf_index, db_bcf, db_bcf_index, sampleName)
        EchtvarAnnotate(SubsetIntersection.out.subset_bcf, sampleName)
        VEPAnnotate(EchtvarAnnotate.out.ev_bcf, EchtvarAnnotate.out.ev_bcf_index, vep_cache, reference, sampleName)
        MergeAndSort(
        AnnotateIntersection.out.annotated_bcf,
        AnnotateIntersection.out.annotated_bcf_index,
        VEPAnnotate.out.vep_bcf,
         VEPAnnotate.out.vep_bcf_index,
         outputDir,
         sampleName)
        // CreateParquet(MergeAndSort.out.final_bcf, MergeAndSort.out.final_bcf_index, sampleName, outputDir, projectPath)
    }

    // Always capture tool versions
    CaptureToolVersions(sampleName, outputDir)
}

