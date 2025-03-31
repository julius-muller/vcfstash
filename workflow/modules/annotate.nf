process VEPAnnotate {
		scratch false
		debug true

		input:
		path input_bcf
		path input_bcf_index
		val vep_cache_dir
		val reference_genome
		val sample_name

		output:
		path "${sample_name}_vep.bcf", emit: vep_bcf
		path "${sample_name}_vep.bcf.csi", emit: vep_bcf_index

		script:
		// Join VEP options from config into a single string
		def vep_opts = params.vep_options.join(' ')
		"""

		# Get list of chromosomes
		mkdir -p chr_tmp
		chromosomes=(\$(bcftools index -s ${input_bcf} | cut -f1))

		# Process each chromosome in parallel
		echo "Processing \${#chromosomes[@]} chromosomes..."
		for chr in "\${chromosomes[@]}"; do
			while [ \$(jobs -p | wc -l) -ge ${params.vep_max_chr_parallel} ]; do
				sleep 5
			done

			(
				echo "Starting chromosome \$chr"
				bcftools view -r \$chr ${input_bcf} | \
				${params.vep_cmd} \
					${vep_opts} \
					--format vcf \
					--cache \
					--dir_cache ${vep_cache_dir} \
					--fa ${reference_genome} \
					--offline \
					--buffer_size ${params.vep_buffer} \
					--fork ${params.vep_max_forks} \
					--vcf \
					--no_stats \
					-i stdin \
					-o stdout | \
				bcftools view -Ob -o chr_tmp/\${chr}_vep.bcf --write-index

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

process EchtvarAnnotate {
	scratch true

	input:
	path subset_bcf
	val sample_name

	output:
	path "${sample_name}_norm_sub_ev.bcf", emit: ev_bcf
	path "${sample_name}_norm_sub_ev.bcf.csi", emit: ev_bcf_index

	script:
	"""
	${params.echtvar_cmd} anno -e ${params.echtvar_gnomad_genome} \
				 -e ${params.echtvar_gnomad_exome} \
				 -e ${params.echtvar_clinvar} \
				 $subset_bcf ${sample_name}_norm_sub_ev.bcf
	bcftools index ${sample_name}_norm_sub_ev.bcf --threads ${task.cpus}
	"""
}

process EchtvarAnnotateChr {
	scratch true
	memory '10 GB'

	input:
	path subset_bcf
	path subset_bcf_index
	val sample_name

	output:
	path "${sample_name}_norm_sub_ev.bcf", emit: ev_bcf
	path "${sample_name}_norm_sub_ev.bcf.csi", emit: ev_bcf_index

	script:
	"""
	# Split by chromosome first
	mkdir -p split_chr
	for chr in \$(bcftools index -s ${subset_bcf} | cut -f1); do
		echo "Processing chromosome \$chr"
		bcftools view ${subset_bcf} \$chr -Ob -o split_chr/\${chr}.bcf --threads 1 --write-index
	done

	# Process each chromosome
	for chr_file in split_chr/*.bcf; do
		if [ -f "\$chr_file" ]; then
			echo "Annotating \$chr_file"
			${params.echtvar_cmd} anno \
				-e ${params.echtvar_gnomad_genome} \
				-e ${params.echtvar_gnomad_exome} \
				-e ${params.echtvar_clinvar} \
				"\$chr_file" "\${chr_file%.bcf}_ev.bcf"
		fi
	done

	# Merge results
	bcftools concat split_chr/*_ev.bcf -Ob -o ${sample_name}_norm_sub_ev.bcf --write-index --threads ${task.cpus}
	"""
}

workflow ANNOTATE {
    take:
    input_bcf
    input_bcf_index
    vep_cache
    reference
    sample_name

    main:
    // Choose the appropriate echtvar process based on db_mode
    if (params.db_mode) {
        echtvar_result = EchtvarAnnotateChr(input_bcf, input_bcf_index, sample_name)
    } else {
        echtvar_result = EchtvarAnnotate(input_bcf, sample_name)
    }

    // Run VEP annotation on the echtvar result
    vep_result = VEPAnnotate(
        echtvar_result.ev_bcf,
        echtvar_result.ev_bcf_index,
        vep_cache,
        reference,
        sample_name
    )

    emit:
    vep_bcf = vep_result.vep_bcf
    vep_bcf_index = vep_result.vep_bcf_index

}