process RunAnnotation {
    debug true
    scratch false

    input:
    path input_bcf
    path input_bcf_index

    output:
    path "vcfstash_annotated.bcf", emit: annotated_bcf
    path "vcfstash_annotated.bcf.csi", emit: annotated_bcf_index
    path "vcfstash_annotated.log", emit: annotated_bcf_log
    path "auxiliary/*", emit: aux_files

    script:
    """
	# Execute the annotation command with strict error handling
	set -euo pipefail

    echo "[`date`] VCFSTASH ANNOTATION PROCESS STARTED" | tee vcfstash_annotated.log
    echo "[`date`] ==========================================" | tee -a vcfstash_annotated.log
    echo "[`date`] Task ID: ${task.process} (${task.index})" | tee -a vcfstash_annotated.log
    echo "[`date`] Working directory: ${task.workDir}" | tee -a vcfstash_annotated.log
    echo "[`date`] Input BCF file: ${input_bcf}" | tee -a vcfstash_annotated.log
    echo "[`date`] Input BCF index: ${input_bcf_index}" | tee -a vcfstash_annotated.log
    echo "[`date`] Allocated resources:" | tee -a vcfstash_annotated.log
    echo "[`date`]   - CPUs: ${task.cpus}" | tee -a vcfstash_annotated.log
    echo "[`date`]   - Memory: ${task.memory}" | tee -a vcfstash_annotated.log
    echo "[`date`]   - Time: ${task.time}" | tee -a vcfstash_annotated.log
    echo "[`date`] ==========================================" | tee -a vcfstash_annotated.log

	mkdir -p auxiliary

    # Define key variables locally
    export INPUT_BCF="${input_bcf}"
    export OUTPUT_BCF="vcfstash_annotated.bcf"
    export AUXILIARY_DIR="\$PWD/auxiliary"

    # Create the annotation command from params
    CMD="${params.annotation_cmd}"

    echo "Checking tool version..."
    echo "Tool version command: ${params.tool_version_command}"
    TOOL_VERSION=`${params.tool_version_command}`
    if [ -z "\$TOOL_VERSION" ]; then
        echo "Error: Unable to determine tool version. Output: \$TOOL_VERSION" >&2
        exit 1
    fi
    echo "Tool version detected: \$TOOL_VERSION"

	if [[ "\$TOOL_VERSION" != "${params.required_tool_version}"* ]]; then
		echo "[`date`] ERROR: Tool version mismatch. Found \$TOOL_VERSION but required ${params.required_tool_version}" | tee -a vcfstash_annotated.log
		exit 1
	fi


	{
		echo "[`date`] Running annotation command: \$CMD"
		eval "\$CMD"
	} 2>&1 | tee -a vcfstash_annotated.log

	CMD_EXIT=\${PIPESTATUS[1]:-0}
	if [ "\$CMD_EXIT" -ne 0 ]; then
		echo "[`date`] ERROR: Annotation command failed with exit code \$CMD_EXIT" | tee -a vcfstash_annotated.log
		exit \$CMD_EXIT
	fi

    # Check for success
    if [ ! -f "\${OUTPUT_BCF}" ]; then
        echo "[`date`] ERROR: Annotation failed! Output file not created." | tee -a vcfstash_annotated.log
        exit 1
    fi

    # Check for required INFO tag
    if [ ! -z "${params.must_contain_info_tag}" ]; then
        echo "[`date`] Checking for required INFO tag: ${params.must_contain_info_tag}" | tee -a vcfstash_annotated.log

        HAS_TAG=`${params.bcftools_cmd} view -h \${OUTPUT_BCF} | grep "##INFO=<ID=${params.must_contain_info_tag},"`

        if [ -z "\$HAS_TAG" ]; then
            echo "[`date`] ERROR: Required INFO tag ${params.must_contain_info_tag} not found in output BCF" | tee -a vcfstash_annotated.log
            exit 1
        fi
    fi

    # Calculate and report variant counts
    INPUT_COUNT=`${params.bcftools_cmd} index -n \${INPUT_BCF}`
    OUTPUT_COUNT=`${params.bcftools_cmd} index -n \${OUTPUT_BCF}`
    echo "[`date`] Annotation complete: Input variants: \${INPUT_COUNT}, Output variants: \${OUTPUT_COUNT}" | tee -a vcfstash_annotated.log
    """
}


workflow ANNOTATE {
    take:
    vcf
    vcf_index

    main:
    RunAnnotation(vcf, vcf_index)

    emit:
    annotated_bcf = RunAnnotation.out.annotated_bcf
    annotated_bcf_index = RunAnnotation.out.annotated_bcf_index
    annotated_bcf_log = RunAnnotation.out.annotated_bcf_log
    aux_files = RunAnnotation.out.aux_files
}
