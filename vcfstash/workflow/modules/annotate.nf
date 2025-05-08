process Annotate {
    debug true
    // Add error strategy to retry on failure
    errorStrategy { task.attempt <= 1 ? 'retry' : 'terminate' }
    maxRetries 1

    // Add resource tracking
    time { params.containsKey('annotation_max_time') ? params.annotation_max_time : '24h' }
    memory { params.containsKey('annotation_memory') ? params.annotation_memory : '16 GB' }

    input:
    path input_bcf
    path input_bcf_index

    output:
    path "vcfstash_annotated.bcf", emit: anno_bcf
    path "vcfstash_annotated.bcf.csi", emit: anno_bcf_index
    path "vcfstash_annotated.log", emit: anno_bcf_log

    script:
    """
    # Create a timestamp function for better logging
    timestamp() {
        date +"[%Y-%m-%d %H:%M:%S]"
    }

    # Setup logging with timestamps and redirections
    {
        echo "\$(timestamp) VCFSTASH ANNOTATION PROCESS STARTED"
        echo "\$(timestamp) =========================================="
        echo "\$(timestamp) Task ID: ${task.process} (${task.index})"
        echo "\$(timestamp) Working directory: ${task.workDir}"
        echo "\$(timestamp) Input BCF file: ${input_bcf}"
        echo "\$(timestamp) Input BCF index: ${input_bcf_index}"
        echo "\$(timestamp) Allocated resources:"
        echo "\$(timestamp)   - CPUs: ${task.cpus}"
        echo "\$(timestamp)   - Memory: ${task.memory}"
        echo "\$(timestamp)   - Time: ${task.time}"
        echo "\$(timestamp) =========================================="

        # Use current directory for all files
        export INPUT_BCF="${input_bcf}"
        export OUTPUT_BCF="vcfstash_annotated.bcf"
        export OUTPUT_DIR="\$PWD"

        echo "\$(timestamp) Using current directory for all output files: \$OUTPUT_DIR"

        # Check input file integrity
        echo "\$(timestamp) Checking input file integrity..."
        if [ ! -f "\$INPUT_BCF" ]; then
            echo "\$(timestamp) ERROR: Input BCF file not found: \$INPUT_BCF"
            exit 1
        fi

        if [ ! -f "\${INPUT_BCF}.csi" ] && [ ! -f "\${INPUT_BCF}.tbi" ]; then
            echo "\$(timestamp) ERROR: Input VCF/BCF index not found: \${INPUT_BCF}.tbi/csi"
            exit 1
        fi

        # Validate the BCF file can be read
        echo "\$(timestamp) Validating input BCF file can be read..."
        HEADER_CHECK=\$(bcftools view -h "\$INPUT_BCF" 2>&1)
        if [ \$? -ne 0 ]; then
            echo "\$(timestamp) ERROR: Cannot read input BCF file. Details:"
            echo "\$HEADER_CHECK"
            exit 1
        fi

        # Get input variant count - redirect stderr to /dev/null to suppress warnings
        input_variants=\$(bcftools index -n "\$INPUT_BCF".csi 2>/dev/null || echo "ERROR_COUNTING")
        if [[ "\$input_variants" == "ERROR_COUNTING" ]]; then
            echo "\$(timestamp) WARNING: Could not count input variants from index. Will try using stats..."
            input_variants=\$(bcftools stats "\$INPUT_BCF" | grep "number of records:" | awk '{print \$4}' || echo "UNKNOWN")
        fi
        echo "\$(timestamp) Input file contains \$input_variants variants"

        # Check tool version
        echo "\$(timestamp) Checking tool version..."
        TOOL_VERSION_OUTPUT=\$(${params.annotation_tool_cmd} ${params.tool_version_command} 2>&1)
        TOOL_VERSION_EXIT_CODE=\$?
        if [ \$TOOL_VERSION_EXIT_CODE -ne 0 ]; then
            echo "\$(timestamp) ERROR: Failed to determine tool version. Exit code: \$TOOL_VERSION_EXIT_CODE"
            echo "\$(timestamp) Command output: \$TOOL_VERSION_OUTPUT"
            exit 1
        fi

        TOOL_VERSION=\$(echo "\$TOOL_VERSION_OUTPUT" | ${params.tool_version_regex})
        if [ -z "\$TOOL_VERSION" ]; then
            echo "\$(timestamp) ERROR: Unable to extract tool version using regex. Output: \$TOOL_VERSION_OUTPUT"
            exit 1
        fi
        echo "\$(timestamp) Tool version detected: \$TOOL_VERSION"

        if [ "\$TOOL_VERSION" != "${params.required_tool_version}" ]; then
            echo "\$(timestamp) ERROR: Tool version mismatch. Found \$TOOL_VERSION but required ${params.required_tool_version}"
            exit 1
        fi

        # Begin annotation
        echo "\$(timestamp) =========================================="
        echo "\$(timestamp) STARTING ANNOTATION COMMAND"
        echo "\$(timestamp) Command: ${params.annotation_cmd}"
        echo "\$(timestamp) =========================================="

        # Wrap the annotation command in time to measure performance
        # Use set -o pipefail to ensure pipeline errors are caught
        set -o pipefail
        START_TIME=\$(date +%s)

        # Execute the annotation command and capture its exit code
        {
            echo "\$(timestamp) ANNOTATION COMMAND OUTPUT BEGIN:"
            echo ""
            ${params.annotation_cmd}
            ANNOTATION_CMD_STATUS=\$?
            echo ""
            echo "\$(timestamp) ANNOTATION COMMAND OUTPUT END"
        } 2>&1

        END_TIME=\$(date +%s)
        DURATION=\$((END_TIME - START_TIME))

        echo "\$(timestamp) Annotation command completed with exit code: \$ANNOTATION_CMD_STATUS (took \$DURATION seconds)"

        if [ \$ANNOTATION_CMD_STATUS -ne 0 ]; then
            echo "\$(timestamp) ERROR: Annotation command failed with exit code \$ANNOTATION_CMD_STATUS"
            exit \$ANNOTATION_CMD_STATUS
        fi

        # Validate output file was created
        if [ ! -f "\$OUTPUT_BCF" ]; then
            echo "\$(timestamp) ERROR: Output file \$OUTPUT_BCF was not created"
            exit 1
        fi

        # Create index if it doesn't exist
        if [ ! -f "\$OUTPUT_BCF.csi" ]; then
            echo "\$(timestamp) Creating index for output BCF..."
            bcftools index -f "\$OUTPUT_BCF"
            if [ \$? -ne 0 ]; then
                echo "\$(timestamp) ERROR: Failed to create index for output BCF"
                exit 1
            fi
        fi

        # Check output variants - redirect stderr to /dev/null to suppress warnings
        output_variants=\$(bcftools index -n "\$OUTPUT_BCF".csi 2>/dev/null || echo "ERROR_COUNTING")
        if [[ "\$output_variants" == "ERROR_COUNTING" ]]; then
            echo "\$(timestamp) WARNING: Could not count output variants from index. Will try using stats..."
            output_variants=\$(bcftools stats "\$OUTPUT_BCF" | grep "number of records:" | awk '{print \$4}' || echo "UNKNOWN")
        fi

        echo "\$(timestamp) ANNOTATION SUMMARY:"
        echo "\$(timestamp) - Input variants: \$input_variants"
        echo "\$(timestamp) - Output variants: \$output_variants"
        echo "\$(timestamp) - Duration: \$DURATION seconds"

        # Check if we have the required INFO tag
        if [[ -n "${params.must_contain_info_tag}" ]]; then
            echo "\$(timestamp) Checking for required INFO tag: ${params.must_contain_info_tag}"
            INFO_TAG_CHECK=\$(bcftools view -h "\$OUTPUT_BCF" | grep "ID=${params.must_contain_info_tag}," || echo "")

            if [ -z "\$INFO_TAG_CHECK" ]; then
                echo "\$(timestamp) ERROR: Required INFO tag ${params.must_contain_info_tag} not found in output BCF"
                echo "\$(timestamp) Header contains these INFO tags:"
                bcftools view -h "\$OUTPUT_BCF" | grep "^##INFO"
                exit 1
            else
                echo "\$(timestamp) Required INFO tag ${params.must_contain_info_tag} found in output BCF"
            fi
        fi

        # Check if output count is too low compared to input (potential data loss)
        if [[ "\$output_variants" != "UNKNOWN" && "\$input_variants" != "UNKNOWN" && "\$input_variants" != "ERROR_COUNTING" ]]; then
            if (( output_variants < input_variants * 90 / 100 )); then
                echo "\$(timestamp) WARNING: Output variant count is significantly lower than input count."
                echo "\$(timestamp) This may indicate data loss during annotation."
                percent=\$((output_variants * 100 / input_variants))
                echo "\$(timestamp) Input: \$input_variants, Output: \$output_variants (\$percent% preserved)"
            fi
        fi
        
        echo "\$(timestamp) =========================================="
        echo "\$(timestamp) ANNOTATION PROCESS COMPLETED SUCCESSFULLY"
        echo "\$(timestamp) =========================================="
        
    } 2>&1 | tee vcfstash_annotated.log
    
    # Modified error check - only look for critical errors, not warnings
    # Check for ERROR: but ignore lines with "ERROR_COUNTING" which is just a marker for our counting function
    grep "ERROR:" vcfstash_annotated.log | grep -v "ERROR_COUNTING" > error_lines.txt || true
    if [ -s error_lines.txt ]; then
        echo "Critical errors found in log file. Failing process."
        cat error_lines.txt
        exit 1
    fi
    """
}

workflow ANNOTATE {
    take:
    input_bcf
    input_bcf_index

    main:
    annotation_result = Annotate(
        input_bcf,
        input_bcf_index
    )

    emit:
    annotated_bcf = annotation_result.anno_bcf
    annotated_bcf_index = annotation_result.anno_bcf_index
    annotated_bcf_log = annotation_result.anno_bcf_log
}