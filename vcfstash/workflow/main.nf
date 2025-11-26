nextflow.enable.dsl = 2

include { NORMALIZE; FILTER_ONLY } from './modules/normalize'
include { ANNOTATE } from './modules/annotate'
include { AnnotateFromDB; FilterMissingAnnotations; MergeAnnotations } from './modules/annotate_cache'
include { MERGE } from './modules/merge'
include { MERGE_VARIANTS } from './modules/merge_variants'
include { UTILS } from './modules/utils'

// Then create a function to copy auxiliary files
def copyAuxiliaryFiles(auxFiles, targetDir) {
    auxFiles.flatten().ifEmpty([]).subscribe { aux_file ->
        if (aux_file && file(aux_file).exists()) {
            def auxiliaryDir = file(targetDir)
            if (!auxiliaryDir.exists()) {
                auxiliaryDir.mkdirs()
            }
            // Copy the file directly to the auxiliary directory without creating subdirectories
            file(aux_file).copyTo("${auxiliaryDir}/${file(aux_file).name}")
        }
    }
}

def validateParams() {

    def validModes = ['annotate', 'stash-add', 'stash-init', 'stash-annotate', 'annotate-nocache', 'annotate-isec']
    if (!validModes.contains(params.db_mode)) {
        error "Invalid db_mode value: ${params.db_mode}. Must be one of: ${validModes.join(', ')}"
    }

    // Required params
    def required = ['output']

    // Input is required except in special cases
    if (!params.input && !params.db_bcf) {
        error "Either input file (--input) or database BCF file (--db_bcf) is required"
    }

    for (param in required) {
        if (params[param] == null) {
            error "Required parameter missing: --${param}"
        }
    }

    // Validate input file exists if specified
    if (params.input && !file(params.input).exists()) {
        error "Input file not found: ${params.input}"
    }

    // Validate db_bcf file exists if specified
    if (params.db_bcf != null && !file(params.db_bcf).exists()) {
       error "Database BCF file not found: ${params.db_bcf}"
    }

    if (params.db_mode == 'annotate' && !params.db_bcf) {
        error "Database BCF file (--db_bcf) is required for annotate mode"
    }

    // Require must_contain_info_tag for annotation-related modes
    if ((params.db_mode == 'stash-annotate' || params.db_mode == 'annotate') &&
        (!params.containsKey('must_contain_info_tag') || params.must_contain_info_tag.trim().isEmpty())) {
        error "must_contain_info_tag parameter is required and cannot be empty for annotation modes"
    }

    // Validate output directory
    def outDir = file(params.output)
    if (!outDir.exists()) {
        outDir.mkdirs()
    }
}

workflow {

    params.db_mode = null
    params.db_bcf = null
	params.auxiliary_dir = "${params.output}/auxiliary"

    validateParams()

    // Extract sample name from the input file
    def sampleName
    if (params.input) {
        sampleName = params.input.tokenize('/')[-1].replace('.vcf', '').replace('.bcf', '')
    } else if (params.db_bcf) {
        sampleName = params.db_bcf.tokenize('/')[-1].replace('.vcf', '').replace('.bcf', '')
    }

    def outputDir = file(params.output)

    // Only require chr_add when normalization is enabled
    def chr_add = null

    // Check if normalization is enabled (default to false if not specified)
    def normalize = params.containsKey('normalize') ? params.normalize.toString().toBoolean() : false

    if (normalize && (params.db_mode == 'stash-init' || params.db_mode == 'stash-add')) {
        // Validate that chr_add is provided when normalization is enabled
        if (!params.containsKey('chr_add') || !params.chr_add) {
            error "chr_add parameter is required when normalization is enabled"
        }

        chr_add = file(params.chr_add)
    }

    if (params.db_mode == 'stash-annotate') {
        // For annotation mode, use db_bcf as the input
        vcf = file(params.db_bcf)

        // Handle index file - need to make sure it's a simple file path
        vcf_index = file("${params.db_bcf}.csi")

        // Direct annotation without normalization
        ANNOTATE(
            vcf,
            vcf_index
        )

        // Publish annotated database
        ANNOTATE.out.annotated_bcf.subscribe { bcf ->
            file(bcf).copyTo("${outputDir}/vcfstash_annotated.bcf")
        }
        ANNOTATE.out.annotated_bcf_index.subscribe { idx ->
            file(idx).copyTo("${outputDir}/vcfstash_annotated.bcf.csi")
        }
        ANNOTATE.out.annotated_bcf_log.subscribe { log ->
            file(log).copyTo("${outputDir}/vcfstash_annotated.bcf.log")
        }

        // Copy all auxiliary files using flatten to handle the collection properly
		copyAuxiliaryFiles(ANNOTATE.out.aux_files, params.auxiliary_dir)

    } else { // in all other modes we're gonna have some sort of input
        vcf = file(params.input)
        vcf_index = file("${params.input}.csi")

        // Only run UTILS process when normalization is enabled
        if (normalize && (params.db_mode == 'stash-init' || params.db_mode == 'stash-add')) {
            UTILS(
                sampleName,
                outputDir,
                vcf,
                chr_add
            )
        }

        // Only apply normalization for stash-init and stash-add, and only if normalize is true
        // Note: normalize is defined above, we're using it here
        def remove_gt = params.db_mode.startsWith('stash-')
        def remove_info = params.db_mode.startsWith('stash-')

        // Define normalized_vcf and normalized_vcf_index variables
        def normalized_vcf
        def normalized_vcf_index
        def normalized_vcf_log

        if (params.db_mode == 'stash-init' || params.db_mode == 'stash-add') {
            if (normalize) {
                // Apply full normalization for stash-init and stash-add if normalize is true
                NORMALIZE(
                    chr_add,
                    vcf,
                    vcf_index,  // Pass the index file explicitly
                    sampleName,
                    remove_gt,
                    remove_info
                )
                normalized_vcf = NORMALIZE.out.norm_bcf
                normalized_vcf_index = NORMALIZE.out.norm_bcf_index
                normalized_vcf_log = NORMALIZE.out.norm_bcf_log
            } else {
                // Apply only GT and INFO removal if normalize is false
                FILTER_ONLY(
                    vcf,
                    vcf_index,
                    sampleName,
                    remove_gt,
                    remove_info
                )
                normalized_vcf = FILTER_ONLY.out.filtered_bcf
                normalized_vcf_index = FILTER_ONLY.out.filtered_bcf_index
                normalized_vcf_log = FILTER_ONLY.out.filtered_bcf_log
            }
        } else {
            // For other modes, use the input VCF directly
            normalized_vcf = vcf
            normalized_vcf_index = vcf_index
            normalized_vcf_log = Channel.empty()
        }

        // stash-add: DATABASE_UPDATE_WORKFLOW - Add new variants to database
        if (params.db_mode == 'stash-add') {
            // If database exists, merge with it (stash-add mode)
            db_bcf = file(params.db_bcf)
            db_bcf_index = file("${params.db_bcf}.csi")
            MERGE_VARIANTS(
                db_bcf,
                db_bcf_index,
                normalized_vcf,
                normalized_vcf_index
            )

            // Publish final merged BCF without annotation
            MERGE_VARIANTS.out.merged_bcf.subscribe { bcf ->
                file(bcf).copyTo("${outputDir}/vcfstash.bcf")
            }
            MERGE_VARIANTS.out.merged_bcf_index.subscribe { idx ->
                file(idx).copyTo("${outputDir}/vcfstash.bcf.csi")
            }

            // Only copy log if normalization was performed
            if ((params.db_mode == 'stash-init' || params.db_mode == 'stash-add') && normalize) {
                normalized_vcf_log.subscribe { log ->
                    file(log).copyTo("${outputDir}/add_${sampleName}.bcf.log")
                }
            }
        } else if (params.db_mode == 'stash-init') {
            // stash-init: First database creation - just use the input (normalized or not)
            // No annotation here
            normalized_vcf.subscribe { bcf ->
                file(bcf).copyTo("${outputDir}/vcfstash.bcf")
            }
            normalized_vcf_index.subscribe { idx ->
                file(idx).copyTo("${outputDir}/vcfstash.bcf.csi")
            }

            // Only copy log if normalization was performed
            if ((params.db_mode == 'stash-init' || params.db_mode == 'stash-add') && normalize) {
                normalized_vcf_log.subscribe { log ->
                    file(log).copyTo("${outputDir}/vcfstash.bcf.log")
                }
            }

        } else if (params.db_mode == 'annotate-nocache') {
            // annotate: DIRECT_ANNOTATION_WORKFLOW - Direct VCF annotation without database comparison
            ANNOTATE(
                normalized_vcf,
                normalized_vcf_index
            )

            // Publish annotated database
            ANNOTATE.out.annotated_bcf.subscribe { bcf ->
                file(bcf).copyTo("${outputDir}/${sampleName}_vst.bcf")
            }
            ANNOTATE.out.annotated_bcf_index.subscribe { idx ->
                file(idx).copyTo("${outputDir}/${sampleName}_vst.bcf.csi")
            }
            ANNOTATE.out.annotated_bcf_log.subscribe { log ->
                file(log).copyTo("${outputDir}/${sampleName}_vst.bcf.log")
            }

            copyAuxiliaryFiles(ANNOTATE.out.aux_files, params.auxiliary_dir)

        } else if (params.db_mode == 'annotate') {
            // annotate: SAMPLE_ANALYSIS_WORKFLOW - Sample comparison against database using bcftools annotate
            db_bcf = file(params.db_bcf)
            db_bcf_index = file("${params.db_bcf}.csi")

            // Step 1: Annotate input with database INFO fields
            AnnotateFromDB(
                normalized_vcf,
                normalized_vcf_index,
                db_bcf,
                db_bcf_index,
                sampleName
            )

            // Step 2: Filter to get variants with missing CSQ annotations
            FilterMissingAnnotations(
                AnnotateFromDB.out.annotated_bcf,
                AnnotateFromDB.out.annotated_bcf_index,
                sampleName
            )

            // Step 3: Run standard ANNOTATE workflow on missing annotations
            ANNOTATE(
                FilterMissingAnnotations.out.filtered_bcf,
                FilterMissingAnnotations.out.filtered_bcf_index
            )

            // Step 4: Merge newly annotated variants back into original file
            MergeAnnotations(
                AnnotateFromDB.out.annotated_bcf,
                AnnotateFromDB.out.annotated_bcf_index,
                ANNOTATE.out.annotated_bcf,
                ANNOTATE.out.annotated_bcf_index,
                sampleName
            )

            // Publish the results
            MergeAnnotations.out.merged_bcf.subscribe { bcf ->
                file(bcf).copyTo("${outputDir}/${sampleName}_vst.bcf")
            }

            MergeAnnotations.out.merged_bcf_index.subscribe { idx ->
                file(idx).copyTo("${outputDir}/${sampleName}_vst.bcf.csi")
            }
            ANNOTATE.out.annotated_bcf_log.subscribe { log ->
                file(log).copyTo("${outputDir}/${sampleName}_vst.bcf.log")
            }
            // Only copy log if normalization was performed
            if ((params.db_mode == 'stash-init' || params.db_mode == 'stash-add') && normalize) {
                normalized_vcf_log.subscribe { log ->
                    file(log).copyTo("${outputDir}/${sampleName}_norm.bcf.log")
                }
            }

            // Copy all auxiliary files using flatten to handle the collection properly
			copyAuxiliaryFiles(ANNOTATE.out.aux_files, params.auxiliary_dir)
        }
    }
}
