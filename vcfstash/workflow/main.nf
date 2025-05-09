nextflow.enable.dsl = 2

include { NORMALIZE } from './modules/normalize'
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

    chr_add = file(params.chr_add)
    reference = file(params.reference)

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

    } else { // in all other modes we're gonna have some sort of input that requires normalization
        vcf = file(params.input)
        vcf_index = file("${params.input}.csi")

        UTILS(
            sampleName,
            outputDir,
            vcf,
            chr_add,
            reference
        )

        def remove_gt = params.db_mode.startsWith('stash-')
        def remove_info = params.db_mode.startsWith('stash-')

        NORMALIZE(
            chr_add,
            vcf,
            vcf_index,  // Pass the index file explicitly
            reference,
            sampleName,
            false,
            remove_gt,
            remove_info
        )

        // stash-add: DATABASE_UPDATE_WORKFLOW - Add new variants to database
        if (params.db_mode == 'stash-add') {
            // If database exists, merge with it (stash-add mode)
            db_bcf = file(params.db_bcf)
            db_bcf_index = file("${params.db_bcf}.csi")
            MERGE_VARIANTS(
                db_bcf,
                db_bcf_index,
                NORMALIZE.out.norm_bcf,
                NORMALIZE.out.norm_bcf_index
            )

            // Publish final merged BCF without annotation
            MERGE_VARIANTS.out.merged_bcf.subscribe { bcf ->
                file(bcf).copyTo("${outputDir}/vcfstash.bcf")
            }
            MERGE_VARIANTS.out.merged_bcf_index.subscribe { idx ->
                file(idx).copyTo("${outputDir}/vcfstash.bcf.csi")
            }
            NORMALIZE.out.norm_bcf_log.subscribe { log ->
                file(log).copyTo("${outputDir}/add_${sampleName}.bcf.log")
            }
        } else if (params.db_mode == 'stash-init') {
            // stash-init: First database creation - just normalize the input
            // No annotation here
            NORMALIZE.out.norm_bcf.subscribe { bcf ->
                file(bcf).copyTo("${outputDir}/vcfstash.bcf")
            }
            NORMALIZE.out.norm_bcf_index.subscribe { idx ->
                file(idx).copyTo("${outputDir}/vcfstash.bcf.csi")
            }
            NORMALIZE.out.norm_bcf_log.subscribe { log ->
                file(log).copyTo("${outputDir}/vcfstash.bcf.log")
            }

        } else if (params.db_mode == 'annotate-nocache') {
            // annotate: DIRECT_ANNOTATION_WORKFLOW - Direct VCF annotation without database comparison
            ANNOTATE(
                NORMALIZE.out.norm_bcf,
                NORMALIZE.out.norm_bcf_index
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
                NORMALIZE.out.norm_bcf,
                NORMALIZE.out.norm_bcf_index,
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
            NORMALIZE.out.norm_bcf_log.subscribe { log ->
                file(log).copyTo("${outputDir}/${sampleName}_norm.bcf.log")
            }

            // Copy all auxiliary files using flatten to handle the collection properly
			copyAuxiliaryFiles(ANNOTATE.out.aux_files, params.auxiliary_dir)
        }
    }
}