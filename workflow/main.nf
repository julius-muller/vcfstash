nextflow.enable.dsl = 2

include { NORMALIZE } from './modules/normalize'
include { INTERSECT } from './modules/intersect'
include { ANNOTATE } from './modules/annotate'
include { ANNOTATE_INFO } from './modules/annotate_info'
include { MERGE } from './modules/merge'
include { MERGE_VARIANTS } from './modules/merge_variants'
include { UTILS } from './modules/utils'

def validateParams() {

    def validModes = ['annotate', 'stash-add', 'stash-init', 'stash-annotate', 'annotate-nocache', 'annotate-isec']
    if (!validModes.contains(params.db_mode)) {
        error "Invalid db_mode value: ${params.db_mode}. Must be one of: ${validModes.join(', ')}"
    }

        // When in annotation mode, set db_bcf to input if not provided
    if (params.db_mode.startsWith('stash-annotate') && params.input) {
        params.db_bcf = params.input
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


    // Validate output directory
    def outDir = file(params.output)
    if (!outDir.exists()) {
        outDir.mkdirs()
    }

}

workflow {
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
    vep_cache = file(params.vep_cache)

    if (params.db_mode == 'stash-annotate') {
        // For annotation mode, use db_bcf as the input
        vcf = file(params.db_bcf)
        vcf_index = file("${params.db_bcf}.csi")

        // Direct annotation without normalization
        ANNOTATE(
			vcf,
			vcf_index,
            vep_cache,
            reference,
            sampleName
		)

        // Publish annotated database
        ANNOTATE.out.vep_bcf.subscribe { bcf ->
            file(bcf).copyTo("${outputDir}/vepstash_annotated.bcf")
        }
        ANNOTATE.out.vep_bcf_index.subscribe { idx ->
            file(idx).copyTo("${outputDir}/vepstash_annotated.bcf.csi")
        }
    } else { // in all other modes were gonna have some sort of input that requires normalization
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
                file(bcf).copyTo("${outputDir}/vepstash.bcf")
            }
            MERGE_VARIANTS.out.merged_bcf_index.subscribe { idx ->
                file(idx).copyTo("${outputDir}/vepstash.bcf.csi")
            }
        } else if (params.db_mode == 'stash-init') {
            // stash-init: First database creation - just normalize the input
            // No annotation here
            NORMALIZE.out.norm_bcf.subscribe { bcf ->
                file(bcf).copyTo("${outputDir}/vepstash.bcf")
            }
            NORMALIZE.out.norm_bcf_index.subscribe { idx ->
                file(idx).copyTo("${outputDir}/vepstash.bcf.csi")
            }
		} else if (params.db_mode == 'annotate-nocache') {
			// annotate: DIRECT_ANNOTATION_WORKFLOW - Direct VEP annotation without database comparison
			ANNOTATE(
				NORMALIZE.out.norm_bcf,
				NORMALIZE.out.norm_bcf_index,
				vep_cache,
				reference,
				sampleName
			)

            // Publish annotated database
            ANNOTATE.out.vep_bcf.subscribe { bcf ->
                file(bcf).copyTo("${outputDir}/${sampleName}_vst.bcf")
            }
            ANNOTATE.out.vep_bcf_index.subscribe { idx ->
                file(idx).copyTo("${outputDir}/${sampleName}_vst.bcf.csi")
            }

		} else if (params.db_mode == 'annotate-isec') {
			// annotate: DEPRECATED isec based mode SAMPLE_ANALYSIS_WORKFLOW - Sample comparison against database
			db_bcf = file(params.db_bcf)
			db_bcf_index = file("${params.db_bcf}.csi")
			INTERSECT(
				NORMALIZE.out.norm_bcf,
				NORMALIZE.out.norm_bcf_index,
				db_bcf,
				db_bcf_index,
				sampleName
			)

			ANNOTATE(
				INTERSECT.out.subset_bcf,
				INTERSECT.out.subset_bcf_index,
				vep_cache,
				reference,
				sampleName
			)

			MERGE_VARIANTS(
				INTERSECT.out.annotated_bcf,
				INTERSECT.out.annotated_bcf_index,
				ANNOTATE.out.vep_bcf,
				ANNOTATE.out.vep_bcf_index
			)

				// Publish final merged BCF
				MERGE_VARIANTS.out.merged_bcf.subscribe { bcf ->
					file(bcf).copyTo("${outputDir}/${sampleName}_vst.bcf")
				}
				MERGE_VARIANTS.out.merged_bcf_index.subscribe { idx ->
					file(idx).copyTo("${outputDir}/${sampleName}_vst.bcf.csi")
				}
		} else if (params.db_mode == 'annotate') {
			// annotate: SAMPLE_ANALYSIS_WORKFLOW - Sample comparison against database using bcftools annotate
            db_bcf = file(params.db_bcf)
            db_bcf_index = file("${params.db_bcf}.csi")

            // Run the new ANNOTATE_INFO workflow
            ANNOTATE_INFO(
                NORMALIZE.out.norm_bcf,
                NORMALIZE.out.norm_bcf_index,
                db_bcf,
                db_bcf_index,
                vep_cache,
                reference,
                sampleName
            )

            // Publish the results
            ANNOTATE_INFO.out.merged_bcf.subscribe { bcf ->
                file(bcf).copyTo("${outputDir}/${sampleName}_vst.bcf")
            }
            ANNOTATE_INFO.out.merged_bcf_index.subscribe { idx ->
                file(idx).copyTo("${outputDir}/${sampleName}_vst.bcf.csi")
            }

		}
	}
}


