nextflow.enable.dsl = 2

include { NORMALIZE } from './modules/normalize'
include { INTERSECT } from './modules/intersect'
include { ANNOTATE } from './modules/annotate'
include { MERGE } from './modules/merge'
include { UTILS } from './modules/utils'

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

    // Validate inputs first
    ValidateInputs(vcf, chr_add, reference, vep_cache)
    ValidateInputs.out.validated.ifEmpty { error "Input validation failed" }

    if (params.db_mode) {
        // Database annotation mode - skip normalization and go directly to annotation
        ANNOTATE(vcf, vcf_index, vep_cache, reference, sampleName)
        // Publish final database BCF
        ANNOTATE.out.vep_bcf.subscribe { bcf ->
            file(bcf).copyTo("${outputDir}/${sampleName}.bcf")
        }
        ANNOTATE.out.vep_bcf_index.subscribe { idx ->
            file(idx).copyTo("${outputDir}/${sampleName}.bcf.csi")
        }
    } else if (params.db_bcf == false) {
        // Direct VEP annotation without database
        NORMALIZE(chr_add, vcf, reference, sampleName)
        ANNOTATE(NORMALIZE.out.norm_bcf, NORMALIZE.out.norm_bcf_index, vep_cache, reference, sampleName)
        MERGE(ANNOTATE.out.vep_bcf, ANNOTATE.out.vep_bcf_index, ANNOTATE.out.ev_bcf, ANNOTATE.out.ev_bcf_index, outputDir, sampleName, projectPath)
    } else {
        // Sample annotation mode - full workflow with normalization
        NORMALIZE(chr_add, vcf, reference, sampleName)
        db_bcf = file(params.db_bcf)
        db_bcf_index = file("${params.db_bcf}.csi")
        INTERSECT(NORMALIZE.out.norm_bcf, NORMALIZE.out.norm_bcf_index, db_bcf, db_bcf_index, sampleName)
        ANNOTATE(INTERSECT.out.subset_bcf, INTERSECT.out.subset_bcf_index, vep_cache, reference, sampleName)
        MERGE(INTERSECT.out.annotated_bcf, INTERSECT.out.annotated_bcf_index, ANNOTATE.out.vep_bcf, ANNOTATE.out.vep_bcf_index, outputDir, sampleName, projectPath)
    }

    // Capture tool versions at the end
    CaptureToolVersions(sampleName, outputDir)
}

