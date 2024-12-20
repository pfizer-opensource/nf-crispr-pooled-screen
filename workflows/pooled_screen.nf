/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate parameters against nextflow_schema.json, including requiring mandatary params
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate required path parameters and convert to Path objects
Path ch_fastq_dir = PooledUtils.validatePathParam(params.fastq_dir, 'fastq_dir', true)
Path ch_metadata = PooledUtils.validatePathParam(params.metadata, 'metadata')
Path ch_library = PooledUtils.validatePathParam(params.library, 'library')

// Validate optional path parameters and convert to Path objects
// Note: if the parameters are not provided, then the channel needs to be defined as []
// to avoid errors when it is passed as a module input
if (params.control_guides) {
    ch_control_guides = PooledUtils.validatePathParam(params.control_guides, 'control_guides')
} else {
    ch_control_guides = []
}
if (params.annotate_dictionary) {
    ch_annotate_dict = PooledUtils.validatePathParam(params.annotate_dictionary, 'annotate_dictionary')
} else {
    ch_annotate_dict = []
}

if (params.run_mageck_mle) {
    if (params.design_matrix) {
        // params.design_matrix is a comma-separated list of paths, where each path can
        // be prefixed with an optional analysis name associated with the design matrix.
        // A colon is used as delimiter, but since that is also part of S3 (and other)
        // URLs, no splitting is done if '://' is present. This overall approach could
        // be improved on in the future to avoid edge cases with colons in paths
        ch_design_matrix = Channel
            .fromList( params.design_matrix.split(',') as List )
            .map {
                def parts = it.split(':', 2)
                if (it.contains('://') || parts.size() == 1) {
                    [
                        analysis_name: '',
                        file: PooledUtils.validatePathParam(it, 'design_matrix')
                    ]
                } else {
                    [
                        analysis_name: parts[0],
                        file: PooledUtils.validatePathParam(parts[1], 'design_matrix')
                    ]
                }
            }
    } else {
        exit 1, 'Design matrix required to run MAGeCK MLE!'
    }
} else {
    ch_design_matrix = null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GUIDE_COUNTS } from '../subworkflows/local/guide_counts.nf'
include { GUIDE_COUNT_QC } from '../subworkflows/local/guide_count_qc.nf'
include { LOAD_METADATA } from '../subworkflows/local/metadata.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MAGECK_MLE } from '../modules/nf-core/mageck/mle/main'
include { MAGECK_TEST } from '../modules/nf-core/mageck/test/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow POOLED_SCREEN {

    ch_versions = Channel.empty()

    LOAD_METADATA (ch_metadata, params.metadata_format, ch_fastq_dir, params.prefix)
    ch_metadata_yaml = LOAD_METADATA.out.meta_yaml
    ch_versions = ch_versions.mix(LOAD_METADATA.out.versions)

    GUIDE_COUNTS (
        ch_metadata_yaml,
        ch_fastq_dir,
        ch_library,
        ch_control_guides,
        params.prefix
    )
    count_tables = GUIDE_COUNTS.out.count_tables
    ch_versions = ch_versions.mix(GUIDE_COUNTS.out.versions)

    GUIDE_COUNT_QC (count_tables, params.prefix, ch_control_guides, ch_annotate_dict)
    ch_versions = ch_versions.mix(GUIDE_COUNT_QC.out.versions)

    if (ch_design_matrix) {
        // Note: The non-normalized count files are used as input for MLE, and the
        // normalization will be done within MLE based on the value of the normalization_method
        // parameter (modules.config is used to pass the parameter to the MAGECK_MLE process)
        mle_inputs = count_tables
                .combine(ch_design_matrix)
                .multiMap { representation, count_file, count_yaml, design_matrix ->
                    def prefix = params.prefix
                    prefix += (design_matrix.analysis_name ? ".${design_matrix.analysis_name}" : '')
                    prefix += (representation ? ".${representation}X" : '')
                    sample: [[id: prefix, representation: representation], count_file]
                    design_matrix: design_matrix.file
                }
        MAGECK_MLE (
            mle_inputs.sample,
            mle_inputs.design_matrix,
            ch_control_guides
        )
        ch_versions = ch_versions.mix(MAGECK_MLE.out.versions)
    }

    if (params.run_mageck_test) {
        // Note: The non-normalized count files are used as input for MAGeCK test, and the
        // normalization will be done within MAGeCK test based on the value of the normalization_method
        // parameter (modules.config is used to pass the parameter to the MAGECK_TEST process)
        // TODO: Checks whether there are contrasts to run should be moved to metadata processing
        mageck_test_inputs = count_tables.map { representation, count_file, count_yaml ->
            contrasts = PooledUtils.createMageckTestContrasts(count_yaml)
            if (! contrasts) {
                error 'At least one contrast required to run MAGeCK test!'
            }
            contrasts.collect {
                [
                    meta: [
                        id: "${params.prefix}.${it.group}.vs.${it.refGroup}.${representation}X",
                        contrast: it.contrast,
                        reference: it.refSamples.join(","),
                        treatment: it.samples.join(","),
                        representation: representation
                    ],
                    count_table: count_file
                ]
            }
        }
        MAGECK_TEST(mageck_test_inputs.flatten(), ch_control_guides)
        ch_versions = ch_versions.mix(MAGECK_TEST.out.versions)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow.onComplete {
    if (params.email) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }

    NfcoreTemplate.summary(workflow, params, log)
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
