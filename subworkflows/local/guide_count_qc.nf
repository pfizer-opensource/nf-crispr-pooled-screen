//
// Subworkflow to produce QC analysis of guide count tables
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CONCAT_CORRELATION_TABLES } from '../../modules/local/concat_correlation_tables'
include { GUIDE_COUNT_QC_MODULE } from '../../modules/local/guide_count_qc_module'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO GENERATE GUIDE COUNT TABLES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GUIDE_COUNT_QC {

    take:
    count_tables   // Channel: Channel of guide count tables
    prefix         // String: Prefix for output filenames
    control_guides // Path: Path to control guides file (optional)
    annotate_dict  // Path: Path to annotation dictionary file (optional)

    main:

    ch_versions = Channel.empty()

    GUIDE_COUNT_QC_MODULE (count_tables, prefix, control_guides, annotate_dict)
    ch_versions = ch_versions.mix(GUIDE_COUNT_QC_MODULE.out.versions)

    // Concatenate the replicate correlations into a single file
    CONCAT_CORRELATION_TABLES (
        GUIDE_COUNT_QC_MODULE.out.correlationtables.collect(),
        prefix
    )

    emit:
    versions = ch_versions
}