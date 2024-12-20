//
// Subworkflow to generate guide count tables
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COUNTS_BY_REPRESENTATION } from '../../modules/local/counts_by_representation'
include { MAGECK_REPORT_HTML } from '../../modules/local/guide_generate_mageck_report'
include { MAGECK_COUNT } from '../../modules/nf-core/mageck/count/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO GENERATE GUIDE COUNT TABLES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GUIDE_COUNTS {

    take:
    metadata_yaml  // File Channel: Path to metadata YAML file
    fastq_dir      // Path: Path to fastq directory
    library        // Path: Path to guide library file
    control_guides // Path: Path to control guides file (optional)
    prefix         // String: Prefix for output filenames

    main:

    ch_versions = Channel.empty()

    // store sample names and fastq files in a map with the following structure:
    // [ meta: [id: <comma-separated string of sample names>],
    //   fastqs_r1: [<list of read1 fastq files>],
    //   fastqs_r2: [<list of read2 fastq files if present, or empty list>],
    //   single_end: [boolean] ]
    metadata_yaml
        .map { PooledUtils.getSamplesFastqMap(it, fastq_dir) }
        .set { samples }

    switch (params.count_program) {
        case "mageck":
            if (! samples.single_end) {
                error 'MAGECK does not support multiple read files per sample. Please use a different count program.'
            }

            MAGECK_COUNT (
                samples.subMap(["meta", "fastqs_r1"]),
                library,
                control_guides
            )
            ch_versions = ch_versions.mix(MAGECK_COUNT.out.versions)

            MAGECK_REPORT_HTML (
                MAGECK_COUNT.out.report,
                prefix
            )
            ch_versions = ch_versions.mix(MAGECK_REPORT_HTML.out.versions)

            count = MAGECK_COUNT.out.count
            break

        default:
            error "Unrecognized count program: ${params.count_program}"
    }

    COUNTS_BY_REPRESENTATION (count, metadata_yaml, prefix)
    ch_versions = ch_versions.mix(COUNTS_BY_REPRESENTATION.out.versions)

    COUNTS_BY_REPRESENTATION.out.counts.flatten()
    .map { [it.getBaseName(), it]}
    .set { counts_by_repr }

    COUNTS_BY_REPRESENTATION.out.meta_yamls.flatten()
    .map { [it.getBaseName(), it] }
    .set { meta_by_repr }

    COUNTS_BY_REPRESENTATION.out.repr_name_map
    .map { PooledUtils.getRepresentationBasenames(it) }
    .flatten()
    .collate(2)
    .set { repr_basenames }

    // Create tuples for each representation value and the corresponding file paths
    // to execute the analyses for each representation in parallel
    count_tables = repr_basenames
        .join(counts_by_repr)
        .join(meta_by_repr)
        // Remove basename variable, which was just a join key, from the tuple
        .map { _basename, representation, counts, metadata ->
            [representation, counts, metadata]
        }

    emit:
    count_tables
    versions = ch_versions
}