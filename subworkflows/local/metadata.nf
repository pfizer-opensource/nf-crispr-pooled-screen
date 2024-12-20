//
// Subworkflow to handle metadata processing
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PARSE_METADATA_CSV } from '../../modules/pfizer-rd/parse_metadata/csv/main'
include { PARSE_METADATA_ILLUMINA } from '../../modules/pfizer-rd/parse_metadata/illumina/main'
include { VALIDATE_METADATA } from '../../modules/local/validate_metadata.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO LOAD METADATA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LOAD_METADATA {

    take:
    metadata  // Path: Path to input metadata file
    format    // String: Format of metadata file (ie, illumina, csv, etc)
    fastq_dir // Path: Path to fastq directory
    prefix    // String: Prefix for output filenames

    main:

    switch (format) {
        case "illumina":
            parse_metadata = PARSE_METADATA_ILLUMINA(metadata)
            break
        case "csv":
            parse_metadata = PARSE_METADATA_CSV(metadata)
            break
        default:
            error "Unrecognized metadata format: ${format}"
    }

    VALIDATE_METADATA(parse_metadata.metadata_yaml, fastq_dir, prefix)

    meta_yaml = VALIDATE_METADATA.out.meta_yaml
    versions = parse_metadata.versions.mix(VALIDATE_METADATA.out.versions)


    emit:
    meta_yaml
    versions
}