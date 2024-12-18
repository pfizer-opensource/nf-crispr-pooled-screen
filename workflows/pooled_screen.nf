/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist

// Check mandatory parameters
if (params.fastq_dir) {
    try {
        ch_input = file(params.fastq_dir, checkIfExists: true)
    } catch (java.nio.file.NoSuchFileException e) {
        // return the error message if the path does not exist
        exit 1, "The folder " + e.getMessage() + " does not exist"
    } catch (java.nio.file.NotDirectoryException e) {
        // return the error message if the path is not a directory
        exit 1, "The path " + e.getMessage() + " is not a directory!"
    } catch (java.nio.file.AccessDeniedException e) {
        // return the error message if the path cannot be accessed
        exit 1, "The folder " + e.getMessage() + " cannot be accessed!"
    } catch (e) {
        // return everything else
        exit 1, "The fastq_dir folder definition failed with following reason " + e
    }
} else {
    exit 1, 'Input fastq_dir not specified!'
}

if (params.metadata) {
    try {
        ch_metadata = file(params.metadata, checkIfExists: true)
    } catch (java.nio.file.NoSuchFileException e) {
        // return the error message if the path does not exist
        exit 1, "The file " + e.getMessage() + " does not exist or cannot be accessed!"
    } catch (java.nio.file.AccessDeniedException e) {
        // return the error message if the path cannot be accessed
        exit 1, "The file " + e.getMessage() + " cannot be accessed!"
    } catch (e) {
        // return everything else
        exit 1, "The metadata file definition failed with following reason " + e
    }
} else {
    exit 1, 'Metadata parameter not specified!'
}
if (params.metadata_format) {
    ch_table_format = params.metadata_format
} else {
    exit 1, 'Metadata format not specified!'
}

if (params.library) {
    try {
        ch_library = file(params.library, checkIfExists: true)
    } catch (java.nio.file.NoSuchFileException e) {
        // return the error message if the path does not exist
        exit 1, "The file " + e.getMessage() + " does not exist or cannot be accessed!"
    } catch (java.nio.file.AccessDeniedException e) {
        // return the error message if the path cannot be accessed
        exit 1, "The file " + e.getMessage() + " cannot be accessed!"
    } catch (e) {
        // return everything else
        exit 1, "The library file definition failed with following reason " + e
    }
} else {
    exit 1, 'Library parameter not specified!'
}

if (params.prefix) {
    ch_prefix = params.prefix
} else {
    exit 1, 'Prefix parameter not specified!'
}

if (params.control_guides) {
    try {
        ch_control_guides = file(params.control_guides, checkIfExists: true)
    } catch (java.nio.file.NoSuchFileException e) {
        // return the error message if the path does not exist
        exit 1, "The file " + e.getMessage() + " does not exist or cannot be accessed!"
    } catch (java.nio.file.AccessDeniedException e) {
        // return the error message if the path cannot be accessed
        exit 1, "The file " + e.getMessage() + " cannot be accessed!"
    } catch (e) {
        // return everything else
        exit 1, "The file definition failed with following reason " + e
    }
} else {
    // Since the control_guides file is optional, if it is not provided, then need to
    // define the channel as [] to avoid errors when it is passed as a module input
    ch_control_guides = []
}

// Since annotate_dictionary is optional, one needs to define as [] to avoid
// the absence of the channel ch_annotate_dict
if (params.annotate_dictionary) {
    try {
        ch_annotate_dict = file(params.annotate_dictionary, checkIfExists: true)
    } catch (java.nio.file.NoSuchFileException e) {
        // return the error message if the path does not exist
        exit 1, "The file " + e.getMessage() + " does not exist or cannot be accessed!"
    } catch (java.nio.file.AccessDeniedException e) {
        // return the error message if the path cannot be accessed
        exit 1, "The file " + e.getMessage() + " cannot be accessed!"
    } catch (e) {
        // return everything else
        exit 1, "The annotate dictionary file definition failed with following reason " + e
    }
} else {
    ch_annotate_dict = []
}

if (params.run_mageck_mle) {
    if (params.design_matrix) {
        try {
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
                        [ analysis_name: '', file: file(it, checkIfExists: true) ]
                    } else {
                        [ analysis_name: parts[0], file: file(parts[1], checkIfExists: true) ]
                    }
                }
        } catch (java.nio.file.NoSuchFileException e) {
            exit 1, "The file " + e.getMessage() + " does not exist or cannot be accessed!"
        } catch (java.nio.file.AccessDeniedException e) {
            exit 1, "The file " + e.getMessage() + " cannot be accessed!"
        } catch (e) {
            exit 1, "failed with following reason " + e
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

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { GUIDE_COUNTS_REP } from '../subworkflows/local/guide_counts_rep_subworkflow'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PARSE_METADATA_CSV } from '../modules/local/parse_metadata_csv.nf'
include { PARSE_METADATA_ILLUMINA } from '../modules/local/parse_metadata_illumina.nf'
include { VALIDATE_METADATA } from '../modules/local/validate_metadata.nf'
include { MAGECK_REPORT_HTML } from '../modules/local/guide_generate_mageck_report'
include { GUIDE_CORR_REP_TOTAL } from '../modules/local/guide_correlation_total_table.nf'
include { GUIDE_COUNT_QC_MODULE } from '../modules/local/guide_count_qc_module.nf'

//
// MODULE: Installed directly from nf-core/modules
//
include { MAGECK_COUNT 	      	      } from '../modules/pfizer/mageck_pfizer/count/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MAGECK_MLE } from '../modules/nf-core/mageck/mle/main'
include { MAGECK_TEST                 } from '../modules/nf-core/mageck/test/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow POOLED_SCREEN {

    ch_versions = Channel.empty()
    def PARSED_METADATA_YAML
    def PARSED_METADATA_VERSIONS
    if (ch_table_format == 'illumina') {
        PARSE_METADATA_ILLUMINA(ch_metadata, ch_prefix)
        PARSED_METADATA_YAML = PARSE_METADATA_ILLUMINA.out.meta_yaml
        PARSED_METADATA_VERSIONS = PARSE_METADATA_ILLUMINA.out.versions
    } else {
        // if it is not illumina, it is assumed to be csv cause it would have broken earlier if nothing is specified
        PARSE_METADATA_CSV(ch_metadata, ch_prefix)
        PARSED_METADATA_YAML = PARSE_METADATA_CSV.out.meta_yaml
        PARSED_METADATA_VERSIONS = PARSE_METADATA_CSV.out.versions
    }

    def VALIDATED_METADATA_YAML
    def VALIDATED_METADATA_VERSIONS
    VALIDATE_METADATA(PARSED_METADATA_YAML, ch_input, ch_prefix)
    VALIDATED_METADATA_YAML = VALIDATE_METADATA.out.meta_yaml
    VALIDATED_METADATA_VERSIONS = VALIDATE_METADATA.out.versions

    // check if validated metadata YAML has multiple read entries
    VALIDATED_METADATA_YAML
        .map { PooledUtils.multipleRead(it) }
        .set { multi_read }

    if (multi_read == true && params.count_program == 'mageck') {
        exit 1, 'MAGECK does not support multiple read files per sample. Please use a different count program.'
    }

    if (params.count_program == 'mageck') {
        // store labels and fastqfiles in channels from PooledUtils.get_labels_fastq_list
        VALIDATED_METADATA_YAML
            .map { PooledUtils.getSampleLabelFastqMapForMageck(it) }
            .multiMap{ it ->
                labels: it[0]
                fastq_list: it[1]
                grp: it[2]
                is_ref: it[3]
            }
            .set { result }
        labels = result.labels
        fastq_list = result.fastq_list

        ch_versions = ch_versions.mix(PARSED_METADATA_VERSIONS, VALIDATED_METADATA_VERSIONS)
        // TODO: make sure the multiple read_* is comma separated when building the fastq file input arg
        MAGECK_COUNT (
            ch_input,
            labels,
            fastq_list,
            ch_library,
            ch_prefix
        )
        ch_versions = ch_versions.mix(MAGECK_COUNT.out.versions)

        MAGECK_REPORT_HTML (
            MAGECK_COUNT.out.report,
            ch_prefix
        )

        ch_versions = ch_versions.mix(MAGECK_REPORT_HTML.out.versions)

        GUIDE_COUNTS_REP (
            MAGECK_COUNT.out.count,
            VALIDATED_METADATA_YAML,
            ch_prefix
        )
        ch_versions = ch_versions.mix(GUIDE_COUNTS_REP.out.versions)
    }
    // Ensure the same order of labels, file names and file paths
    // while creating a list of channels to execute the analysis in parallel
    count_tables = GUIDE_COUNTS_REP.out.label_list
                .flatten()
                .merge(GUIDE_COUNTS_REP.out.label_files.flatten())
                .merge(GUIDE_COUNTS_REP.out.rep_count_files.flatten())
                .merge(GUIDE_COUNTS_REP.out.rep_count_files_yml.flatten())

    if (ch_design_matrix){
        // Note: The non-normalized count files are used as input for MLE, and the
        // normalization will be done within MLE based on the value of the normalization_method
        // parameter (modules.config is used to pass the parameter to the MAGECK_MLE process)
        mle_inputs = count_tables
                .combine(ch_design_matrix)
                .multiMap { label, _, count_file, count_yaml, design_matrix ->
                    def prefix = ch_prefix
                    prefix += (design_matrix.analysis_name ? ".${design_matrix.analysis_name}" : '')
                    prefix += (label ? ".${label}X" : '')
                    sample: [[id: prefix, representation: label], count_file]
                    design_matrix: design_matrix.file
                }
        MAGECK_MLE(mle_inputs.sample, mle_inputs.design_matrix)
    }

    GUIDE_COUNT_QC_MODULE (count_tables, ch_prefix, ch_control_guides, ch_annotate_dict)
    ch_versions = ch_versions.mix(GUIDE_COUNT_QC_MODULE.out.versions)

    // merging the *correlation_representation* into one file using collectFile
    GUIDE_CORR_REP_TOTAL(GUIDE_COUNT_QC_MODULE.out.correlationtables.toList(), ch_prefix)
    ch_versions = ch_versions.mix(GUIDE_CORR_REP_TOTAL.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    if (params.run_mageck_test) {
        def mageck_test_inputs = []
        // Note: The non-normalized count files are used as input for MAGeCK test, and the
        // normalization will be done within MAGeCK test based on the value of the normalization_method
        // parameter (modules.config is used to pass the parameter to the MAGECK_TEST process)
        // TODO: Checks whether there are contrasts to run should be moved to metadata processing
        if (ch_table_format == 'illumina') {
            mageck_test_inputs = count_tables.map{label, _, count_file, count_yaml ->
                argMap = PooledUtils.createMageckTestArgumentListIllumina(count_yaml)
                if (! argMap) {
                    exit 1, 'At least one contrast required to run MAGeCK test!'
                }
                metaList = [];
                for(args : argMap){
                    metaList.add([
                        shouldRun: args.run_mageck_test,
                        contrast: args.contrast_prefix,
                        reference : args.c.join(","),
                        treatment : args.t.join(","),
                        representation: label,
                        prefix : "${ch_prefix}.${args.group}.vs.${args.controlGroup}.${label}X",
                        count_table : count_file
                    ])
                }
                metaList
            }
        } else {
            // if it is not illumina, it is assumed to be csv cause it would have broken earlier if nothing is specified
            mageck_test_inputs = count_tables.map{label, _, count_file, count_yaml ->
                argMap = PooledUtils.createMageckTestArgumentMap(count_yaml)
                if (! argMap.run_mageck_test) {
                    exit 1, 'At least one contrast required to run MAGeCK test!'
                }
                metaList = [];
                i=0;
                for(tArg : argMap.t){
                    metaList.add([
                        shouldRun: argMap.run_mageck_test,
                        contrast: '',
                        reference : argMap.c,
                        treatment : tArg.join(","),
                        representation: label,
                        prefix : "${ch_prefix}.${argMap.group[i++]}.vs.${argMap.controlGroup}.${label}X",
                        count_table : count_file
                    ])
                }
                metaList
            }
        }
        MAGECK_TEST(mageck_test_inputs.flatten())
    }
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
