/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist

// Check mandatory parameters
def checkPathParamList = [ params.fastq_dir, params.outdir, params.sample_input_meta_table, params.library ]
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

if (params.sample_input_meta_table) {
    try {
        ch_sample_input_meta_table = file(params.sample_input_meta_table, checkIfExists: true)
    } catch (java.nio.file.NoSuchFileException e) {
        // return the error message if the path does not exist
        exit 1, "The file " + e.getMessage() + " does not exist or cannot be accessed!"
    } catch (java.nio.file.AccessDeniedException e) {
        // return the error message if the path cannot be accessed
        exit 1, "The file " + e.getMessage() + " cannot be accessed!"
    } catch (e) {
        // return everything else
        exit 1, "The sample_input_meta_table file definition failed with following reason " + e
    }
} else {
    exit 1, 'Input sample_input_meta_table not specified!'
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
        exit 1, "The sample_input_meta_table file definition failed with following reason " + e
    }
} else {
    exit 1, 'Input library not specified!'
}
if (params.prefix) { ch_prefix = params.prefix }

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
        exit 1, "The sample_input_meta_table file definition failed with following reason " + e
    }
} else {
    ch_annotate_dict = []
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
include { SAMPLE_META_MODULE } from '../modules/local/sample_meta_module.nf'
include { MAGECK_REPORT_HTML } from '../modules/local/guide_generate_mageck_report'
include { GUIDE_CORR_REP_TOTAL } from '../modules/local/guide_correlation_total_table.nf'
include { GUIDE_COUNT_QC_MODULE } from '../modules/local/guide_count_qc_module.nf'

//
// MODULE: Installed directly from nf-core/modules
//
include { MAGECK_COUNT 	      	      } from '../modules/pfizer/mageck_pfizer/count/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow POOLED_SCREEN {

    ch_versions = Channel.empty()

    SAMPLE_META_MODULE (
        ch_input,
        ch_sample_input_meta_table,
        ch_prefix
    )
    // check if SAMPLE_META.out.sample_meta_yaml has multiple read entries
    SAMPLE_META_MODULE.out.sample_meta_yaml
        .map { PooledUtils.multipleRead(it) }
        .set { multi_read }

    if (multi_read == true && params.count_program == 'mageck') {
        exit 1, 'MAGECK does not support multiple read files per sample. Please use a different count program.'
    }

    if (params.count_program == 'mageck') {
        // store labels and fastqfiles in channels from PooledUtils.get_labels_fastq_list

        SAMPLE_META_MODULE.out.sample_meta_yaml
            .map { PooledUtils.getSampleLabelFastqMapForMageck(it) }
            .multiMap{ it ->
                labels: it[0]
                fastq_list: it[1]
            }
            .set { result }
        labels = result.labels
        fastq_list = result.fastq_list

        ch_versions = ch_versions.mix(SAMPLE_META_MODULE.out.versions)
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
            SAMPLE_META_MODULE.out.sample_meta_yaml,
            ch_prefix
        )
        ch_versions = ch_versions.mix(GUIDE_COUNTS_REP.out.versions)
    }

    // Ensure the same order of labels, file names and file paths
    // while creating a list of channels to execute the analysis in parallel

    my_inputs = GUIDE_COUNTS_REP.out.label_list
                .flatten()
                .merge(GUIDE_COUNTS_REP.out.label_files.flatten()
                .merge(GUIDE_COUNTS_REP.out.rep_count_files.flatten())
                .merge(GUIDE_COUNTS_REP.out.rep_count_files_yml.flatten()))

    GUIDE_COUNT_QC_MODULE (my_inputs, ch_prefix, ch_annotate_dict)
    ch_versions = ch_versions.mix(GUIDE_COUNT_QC_MODULE.out.versions)

    // merging the *correlation_representation* into one file using collectFile
    GUIDE_CORR_REP_TOTAL(GUIDE_COUNT_QC_MODULE.out.correlationtables.toList(), ch_prefix)
    ch_versions = ch_versions.mix(GUIDE_CORR_REP_TOTAL.out.versions)

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
