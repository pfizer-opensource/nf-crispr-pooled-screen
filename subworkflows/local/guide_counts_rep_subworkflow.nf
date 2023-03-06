//
// Generate sample metadata yaml file
//
include { GUIDE_COUNT_REP_MODULE } from '../../modules/local/guide_count_rep_module.nf'


workflow GUIDE_COUNTS_REP {

    take:
    count_file
    sample_meta_yaml
    prefix

    main:

    GUIDE_COUNT_REP_MODULE ( count_file, sample_meta_yaml ,prefix )
    // TODO: It might not be needed to have the list of labels and files
    GUIDE_COUNT_REP_MODULE.out.out_file_map.map {
        PooledUtils.getCountByRepFilesLabel(it)}
    .multiMap{ it ->
        labels: it[0]
        count_files: it[1]
    }
    .set { label_count_file }

    rep_count_files = GUIDE_COUNT_REP_MODULE.out.rep_count_files
    rep_count_files_yml = GUIDE_COUNT_REP_MODULE.out.rep_count_files_yml
    label_list = label_count_file.labels
    label_files = label_count_file.count_files
    console_output = GUIDE_COUNT_REP_MODULE.out.console_output

    emit:
    label_list
    label_files
    rep_count_files
    rep_count_files_yml
    console_output
    versions = GUIDE_COUNT_REP_MODULE.out.versions // channel: [ versions.yml ]
}


