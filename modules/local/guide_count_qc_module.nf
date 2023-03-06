process GUIDE_COUNT_QC_MODULE {
    tag "GuideCounts"
    label 'process_low'

    container "artifacts.example.com/nextflow/functional_genomics:0.0.3"

    input:
    tuple val(labels_list), val(count_files_list), path(rep_count_files), path(rep_count_files_yml)
    val prefix
    path annotate_dict



    output:
    path '*rep_correlation_*.tsv', emit: correlationtables
    path '*normalized_count_*.tsv', emit: normalizedtables
    path '*data_long_shape_*.tsv', emit: longshapetables
    path '*.pdf', emit: pdf, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when
    // check if annotate_dict is not null and define args

    script: // This script is bundled with the pipeline, in nf/pooled_screen/bin/

    def args = annotate_dict ? "--annotate_dict=${annotate_dict}" : ""
    """
    # Making sure the labels and files match
    if [[ ! "${rep_count_files_yml}" == *"${labels_list}"* ]]; then
        echo "ERROR: ${labels_list} doesn't match ${rep_count_files_yml}"
        exit 1
    fi

    if [[ ! "${rep_count_files}" == *"${labels_list}"* ]]; then
        echo "ERROR: ${labels_list} doesn't match ${rep_count_files}"
        exit 1
    fi

    guide_count_QC.py \\
        --repr_lis=${labels_list} \\
        --rep_count_files=${rep_count_files} \\
        --rep_count_files_yml=${rep_count_files_yml} \\
        --prefix=${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
