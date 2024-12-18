process GUIDE_COUNT_QC_MODULE {
    tag "GuideCounts"
    label 'process_low'
    label 'image_crispr_pooled_screen'

    input:
    tuple val(representation), val(count_file_label), path(count_file), path(samples_yml)
    val prefix
    path control_guides
    path annotate_dict

    output:
    path '*.count_normalized.txt', emit: normalizedtables
    path '*.guide_data_long.tsv', emit: longshapetables
    path '*.replicates_cor.tsv', emit: correlationtables
    path '*.replicates_cor.pdf', emit: pdf, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf/pooled_screen/bin/

    def args = task.ext.args ?: ''
    args += control_guides ? " --control_guides=${control_guides}" : ""
    args += annotate_dict ? " --annotate_dict=${annotate_dict}" : ""
    """
    # Making sure the labels and files match
    if [[ ! "${samples_yml}" == *"${representation}"* ]]; then
        echo "ERROR: ${representation} doesn't match ${samples_yml}"
        exit 1
    fi

    if [[ ! "${count_file}" == *"${representation}"* ]]; then
        echo "ERROR: ${representation} doesn't match ${count_file}"
        exit 1
    fi

    guide_count_QC.py \\
        --count_file="${count_file}" \\
        --representation="${representation}" \\
        --samples_yml="${samples_yml}" \\
        --prefix="${prefix}." \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
