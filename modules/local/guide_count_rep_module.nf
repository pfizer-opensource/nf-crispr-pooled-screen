process GUIDE_COUNT_REP_MODULE {
    tag "CountRepresentations"
    label 'process_low'
    label 'image_crispr_pooled_screen'

    input:
    path count_file
    path sample_meta_yaml
    val prefix

    output:
    stdout emit: console_output
    path '*.output_file_map.yml', emit: out_file_map
    path "*.count.txt", emit: rep_count_files
    path "*.count.yml", emit: rep_count_files_yml
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf/pooled_screen/bin/
    """
    guide_counts_rep.py \\
        --count_file="${count_file}" \\
        --sample_meta_yaml="${sample_meta_yaml}" \\
        --prefix="${prefix}." \\
        --output_file_map="${prefix}.output_file_map.yml"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
