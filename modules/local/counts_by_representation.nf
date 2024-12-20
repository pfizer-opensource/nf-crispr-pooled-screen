process COUNTS_BY_REPRESENTATION {
    tag "CountRepresentations"
    label 'process_low'
    label 'image_crispr_pooled_screen'

    input:
    tuple val(meta), path(count_file)
    path metadata_yaml
    val prefix

    output:
    stdout emit: console_output
    path '*.output_file_map.yml', emit: repr_name_map
    path "*.count.txt", emit: counts
    path "*.count.yml", emit: meta_yamls
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf/pooled_screen/bin/
    """
    guide_count_by_representation.py \\
        --count_file="${count_file}" \\
        --sample_meta_yaml="${metadata_yaml}" \\
        --prefix="${prefix}." \\
        --output_file_map="${prefix}.output_file_map.yml"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
