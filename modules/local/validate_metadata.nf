process VALIDATE_METADATA {
    tag "ValidateMetadata"
    label 'process_low'
    label 'image_crispr_pooled_screen'

    input:
        path metadata_yaml
        path fastq_dir
        val prefix

    output:
        path "${prefix}.metadata.yaml", emit: meta_yaml
        path "versions.yml",            emit: versions

    shell:
        '''
        # The parsed metadata YAML input needs to be pre-processed to
        # update the legacy metadata template tags with those that will
        # work with the new pipeline
        sed --regexp-extended '
            s/\\{condition_(\\w+)\\}/{conditions[\\1][value]}/g
            s/\\{reference_(\\w+)\\}/{reference[\\1]}/g
        ' "!{metadata_yaml}" \
        | metadata-validator.py "!{fastq_dir}" \
        > "!{prefix}.metadata.yaml"

        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            python: $(python --version | sed 's/Python //g')
        END_VERSIONS
        '''
}
