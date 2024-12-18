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
        < "!{metadata_yaml}" \
        metadata-validator.py "!{fastq_dir}" \
        > "!{prefix}.metadata.yaml"

        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            python: $(python --version | sed 's/Python //g')
        END_VERSIONS
        '''
}
