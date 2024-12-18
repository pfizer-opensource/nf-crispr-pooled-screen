process PARSE_METADATA_ILLUMINA {
    tag "ParseIlluminaMetadata"
    label 'process_low'
    label 'image_crispr_pooled_screen'

    input:
        path metadata
        val prefix

    output:
        path "${prefix}.parsed-metadata.yaml", emit: meta_yaml
        path "versions.yml",                   emit: versions

    shell:
        '''
        < "!{metadata}" \
        illumina-sample-sheet-parser.py "!{task.ext.args}" \
        > "!{prefix}.parsed-metadata.yaml"

        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            python: $(python --version | sed 's/Python //g')
        END_VERSIONS
        '''
}
