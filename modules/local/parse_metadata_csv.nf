process PARSE_METADATA_CSV {
    tag "ParseCSVSampleMetadata"
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
        # The raw CSV input needs to be pre-processed to meet the
        # expectations of the Illumina Sample Sheet parser
        sed --regexp-extended '
            # Insert Illumina section
            1i [Samples]

            # Update template tags
            2, $ {
                s/\\{condition_(\\w+)\\}/{conditions[\\1][value]}/g
                s/\\{reference_(\\w+)\\}/{reference[\\1]}/g
            }
        ' "!{metadata}" \
        | illumina-sample-sheet-parser.py "!{task.ext.args}" \
        > "!{prefix}.parsed-metadata.yaml"

        cat <<-END_VERSIONS > versions.yml
        "!{task.process}":
            python: $(python --version | sed 's/Python //g')
        END_VERSIONS
        '''
}
