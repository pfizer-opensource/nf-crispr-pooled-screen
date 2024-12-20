process PARSE_METADATA_CSV {
    tag "ParseMetadataCSV"
    label 'process_low'

    container "artifacts.example.com/functional-genomics/parse_metadata:latest"

    input:
        path metadata_csv

    output:
        path "parsed-metadata.yaml", emit: metadata_yaml

        path "versions.yml",         emit: versions
        path "parse-metadata.log",   emit: log

    shell:
        schema = "${projectDir}/${task.ext.schema}"
        allow_ambiguous_schema = task.ext.allow_ambiguous_schema ?: false
        debug = task.ext.debug ?: false

        '''
        set -o pipefail

        declare -a ARGS=(
            --output "parsed-metadata.yaml"
        )

        if [[ "!{allow_ambiguous_schema}" == "true" ]]; then
            ARGS+=(--allow-ambiguous-schema)
        fi

        if [[ "!{debug}" == "true" ]]; then
            export DEBUG=1
        fi

        {
            # Use [yq] to extract the section mapping from the schema
            # file; for raw CSV data, we only need one section, so go
            # for the first one present in the schema. This way, we can
            # convert the raw CSV input into an Illumina Sample Sheet
            # that is conducive to the supplied schema and reuse the
            # normal Illumina parser.
            #
            # [yq]: https://mikefarah.gitbook.io/yq

            yq eval --unwrapScalar '"[\\(to_entries[0].value | (.section.name // .section))]"' "!{schema}"
            cat "!{metadata_csv}"
        } \
        | illumina-sample-sheet-parser "${ARGS[@]}" "!{schema}" 2>&1 \
        | tee "parse-metadata.log"

        cat <<-VERSIONS >versions.yml
        "!{task.process}":
            python: $(python --version | sed "s/Python //")
            illumina-sample-sheet-parser: $(illumina-sample-sheet-parser --version)
        VERSIONS
        '''
}
