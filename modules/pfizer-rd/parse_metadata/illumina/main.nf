process PARSE_METADATA_ILLUMINA {
    tag "ParseMetadataIllumina"
    label 'process_low'

    container "artifacts.example.com/functional-genomics/parse_metadata:latest"

    input:
        path metadata_illumina

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
            --sample-sheet "!{metadata_illumina}"
            --output "parsed-metadata.yaml"
        )

        if [[ "!{allow_ambiguous_schema}" == "true" ]]; then
            ARGS+=(--allow-ambiguous-schema)
        fi

        if [[ "!{debug}" == "true" ]]; then
            export DEBUG=1
        fi

        illumina-sample-sheet-parser "${ARGS[@]}" "!{schema}" 2>&1 \
        | tee "parse-metadata.log"

        cat <<-VERSIONS >versions.yml
        "!{task.process}":
            python: $(python --version | sed "s/Python //")
            illumina-sample-sheet-parser: $(illumina-sample-sheet-parser --version)
        VERSIONS
        '''
}
