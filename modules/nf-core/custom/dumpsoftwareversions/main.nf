process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the functional_genomics container
    
    container "artifacts.example.com/nextflow/functional_genomics:0.0.3"

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}
