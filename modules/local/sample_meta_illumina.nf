process SAMPLE_META_MODULE_ILLUM {
    tag "GenerateSampleMeta"
    label 'process_low'

    container "artifacts.example.com/nextflow/functional_genomics:0.0.3"

    input:
    path fastq_dir
    path metadata
    val prefix

    output:
    path "*.metadata.yaml"   , emit: sample_meta_yaml
    path "versions.yml"      , emit: versions

    script: // This script is bundled with the pipeline, in nf/pooled_screen/bin/
    """
    illum_meta_generation.py \\
        --fastq_dir="${fastq_dir}" \\
        --metadata="${metadata}" \\
        --output="${prefix}.metadata.yaml"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
