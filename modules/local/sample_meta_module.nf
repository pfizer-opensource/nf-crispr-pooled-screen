process SAMPLE_META_MODULE {
    tag "GenerateSampleMeta"
    label 'process_low'

    container "artifacts.example.com/nextflow/functional_genomics:0.0.3"

    input:
    path fastq_dir
    path input_sample_meta_table
    val prefix

    output:
    path "*sample_meta.yml"   , emit: sample_meta_yaml
    path "versions.yml"      , emit: versions

    script: // This script is bundled with the pipeline, in nf/pooled_screen/bin/
    """
    sample_meta_generation.py \\
        --fastq_dir=${fastq_dir} \\
        --sample_input_meta_table=${input_sample_meta_table} \\
        --output=${prefix}sample_meta.yml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
