process GUIDE_CORR_REP_TOTAL {
    tag "GuideCountsCorrTable"
    label 'process_low'

    container "artifacts.example.com/nextflow/functional_genomics:0.0.3"

    input:
    path (correlation_representation)
    val prefix

    output:
    path '*.replicates_cor.tsv', emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script: // This script is bundled with the pipeline, in nf/pooled_screen/bin/
    def correlation_file_list = correlation_representation.join(",")
    """
    guide_count_correlation_representation_total.py \\
        --correlation_file_list "${correlation_file_list}" \\
        --out_file "${prefix}.replicates_cor.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
