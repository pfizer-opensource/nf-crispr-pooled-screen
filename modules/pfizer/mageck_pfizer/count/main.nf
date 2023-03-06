process MAGECK_COUNT {
    tag "MageckCount"
    label 'process_low'

    container "artifacts.example.com/nextflow/mageck:0.5.9.4"

    input:
    path fastqfolder
    val  sample_label
    val  fastqlist
    path library
    val  prefix

    output:
    path "*.count.txt", emit: count
    path "*.count_normalized.txt", emit: norm
    path "*.countsummary.txt", emit: summary
    path "*.count_report.Rmd", emit: report
    path "versions.yml"           , emit: versions
    path "mageck_count.log", emit: log


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = params.mageck_count_options ?: ''
    def norm_mode = params.mageck_count_normalization_method ?: "none"
    """
    mageck count \\
        --fastq ${fastqlist.join(" ")} \\
        --sample-label '${sample_label.join(",")}' \\
        --list-seq ${library} \\
        --output-prefix '${prefix}' \\
        --norm-method ${norm_mode} \\
        $args

    cp .command.log mageck_count.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mageck: \$(mageck -v)
    END_VERSIONS
    """
}
