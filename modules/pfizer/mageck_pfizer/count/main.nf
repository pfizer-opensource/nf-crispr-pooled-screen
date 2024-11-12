process MAGECK_COUNT {
    tag "MageckCount"
    label 'process_low'

<<<<<<< HEAD
    container "artifacts.example.com/nextflow/mageck:0.5.9.4"
=======
    container "artifacts.example.com/functional-genomics/mageck_pfizer:latest"
>>>>>>> 05f7497 (Release v1.0)

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
    def args = task.ext.args ?: ''
    """
    mageck count \\
        --fastq ${fastqlist.join(" ")} \\
        --sample-label '${sample_label.join(",")}' \\
        --list-seq ${library} \\
        --output-prefix '${prefix}' \\
        $args

    cp .command.log mageck_count.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mageck: \$(mageck -v)
    END_VERSIONS
    """
}
