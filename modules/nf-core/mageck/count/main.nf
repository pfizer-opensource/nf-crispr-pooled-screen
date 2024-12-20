process MAGECK_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mageck:0.5.9.5--py39h1f90b4d_3':
        'biocontainers/mageck:0.5.9.5--py39h1f90b4d_3' }"

    input:
    tuple val(meta), path(inputfile)
    path(library)
    path(control_sgrna)

    output:
    tuple val(meta), path("*.count.txt")           , emit: count
    tuple val(meta), path("*.count_normalized.txt"), emit: norm
    tuple val(meta), path("*.countsummary.txt")    , emit: summary
    tuple val(meta), path("*.count_report.Rmd")    , emit: report
    tuple val(meta), path("*.log")                 , emit: logs
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    args += control_sgrna ? " --control-sgrna=${control_sgrna}" : ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def is_fastq = "$inputfile".endsWith(".fastq.gz") || "$inputfile".endsWith(".fq.gz")
    if (is_fastq) {
        // For FASTQ input, `inputfile` is a list of FASTQ files and `meta.id` is a
        // comma-separated string of sample names in the same order. FASTQ files with the
        // same sample name (ie, technical replicates) are represented as a comma-
        // separated string of filenames, then filenames for different samples are space-
        // separated. `meta.id` is reduced to a comma-separated string of unique names.
        def samples = [meta.id.split(','), inputfile]
            .transpose()
            .groupBy { it[0] }
            .collectEntries { k, v -> [(k): v.flatten().tail().join(",")] }
        inputfile = samples.values().join(" ")
        meta.id = samples.keySet().join(",")
    }
    def input_file = ("$inputfile".endsWith(".fastq.gz") || "$inputfile".endsWith(".fq.gz")) ? "--fastq ${inputfile}" : "-k ${inputfile}"
    def sample_label = ("$inputfile".endsWith(".fastq.gz") || "$inputfile".endsWith(".fq.gz")) ? "--sample-label ${meta.id}" : ''

    """
    mageck \\
        count \\
        $args \\
        -l $library \\
        -n $prefix \\
        $sample_label \\
        $input_file \\

    mv "${prefix}.log" mageck_count.log


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mageck: \$(mageck -v)
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_file = ("$inputfile".endsWith(".fastq.gz")) ? "--fastq ${inputfile}" : "-k ${inputfile}"
    def sample_label = ("$inputfile".endsWith(".fastq.gz") || "$inputfile".endsWith(".fq.gz")) ? "--sample-label ${meta.id}" : ''
    """
    touch ${prefix}.count.txt
    touch ${prefix}.count_normalized.txt
    touch ${prefix}.countsummary.txt
    touch ${prefix}.count_report.Rmd
    touch ${prefix}.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mageck: \$(mageck -v)
    END_VERSIONS
    """
}
