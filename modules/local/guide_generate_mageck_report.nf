process MAGECK_REPORT_HTML {
    tag "MageckReport"
    label 'process_low'

    container "artifacts.example.com/nextflow/mageck_report_html:0.0.3"

    input:
    path report_rmd_file
    val prefix

    output:
    path '*.html', emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf/pooled_screen/bin/
    """
    # the render command seems like not to like links, which is how nextflow
    # works. So we need to resolve the path to the report file
    export report=\$(realpath ${report_rmd_file})
    Rscript -e "library(ggplot2);\\
                library(tidyr);\\
                rmarkdown::render(\\"\${report}\\",
                    output_file=\\"${prefix}mageck_report.html\\",
                    output_dir=getwd(),
                    quiet=TRUE)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | sed -n '1s/.* \\([0-9].[0-9].[0-9]\\).*/1/p')
    END_VERSIONS
    """
}
