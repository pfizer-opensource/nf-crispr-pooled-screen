process CONCAT_CORRELATION_TABLES {
    tag "ConcatCorrTables"
    label 'process_low'
    label 'image_crispr_pooled_screen'

    input:
    path corr_tables
    val prefix

    output:
    path '*.replicates_cor.tsv', emit: combined_table

    when:
    task.ext.when == null || task.ext.when


    shell:
    // Ensure corr_tables is a list, if we only have one
    // NOTE The `arity` parameter, introduced in NF 23.09, can be
    // used to supersede this manual check.
    if (corr_tables !instanceof List) { corr_tables = [ corr_tables ] }
    
    // Sort by file basename for stable output
    corr_tables.sort {a, b -> a.baseName <=> b.baseName }

    '''
    # Output first file in full
    cp "!{corr_tables.head()}" "!{prefix}.replicates_cor.tsv"

    # Concatenate subsequent files, without the header
    declare -a TAIL=(!{corr_tables.tail().join(" ")})
    for FILE in "${TAIL[@]}"; do
        sed 1d "$FILE" >> "!{prefix}.replicates_cor.tsv"
    done
    '''
}