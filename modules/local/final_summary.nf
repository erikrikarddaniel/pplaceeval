process FINAL_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    container "pplaceeval-r"

    input:
    tuple val(meta), path(ss_summaries)

    output:
    tuple val(meta), path("*.overallsum.tsv"), emit: summary
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(readr))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(tidyr))
    suppressPackageStartupMessages(library(stringr))

    c('${ss_summaries.join("', '")}') %>%
        read_tsv(show_col_types = FALSE) %>%
        group_by(method, rank) %>%
        summarise(mean = mean(n), .groups = 'drop') %>%
        pivot_wider(names_from = c(method, rank), values_from = mean) %>%
        write_tsv('${prefix}.overallsum.tsv')

    # Output versions.yml
    writeLines(
        c(
            "\\"${task.process}\\":", 
            paste0("    tidyverse: ", packageVersion("tidyverse")) 
        ), 
        "versions.yml"
    )
    """
}
