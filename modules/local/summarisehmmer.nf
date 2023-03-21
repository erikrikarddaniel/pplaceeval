process SUMMARISEHMMER {
    tag "$meta.id"
    label 'process_single'

    container "pplaceeval-r"

    input:
    tuple val(meta), path(tblouts)

    output:
    tuple val(meta), path("*.hmmclassify.tsv.gz"), emit: best_ranked
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(readr))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(tidyr))
    suppressPackageStartupMessages(library(stringr))

    Sys.glob('*.tbl.gz') %>%
        read_tsv(col_names = c('c'), show_col_types = FALSE, comment = '#') %>%
        separate(
            c, 
            c(
                'target', 'f0', 'profile', 'f1', 'evalue', 'score', 'bias', 
                'domevalue', 'domscore', 'dombias', 
                'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10'
            ), 
            sep = ' +'
        ) %>%
        transmute(
            target, profile = str_remove(profile, '.*class.'),
            score = as.double(score), evalue = as.double(evalue)
        ) %>%
        group_by(target) %>%
        arrange(desc(score), evalue, profile) %>%
        mutate(rank = row_number()) %>%
        ungroup() %>%
        write_tsv('${prefix}.hmmclassify.tsv.gz')

    # Output versions.yml
    writeLines(
        c(
            "\\"${task.process}\\":", 
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    treeio: ", packageVersion("treeio")),
            paste0("    Biostrings: ", packageVersion("Biostrings")),
            paste0("    tidyverse: ", packageVersion("tidyverse")) 
        ), 
        "versions.yml"
    )
    """
}
