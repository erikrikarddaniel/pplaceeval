process SUBSET_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    container "pplaceeval-r"

    input:
    tuple val(meta), path(true_clf), path(pplace_clf), path(hmmer_clf)

    output:
    tuple val(meta), path("*.clfsum.tsv.gz"), emit: summary
    path "versions.yml"           , emit: versions

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

    true_clf   <- read_tsv('$true_clf',   show_col_types = FALSE, col_names = c('protein', 'clh')) %>%
        mutate(true_class = str_remove(clh, '.*;'))

    pplace_clf <- read_tsv('$pplace_clf', show_col_types = FALSE) %>%
        rename(protein = name) %>%
        group_by(protein) %>%
        mutate(rank = rank(desc(LWR))) %>%
        ungroup() %>%
        mutate(pplace_class = str_remove(taxopath, '.*;'))
        
    hmmer_clf  <- read_tsv('$hmmer_clf',  show_col_types = FALSE) %>%
        rename(protein = target, hmmer_class = profile)

    pplace_clf %>%
        left_join(true_clf, by = 'protein', multiple = 'all') %>%
        filter(true_class == pplace_class) %>%
        count(rank) %>%
        mutate(method = 'pplace') %>%
        relocate(method) %>%
        union(
            hmmer_clf %>%
                left_join(true_clf, by = 'protein', multiple = 'all') %>%
                filter(true_class == hmmer_class) %>%
                count(rank) %>%
                mutate(method = 'hmmer') %>%
                relocate(method)
        ) %>%
        mutate(subset = '${prefix}') %>%
        relocate(subset) %>%
        write_tsv('${prefix}.clfsum.tsv.gz')

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
