process SUBSETTREE {
    tag "$meta.id:${sprintf('%05d', meta.replicate)}"
    label 'process_single'

    // Ghada: This is unusual. A normal nf-core pipeline would have biocontainer containers,
    // but this is difficult with R since we need many libraries. In nf-core pipelines one
    // would ask the biocontainer people for a "mulled" container but that takes time (ours
    // as well as other peoples'). I hence decided to create a special Docker container for
    // our R modules. We can add whatever we need to it. Files are in images/pplaceeval-r.
    container "pplaceeval-r"

    input:
    tuple val(meta), path(aln), path(phylo), path(taxonomy) // meta is assumed to contain .id (a name), .seed (random seed), proportion (proportion of labels to extract) and .replicate (a number)

    output:
    tuple val(meta), path("${meta.id}.${sprintf('%05d', meta.replicate)}.ref.alnfasta") , emit: ssrefalnfasta
    tuple val(meta), path("${meta.id}.${sprintf('%05d', meta.replicate)}.query.fasta")  , emit: ssqfasta
    tuple val(meta), path("${meta.id}.${sprintf('%05d', meta.replicate)}.ref.newick")   , emit: ssrefnewick
    tuple val(meta), path("${meta.id}.${sprintf('%05d', meta.replicate)}.ref.tsv")      , emit: sstaxonomy
    tuple val(meta), path("${meta.id}.${sprintf('%05d', meta.replicate)}.class.*.fasta"), emit: ssclassfasta
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = sprintf("%s.%05d", meta.id, meta.replicate)

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(treeio))
    suppressPackageStartupMessages(library(Biostrings))
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(readr))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(tidyr))

    set.seed($meta.seed + $meta.replicate)

    write_fasta <- function(d, f) {
        d %>%
            transmute(s = sprintf(">%s\\n%s", name, sequence)) %>%
            write.table(f, col.names = FALSE, row.names = FALSE, quote = FALSE)
    }

    aln     <- readAAMultipleAlignment("$aln") %>%
        as('AAStringSet') %>%
        data.frame() %>%
        tibble::rownames_to_column('name')
    colnames(aln) <- c('name', 'sequence')

    tree    <- read.newick("$phylo")
    tax     <- read_tsv("$taxonomy", col_names = c('seqname', 'class'), show_col_types = FALSE)

    labels  <- tree\$tip.label
    qlabels <- sample(labels, round(length(labels) * $meta.proportion))

    sstree  <- drop.tip(tree, qlabels)
    ssaln   <- aln %>% filter(name %in% sstree\$tip.label)
    sstax   <- tax %>% filter(! seqname %in% qlabels)
    qseq    <- aln %>% filter(name %in% qlabels) %>% mutate(sequence = str_remove_all(sequence, '-'))

    write_fasta(ssaln, "${prefix}.ref.alnfasta")
    write_fasta(qseq, "${prefix}.query.fasta")
    write.tree(sstree, "${prefix}.ref.newick")
    write.table(sstax, "${prefix}.ref.tsv", row.names = FALSE, col.names = FALSE, sep = '\\t', quote = FALSE)

    # Write unaligned fasta files for classes
    sstax   <- sstax %>% separate(class, c('l0', 'l1', 'l2', 'l3', 'l4'), sep = ';', fill = 'right')
    levels  <- tibble(level = character(), name = character())
    g <- function(tab, level) {
        tab %>%
            filter(!is.na(!!sym(level))) %>%
            group_by(level = !!sym(level)) %>%
            summarise(name = list(seqname)) %>%
            mutate(level = level) %>%
            unnest(name)
    }
    for ( l in c('l0', 'l1', 'l2', 'l3', 'l4') ) {
        levels <- dplyr::union(levels, g(sstax, l))
    }
    seqsets <- levels %>%
        inner_join(ssaln %>% mutate(sequence = str_remove_all(sequence, '-')), by = 'name')

    for ( s in seqsets %>% distinct(level) %>% pull(level) ) {
        write_fasta(
            seqsets %>% filter(level == s) %>% select(name, sequence),
            sprintf("${prefix}.class.%s.fasta", s)
        )
    }

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
