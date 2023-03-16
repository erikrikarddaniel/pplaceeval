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
    tuple val(meta), path("${meta.id}.${sprintf('%05d', meta.replicate)}.ref.alnfasta")  , emit: ssrefalnfasta
    tuple val(meta), path("${meta.id}.${sprintf('%05d', meta.replicate)}.query.alnfasta"), emit: ssqalnfasta
    tuple val(meta), path("${meta.id}.${sprintf('%05d', meta.replicate)}.ref.newick")    , emit: ssrefnewick
    tuple val(meta), path("${meta.id}.${sprintf('%05d', meta.replicate)}.ref.tsv")       , emit: sstaxonomy
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = sprintf("%s.%05d", meta.id, meta.replicate)

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(treeio))
    suppressPackageStartupMessages(library(readr))
    suppressPackageStartupMessages(library(dplyr))

    set.seed($meta.seed + $meta.replicate)

    aln     <- read.fasta("$aln")
    tree    <- read.newick("$phylo")
    tax     <- read_tsv("$taxonomy", col_names = c('seqname', 'class'), show_col_types = FALSE)

    labels  <- tree\$tip.label
    qlabels <- sample(labels, round(length(labels) * $meta.proportion))

    sstree  <- drop.tip(tree, qlabels)
    ssaln   <- aln[sstree\$tip.label]
    sstax   <- tax %>% filter(! seqname %in% qlabels)
    qaln    <- aln[qlabels]

    ape::write.FASTA(ssaln, "${prefix}.ref.alnfasta")
    ape::write.FASTA(qaln, "${prefix}.query.alnfasta")
    write.tree(sstree, "${prefix}.ref.newick")
    write.table(sstax, "${prefix}.ref.tsv", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

    # Output versions.yml
    writeLines(
        c(
            "\\"${task.process}\\":", 
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    treeio: ", packageVersion("treeio")),
            paste0("    ape: ", packageVersion("ape")),
            paste0("    tidyverse: ", packageVersion("tidyverse")) 
        ), 
        "versions.yml"
    )
    """
}
