/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPplaceeval.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SUBSETTREE     } from '../modules/local/subsettree'
include { HMMER_CLASSIFY } from '../subworkflows/local/hmmer_classify'
include { SUBSET_SUMMARY } from '../modules/local/subset_summary'
include { FINAL_SUMMARY  } from '../modules/local/final_summary'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTA_NEWICK_EPANG_GAPPA } from '../subworkflows/nf-core/fasta_newick_epang_gappa/main'

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PPLACEEVAL {

    ch_versions = Channel.empty()

    // 1. Create a number of replicates where the reference phylogeny, taxonomy and alignment is subset,
    // and a set of test sequences are output.
    Channel.fromList(0..(params.n_replicates - 1))
        .map { [ id: params.id, seed: params.seed, proportion: params.proportion, replicate: it ] }
        .combine(Channel.fromPath(params.refseqfile))
        .combine(Channel.fromPath(params.refphylogeny))
        .combine(Channel.fromPath(params.taxonomy))
        .set { ch_subset }
    SUBSETTREE ( ch_subset )
    ch_versions = ch_versions.mix(SUBSETTREE.out.versions)

    SUBSETTREE.out.ssqfasta
        .join(SUBSETTREE.out.ssrefalnfasta)
        .join(SUBSETTREE.out.ssrefnewick)
        .join(SUBSETTREE.out.sstaxonomy)
        .map { [
            meta: [ id: sprintf("%s-%05d", it[0].id, it[0].replicate) ],
            data: [
                alignmethod:  params.alignmethod ? params.alignmethod : 'hmmer',
                queryseqfile: it[1],
                refseqfile:   it[2],
                refphylogeny: it[3],
                taxonomy:     it[4],
                model:        params.model
            ]
        ] }
        .set { ch_ssplace }

    // 2. Place the test sequences in the corresponding reference phylogeny and evaluate
    FASTA_NEWICK_EPANG_GAPPA ( ch_ssplace )
    ch_versions = ch_versions.mix(FASTA_NEWICK_EPANG_GAPPA.out.versions)

    // 3. Create hmm profiles from the subset reference alignment (realign), search test set with
    // the profiles and classify by taking the best hit. Evaluate.
    // Subworkflow, i.e. a collection of module calls.
    SUBSETTREE.out.ssqfasta
        .map { [
            [ id: sprintf("%s-%05d", it[0].id, it[0].replicate) ],
            it[1]
        ] }
        .set { ch_hmmer_classify_queries }
    SUBSETTREE.out.ssclassfasta
        .map { [
            [ id: sprintf("%s-%05d", it[0].id, it[0].replicate) ],
            it[1]
        ] }
        .set { ch_hmmer_classify_classes }
    HMMER_CLASSIFY( ch_hmmer_classify_queries, ch_hmmer_classify_classes )
    ch_versions = ch_versions.mix(HMMER_CLASSIFY.out.versions)

    // 4. Compare output from 2. and 3.
    FASTA_NEWICK_EPANG_GAPPA.out.taxonomy_per_query
        .join(HMMER_CLASSIFY.out.hmmer_summary)
        .combine(Channel.fromPath(params.taxonomy))
        .map { [ it[0], it[3], it[1], it[2] ] }
        .set { ch_subset_summary }

    SUBSET_SUMMARY(ch_subset_summary)
    ch_versions = ch_versions.mix(SUBSET_SUMMARY.out.versions)

    SUBSET_SUMMARY.out.summary
        .collect { it[1] }
        .map { [ [ id: params.id ], it ] }
        .set { ch_final_summary }

    FINAL_SUMMARY(ch_final_summary)
    ch_versions = ch_versions.mix(FINAL_SUMMARY.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml')
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, [])
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
