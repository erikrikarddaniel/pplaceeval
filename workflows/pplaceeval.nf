/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPplaceeval.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
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

include { SUBSETTREE } from '../modules/local/subsettree'

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

/**
ch_pp_data = Channel.of([
    meta: [ id: params.id ],
    data: [
        alignmethod:  params.alignmethod ? params.alignmethod    : 'hmmer',
        queryseqfile: file(params.queryseqfile),
        refseqfile:   file(params.refseqfile),
        hmmfile:      params.hmmfile     ? file(params.hmmfile)  : [],
        refphylogeny: file(params.refphylogeny),
        model:        params.model,
        taxonomy:     params.taxonomy    ? file(params.taxonomy) : []
    ]
])
**/

    // 2. Place the test sequences in the corresponding reference phylogeny and evaluate
    //ch_ssplace.view()
    FASTA_NEWICK_EPANG_GAPPA ( ch_ssplace )
    ch_versions = ch_versions.mix(FASTA_NEWICK_EPANG_GAPPA.out.versions)

    // 3. Create hmm profiles from the subset reference alignment (realign), search test set with
    // the profiles and classify by taking the best hit. Evaluate.
    // Subworkflow, i.e. a collection of module calls.

    // 4. Compare output from 2. and 3.
    // Subworkflow, i.e. a collection of module calls.

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
