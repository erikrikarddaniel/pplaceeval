include { MAFFT           } from '../../modules/nf-core/mafft/main'
include { HMMER_HMMBUILD  } from '../../modules/nf-core/hmmer/hmmbuild/main'
include { HMMER_HMMSEARCH } from '../../modules/nf-core/hmmer/hmmsearch/main'
include { SUMMARISEHMMER  } from '../../modules/local/summarisehmmer'

workflow HMMER_CLASSIFY {

    take:
    ch_queries // channel: [ val(meta), [ queryfasta ] ]
    ch_classes // channel: [ val(meta), [ class_fastas ] ]

    main:

    ch_versions = Channel.empty()

    ch_classes
        .transpose()
        .set { ch_mafft }
    MAFFT(ch_mafft, [])
    ch_versions = ch_versions.mix(MAFFT.out.versions.first())

    HMMER_HMMBUILD(MAFFT.out.fas, [])
    ch_versions = ch_versions.mix(HMMER_HMMBUILD.out.versions.first())

    //ch_queries.view()
    ch_queries
        .cross(HMMER_HMMBUILD.out.hmm)
        .map { [ it[0][0], it[1][1], it[0][1], true, true, [] ] }
        .set { ch_hmmsearch }

    HMMER_HMMSEARCH(ch_hmmsearch)
    ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions.first())

    HMMER_HMMSEARCH.out.target_summary
        .groupTuple()
        .set { ch_summarise }

    SUMMARISEHMMER(ch_summarise)
    ch_versions = ch_versions.mix(SUMMARISEHMMER.out.versions.first())

    emit:
    hmmer_summary = SUMMARISEHMMER.out.best_ranked

    versions      = ch_versions                     // channel: [ versions.yml ]
}

