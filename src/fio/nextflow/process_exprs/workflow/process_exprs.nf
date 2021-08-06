#!/usr/bin/env nextflow
/*
========================================================================================
    PACT RNAseq batch correction
========================================================================================
*/
// CSV input
//Channel.from(ch_input)
//    .splitCsv(sep: '\t', header: true)
//    .map{ row-> tuple(
//        row.sample, 
//        row.gene_counts,
//        row.manifest,
//        row.primary_site
//        ) }
//   // .groupTuple(by: [0])
//    .set { ch_input } 

// PACT RNAseq expression matrix
ch_pact_emat = Channel.fromPath( 
    [
        file(params.pact_gid_tpm),
        file(params.pact_hgnc_tpm),
        file(params.tcga_gtex_map),
        file(params.pact_emat_log)
    ]
).collect()

// Batch correction
ch_batch_correct = Channel.fromPath(
    [
        file(params.xena_hgnc_tpm),
        file(params.bc_emat_log),
        file(params.blacklist)
    ]
).collect()

// Pact/Xena metadata
ch_metadata_tsv = file(params.metadata_tsv)

include { UPDATE_PACT_EMAT } from './update_pact_emat' addParams( params.modules["update_pact_emat"] )
include { BATCH_CORRECT     } from './batch_correct'     addParams( params.modules["batch_correct"] )


workflow PROCESS_EXPRS {
    take:
        ch_input

    main:
    UPDATE_PACT_EMAT (
        ch_input,
        ch_pact_emat,
        ch_metadata_tsv
    )

    // Batch correction
    BATCH_CORRECT (
        ch_input,
        UPDATE_PACT_EMAT.out.pact_emat_updated, // val(meta), path(metadata_tsv_updated), path(pact_emat_log_updated), path(pact_hgnc_tpm_updated)
        ch_batch_correct
        )
    
    emit:
    bc_hgnc_tpm = BATCH_CORRECT.out.bc_hgnc_tpm
}
