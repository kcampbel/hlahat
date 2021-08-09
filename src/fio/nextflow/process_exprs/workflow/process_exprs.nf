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

params.pact_tcga_gtex_map = [:]

// PACT RNAseq expression matrix
ch_pact_emat = Channel.fromPath( 
    [
        file(params.pact_gid_counts),
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

include { STAGE_INPUT      } from '../process/stage_input_pe'      addParams( params.modules["stage_input"] )
include { UPDATE_PACT_EMAT } from '../process/update_pact_emat' addParams( params.modules["update_pact_emat"] )
include { BATCH_CORRECT    } from '../process/batch_correct'    addParams( params.modules["batch_correct"] )

include { create_pe_channel } from '../lib/functions.nf'
workflow PROCESS_EXPRS {
    take:
        ch_input

    main:
    STAGE_INPUT (
        ch_input
    )
    STAGE_INPUT.out.tsv
      .splitCsv(sep: '\t', header: true)
      .map { create_pe_channel(it) }
    //.groupTuple(by: [0])
    //  .view { it }
      .set { ch_pe_input }

    UPDATE_PACT_EMAT (
        //ch_input,
        ch_pe_input,
        ch_pact_emat,
        ch_metadata_tsv
    )

    // Batch correction
    BATCH_CORRECT (
        ch_pe_input,
        UPDATE_PACT_EMAT.out.pact_hgnc_tpm_updated, // tuple val(meta), path(metadata_tsv_updated), path(pact_hgnc_tpm_updated)
        ch_batch_correct
        )
    
    emit:
    pact_gid_counts = UPDATE_PACT_EMAT.out.pact_gid_counts_updated // tuple val(meta), path(metadata_tsv_updated), path(pact_gid_counts_updated)
    bc_hgnc_tpm = BATCH_CORRECT.out.bc_hgnc_tpm_updated // val(meta), path(bc_log_updated), path(bc_tpm_updated)
}
