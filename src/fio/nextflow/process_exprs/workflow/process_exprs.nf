#!/usr/bin/env nextflow
/*
========================================================================================
    PACT Expression Matrix Updater and Batch Correction
========================================================================================
*/
params.pact_tcga_gtex_map = [:]

// PACT RNAseq expression matrix
ch_pact_emat = Channel.fromPath( 
    [
        params.pact_gid_counts,
        params.pact_gid_tpm,
        params.pact_hgnc_tpm,
        params.tcga_gtex_map,
        params.pact_emat_log
    ]
).collect()

// Batch correction
ch_batch_correct = Channel.fromPath(
    [
        params.xena_hgnc_tpm,
        params.bc_emat_log,
        params.blacklist
    ]
).collect()

// Pact/Xena metadata
ch_metadata_tsv = file(params.metadata_tsv)

include { STAGE_INPUT      } from '../process/stage_input'      addParams( params.modules["stage_input"] )
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
      .set { ch_pe_input }

    UPDATE_PACT_EMAT (
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
    emats = BATCH_CORRECT.out.bc_hgnc_tpm_updated
      .join(UPDATE_PACT_EMAT.out.pact_gid_counts_updated)
      .collect { it[3,2,-1] } // tuple path(metadata), path(bc_tpm_hgnc), path(pact_counts_gid)
    pe_input = ch_input // meta, bc_log_updated, bc_tpm_updated, metadata_tsv  
}
