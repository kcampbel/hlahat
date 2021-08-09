#!/usr/bin/env nextflow

params.metadata             = [:] // Xena, PACT metadata tsv
params.bc_tpm_hgnc     = [:] // Batch corrected TPM matrix
params.pact_counts_gid = [:] // PACT raw counts matrix
params.rmd = [:] // tme_report.Rmd

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/
include { STAGE_INPUT } from '../process/stage_input' addParams( params.modules["tme_report"] )
include { TME_REPORT  } from '../process/tme_report'  addParams( params.modules["tme_report"] )

// External files input

//Channel.fromList(
//    [
//        file(params.metadata),
//        file(params.bc_tpm_hgnc),
//        file(params.pact_counts_gid)
//    ])
//    .collect { it }
//    .set { ch_bctpm }

workflow TME {
    take:
        ch_input
        ch_emats

    main:
    STAGE_INPUT (
        ch_input
    )
    TME_REPORT (
        STAGE_INPUT.out.tme_config, // val(meta), path(*.yml)
        ch_emats
    )
}