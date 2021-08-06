#!/usr/bin/env nextflow

params.metadata = [:] // Xena, PACT metadata tsv
params.bctpm = [:] // Batch corrected TPM matrix
params.counts = [:] // PACT raw counts matrix

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/
include { STAGE_INPUT } from '../process/stage_input' addParams( params.modules["tme_report"] )
include { TME_REPORT  } from '../process/tme_report'  addParams( params.modules["tme_report"] )

// External files input
Channel.fromList(
    [
        file(params.metadata),
        file(params.bctpm),
        file(params.counts)
    ])
    .collect { it }
    .set { ch_bctpm }

workflow TME {
    take:
        ch_input

    main:
    STAGE_INPUT (
        ch_input
    )
    TME_REPORT (
        STAGE_INPUT.out.tme_config, // val(meta), path(*.yml)
        ch_bctpm
    )
}