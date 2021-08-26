#!/usr/bin/env nextflow

params.rmd = [:] // tme_report.Rmd

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/
include { STAGE_INPUT } from '../process/stage_input' addParams( params.modules["stage_input_tme"] )
include { TME_REPORT  } from '../process/tme_report'  addParams( params.modules["tme_report"] )

include { create_tme_channel } from '../lib/functions.nf'
workflow TME {
    take:
        ch_input // tuple val(meta), path(manifest), path(input_folder)
        ch_emats // tuple path(metadata), path(bc_tpm_hgnc), path(pact_counts_gid)

    main:
    STAGE_INPUT (
        ch_input
    )
    STAGE_INPUT.out.tsv
      .splitCsv(sep: '\t', header: true)
      .map { create_tme_channel(it) }
      .set { ch_tme_input }

    TME_REPORT (
        ch_tme_input, // val(meta), path(*.yml)
        ch_emats
    )
}