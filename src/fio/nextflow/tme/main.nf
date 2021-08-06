#!/usr/bin/env nextflow
/*
========================================================================================
    Tumor Microenvironment Report
========================================================================================
*/
nextflow.enable.dsl = 2
PIPELINE_NAME = 'TME Report'
VERSION = 0.1

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input 
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
params.email = [:]
params.metadata = [:] // Xena, PACT metadata tsv
params.bctpm = [:] // Batch corrected TPM matrix
params.counts = [:] // PACT raw counts matrix
//params.metadata = '/Users/csmith/git/bioinfo-fio/tme/test_data/meta_pact_xena.tsv'
//params.counts = '/Users/csmith/git/bioinfo-fio/tme/test_data/pact_rsem_raw.exprs.gid.tsv.gz'
//params.bctpm = '/Users/csmith/Documents/References/xena/subsets/PactTcgaTargetGtex_rsem_gene_tpm/test_combat_ColonRectum_hgnc.tsv.gz'

// Header log info
log.info "========================================="
log.info "Tumor microenvironment Report v${VERSION}"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "========================================="

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/
include { STAGE_INPUT } from './process/stage_input' addParams( params.modules["tme_report"] )
include { TME_REPORT } from './process/tme_report' addParams( params.modules["tme_report"] )

// Sample input
include { create_tme_channel } from './tme/lib/functions.nf'
Channel.from(ch_input)
    .splitCsv(sep: '\t', header: true)
    .map { create_tme_channel(it) }
    //.groupTuple(by: [0])
    //.view { it }
    .set { ch_input }

// External files input
Channel.fromList(
    [
        file(params.metadata),
        file(params.bctpm),
        file(params.counts)
    ])
    //.view { it }
    .collect { it }
    .set { ch_bctpm }

workflow TME {
   // STAGE_INPUT (
   //     ch_input
   // )
    TME_REPORT (
        STAGE_INPUT.out.tme_config, // val(meta), path(*.yml)
        ch_bctpm
    )
}

workflow.onComplete {
    log.info("=======================================")
    status="${ workflow.success ? 'OK' : 'FAILED' }"
    message="Pipeline completed at: ${workflow.complete}\nExecution status: ${status}\nWorkdir: ${workflow.workDir}\nPublish dir: ${params.outdir}"
    log.info(message)

    // Email 
    if ( params.email ) {
        log.info("Emailing ${params.email}")
        subject="${PIPELINE_NAME} ${status}"
        ['aws', 'sns', 'publish', '--topic-arn', 'arn:aws:sns:us-west-2:757652839166:scrnaseq-nextflow-pipeline', '--subject', subject, '--message', message, '--region', 'us-west-2'].execute()
    }
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    TME ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
