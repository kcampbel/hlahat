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
params.email = [:]
params.metadata = [:] // Xena, PACT metadata tsv
params.bc_tpm_hgnc = [:] // Batch corrected TPM matrix
params.pact_gid_counts = [:] // PACT raw counts matrix

checkPathParamList = [
    params.input, params.metadata, params.bc_tpm_hgnc, params.pact_gid_counts 
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Conda enviroment location
params.conda_basedir = file(params.condaprefix).getParent() 

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
include { TME } from './workflow/tme'

// Sample input
Channel.from(ch_input)
    .splitCsv(sep: '\t', header: true)
    .set { ch_input } // tuple val(meta), file(manifest), path(input_dir)

// External files input
Channel.fromPath(
    [
        params.metadata,
        params.bc_tpm_hgnc,
        params.pact_gid_counts
    ])
    .collect { it }
    .set { ch_emats }

workflow MAIN {
    TME (
        ch_input,
        ch_emats
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
    MAIN ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
