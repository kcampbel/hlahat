#!/usr/bin/env nextflow
/*
========================================================================================
    PACT Expression Matrix Updater and Batch Correction
========================================================================================
*/
nextflow.enable.dsl = 2
VERSION = 0.1

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/
// Check input path parameters to see if they exist
params.email = [:]
params.tcga_gtex_map = [:]

checkPathParamList = [
    params.input, params.tcga_gtex_map
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Conda enviroment location
params.conda_basedir = file(params.condaprefix).getParent() 

// Header log info
log.info "========================================="
log.info "Process RNAseq expression matrices v${VERSION}"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "========================================="

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/
include { PROCESS_EXPRS } from './workflow/process_exprs'

// CSV input
Channel.from(ch_input)
    .splitCsv(sep: '\t', header: true)
    .set { ch_input } // val(meta), file(manifest), file(input_dir)

workflow MAIN {
    PROCESS_EXPRS (
        ch_input
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
        subject="Process expression ${status}"
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
