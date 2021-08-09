#!/usr/bin/env nextflow
/*
========================================================================================
    PACT FIO Pipeline
========================================================================================
*/
nextflow.enable.dsl = 2
PIPELINE_NAME = 'FIO'
VERSION = 0.1
/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/
params.conda_basedir = file(params.condaprefix).getParent() 

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input 
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
params.email = [:]

// Header log info
log.info "========================================="
log.info "FIO v${VERSION}"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "========================================="

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/
include { PROCESS_EXPRS } from './process_exprs/workflow/process_exprs'
include { TME           } from './tme/workflow/tme'

Channel.from(ch_input)
    .splitCsv(sep: '\t', header: true)
    .set { ch_input }

// PROCESS_EXPRS
//Channel.from(ch_input)
//    .splitCsv(sep: '\t', header: true)
//    .map{ row-> tuple(
//        row.sample, 
//        row.manifest,
//        row.primary_site,
//        row.gene_counts
//        ) }
//   // .groupTuple(by: [0])
//    .set { ch_process_exprs } 


// TME
//Channel.from(ch_input)
//    .splitCsv(sep: '\t', header: true)
// //   .map { create_tme_channel(it) }
//    //.groupTuple(by: [0])
//    //.view { it }
//    .set { ch_tme }

workflow FIO {
    PROCESS_EXPRS (
        ch_input
        //ch_process_exprs
    )

//    bc_hgnc_tpm = PROCESS_EXPRS.out.bc_hgnc_tpm  
//    pact_gid_counts = PROCESS_EXPRS.out.pact_gid_counts
//    Channel.fromList(
//        [
//            file(bc_hgnc_tpm[0]), // metadata
//            file(bc_hgnc_tpm[-1]), // batch corrected tpms emat
//            file(pact_gid_counts[-1]) // PACT raw count emat 
//        ])
//        .collect { it }
//        .set { ch_emats }
//    TME (
//        ch_tme,
//        ch_emats
//    )
}

workflow.onComplete {
    log.info("=======================================")
    status="${ workflow.success ? 'OK' : 'FAILED' }"
    message="Pipeline completed at: ${workflow.complete}\nExecution status: ${status}\nWorkdir: ${workflow.workDir}\nPublish dir: ${params.outdir}"
    log.info(message)

    // Email 
    if ( params.email ) {
        log.info("Emailing ${params.email}")
        subject="FIO ${status}"
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
    FIO ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
