#!/usr/bin/env nextflow
/*
========================================================================================
    PACT RNAseq batch correction
========================================================================================
*/
nextflow.enable.dsl = 2
VERSION = 0.1

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/
// Check input path parameters to see if they exist
checkPathParamList = [
    params.input
    // params.fasta, params.transcript_fasta, params.additional_fasta,
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Header log info
log.info "========================================="
log.info "Process RNAseq expression matrices v${VERSION}"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "========================================="
log.info "$PYTHONPATH"

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/
include { UPDATE_PACT_EMAT } from './workflow/update_pact_emat' addParams( params.modules["update_pact_emat"] )
include { HELLOWORLD       } from './workflow/helloworld'       addParams( params.modules["pact_emat"] )
//include { BATCH_CORRECT     } from './workflow/batch_correct'     addParams( params.modules["batch_correct"] )

Channel.from(ch_input)
    .splitCsv(sep: '\t', header: true)
    .map{ row-> tuple(
        row.sample, 
        row.gene_counts,
        row.manifest
        ) }
    .groupTuple(by: [0])
    .set { ch_input } 

ch_metadata_tsv = file(params.metadata_tsv)

// PACT RNAseq expression matrix

//pact_gid_tpm  = file(params.pact_gid_tpm)
//pact_hgnc_tpm = file(params.pact_hgnc_tpm)
//tcga_gtex_map = file(params.tcga_gtex_map)
//pact_emat_log = file(params.pact_emat_log)

ch_pact_emat = Channel.fromPath( 
    [
        file(params.pact_gid_tpm),
        file(params.pact_hgnc_tpm),
        file(params.tcga_gtex_map),
        file(params.pact_emat_log)
    ]
).collect()

workflow PROCESS_EXPRS {
    // Load workflow input
//    HELLOWORLD (
//        ch_input,
//        ch_pact_emat,
//        ch_metadata_tsv
//    )

    UPDATE_PACT_EMAT (
        ch_input,
        ch_pact_emat,
        ch_metadata_tsv
    )

    // Batch correction
//        ch_batch_correct = Channel.fromPath(
//        [
//            file(params$xena_gid_tpm),
//            file(params$bc_hgnc_tpm),
//            file(params$bc_emat_log),
//            file(params$blacklist)
//        ]
//    )
//    
//    BATCH_CORRECT (
//        ch_input,
//        ch_batch_correct,
//        UPDATE_PACT_EMAT.out.pact_emat_updated
//        )
//    
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
        subject="Process expression ${status}"
        ['aws', 'sns', 'publish', '--topic-arn', 'arn:aws:sns:us-west-2:757652839166:scrnaseq-nextflow-pipeline', '--subject', subject, '--message', message, '--region', 'us-west-2'].execute()
    }
}

//include { foo } from './workflow/foo'

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
    PROCESS_EXPRS ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
