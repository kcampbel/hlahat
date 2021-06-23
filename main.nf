#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/rnaseq
========================================================================================
    Github : https://github.com/nf-core/rnaseq
    Website: https://nf-co.re/rnaseq
    Slack  : https://nfcore.slack.com/channels/rnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================

params.fasta        = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf          = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff          = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed     = WorkflowMain.getGenomeAttribute(params, 'bed12')
params.star_index   = WorkflowMain.getGenomeAttribute(params, 'star')
params.hisat2_index = WorkflowMain.getGenomeAttribute(params, 'hisat2')
params.rsem_index   = WorkflowMain.getGenomeAttribute(params, 'rsem')
params.salmon_index = WorkflowMain.getGenomeAttribute(params, 'salmon')
/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

//params.input_hlahat = false
params.email = null

include { HLA_HAT } from './hlahat/workflow/hla_hat'
include { HLA_CN } from './hlacn/workflow/hla_cn'
workflow IMMUNE_ESCAPE {
    HLA_HAT ()
    HLA_CN (
        hisatgt_hlatypes = HLA_HAT.out.hisatgt_hlatypes
    )
}

process foo {
    output:
      path 'foo.txt'
      stdout emit: verbiage
    script:
      """
      echo 'test' > foo.txt
      echo "thisis a test"
      """
}

workflow.onComplete {
    log.info("=======================================")
    status="${ workflow.success ? 'OK' : 'FAILED' }"
    message="Pipeline completed at: ${workflow.complete}\nExecution status: ${status}\nWorkdir: ${workflow.workDir}\nPublish dir: ${params.outdir}"
    log.info(message)

    // Email 
    if ( params.email ) {
        log.info("Emailing ${params.email}")
        subject="Immune Escape ${status}"
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
    IMMUNE_ESCAPE ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
