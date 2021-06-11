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

include { GENOTYPER } from './workflow/genotyper'

workflow HLA_HAT {
    GENOTYPER ()
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
    HLA_HAT ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
