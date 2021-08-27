#!/usr/bin/env nextflow
/*
========================================================================================
    HLA Allele Copy Number Pipeline
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

checkPathParamList = [
    params.input, params.hisat2_hlatypes, params.genome_fasta
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.hisat2_hlatypes) { ch_hisat2_hlatypes = file(params.hisat2_hlatypes) } else { exit 1, 'Hisat2 hlatypes not specified!' }
if (params.genome_fasta) { ch_genome_fasta = file(params.genome_fasta) } else { exit 1, 'Genome fasta not specified!' }

// Conda enviroment location
params.conda_basedir = file(params.condaprefix).getParent() 

// Header log info
log.info "========================================="
log.info "HLA Allele Copy Number Pipeline v${VERSION}"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "========================================="

include { HLACN } from './workflow/hlacn'

// Sample input tsv
Channel.from(ch_input)
    .splitCsv(sep: '\t', header: true)
    .set { ch_input } // val(specimen_id), file(manifest), file(input_dir)

// Hisat2 genotype input tsv
Channel.from(ch_hisat2_hlatypes)
    .splitCsv(sep: '\t', header: true)
    .map { row -> [ row.specimen_id, [ file(row.hisat2_hlatypes)] ]}
    .set { ch_hisat2_hlatypes } // val(specimen_id), file(hisat2_hla_types)

workflow MAIN {
    HLACN (
        ch_input,
        ch_hisat2_hlatypes
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
        subject="HLACN ${status}"
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
