// HLA-HAT
nextflow.enable.dsl=2

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input_hlahat, params.hisat_prefix
    // params.fasta, params.transcript_fasta, params.additional_fasta,
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input_hlahat) { ch_input = file(params.input_hlahat) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//include { STAGE_FASTQS } from '../workflow/stage_fastqs'
//include { INPUT_CHECK    } from '../workflow/input_check'    addParams( options: [:] )
//include { CAT_FASTQ      } from '../workflow/input_check'    addParams( options: [:] )
include { EXTRACT_READS  } from './process/hisat_genotype'  addParams( params.modules["extract_reads"] )
include { HISAT_GENOTYPE } from './process/hisat_genotype'  addParams( params.modules["hisat_genotype"] )
//ch_hisat_prefix = Channel.fromPath(params.hisat_prefix)

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/ 

workflow HLA_HAT {
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    Channel.from(ch_input)
        .splitCsv(sep: ',', header: true)
        .map{ row-> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
        .groupTuple(by: [0])
        .set { ch_fastq }

//    INPUT_CHECK (
//        ch_input
//    )
//    .map {
//        meta, fastq ->
//            meta.id = meta.id.split('_')[0..-2].join('_')
//            [ meta, fastq ] }
//    .groupTuple(by: [0])
//    .branch {
//        meta, fastq ->
//            single  : fastq.size() == 1
//                return [ meta, fastq.flatten() ]
//            multiple: fastq.size() > 1
//                return [ meta, fastq.flatten() ]
//    }
//    .set { ch_fastq }
//    
////    //
////    // MODULE: Concatenate FastQ files from same sample if required
////    //
//    CAT_FASTQ (
//        ch_fastq.multiple
//    )
//    .mix(ch_fastq.single)
//    .view { it }
//    .set { ch_cat_fastq }

    ch_hisat_prefix = file(params.hisat_prefix)
    EXTRACT_READS (
        ch_fastq,
        ch_hisat_prefix
    )

    HISAT_GENOTYPE (
        EXTRACT_READS.out.reads,
        ch_hisat_prefix
    )
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
    HLA_HAT ()
}
