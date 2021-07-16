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

include { EXTRACT_READS     } from './process/extract_reads'     addParams( params.modules["extract_reads"] )
include { HISAT_GENOTYPE    } from './process/hisat_genotype'    addParams( params.modules["hisat_genotype"] )
include { PATIENT_REFERENCE } from './process/patient_reference' addParams( params.modules["patient_reference"] )
include { ALIGN_HLA_READS   } from './process/align_hla_reads'   addParams( params.modules["align_hla_reads"] )
include { PICARD_MARK_DUPLICATES   } from './process/mark_duplicates'   addParams( params.modules["mark_duplicates"] )

// Metadata input
include { create_fastq_channels } from './lib/functions.nf'
Channel.from(ch_input)
    .splitCsv(sep: '\t', header: true)
    .map { create_fastq_channels(it) }
    .groupTuple(by: [0])
   // .view { it }
    .set { ch_fastq }

// Global params
ch_hisat_prefix = file(params.hisat_prefix)
ch_imgthla = file(params.imgthla)

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/ 

workflow HLA_HAT {

    take:

    main:
    EXTRACT_READS (
        ch_fastq,
        ch_hisat_prefix
    )

    HISAT_GENOTYPE (
        EXTRACT_READS.out.reads,
        ch_hisat_prefix
    )

    PATIENT_REFERENCE (
        HISAT_GENOTYPE.out.hisatgt_report,
        ch_imgthla 
    )

    ch_align_hla_reads = EXTRACT_READS.out.reads
        .map {
            meta, fastqs ->
                specimen_id = meta.specimen_id
                [ specimen_id, meta, fastqs ]   
        }
        .combine(PATIENT_REFERENCE.out.patient_reference, by: 0)
        .map {
            { it [1..-1 ]}
        }
        .view {it} // [meta, fastqs, hla_fasta, hla_bed]

    ALIGN_HLA_READS (
        ch_align_hla_reads
    )

    PICARD_MARK_DUPLICATES (
        ALIGN_HLA_READS.out.hla_bam
    )

    emit:
    hisatgt_hlatypes = PATIENT_REFERENCE.out.top_hlatypes

}

workflow {
    HLA_HAT ()
}