#!/usr/bin/env nextflow
/*
========================================================================================
    HLA Copy Number Workflow
========================================================================================
*/
params.sequenza_solution = 1 // Sequenza solution to evaluate

include { STAGE_INPUT     } from '../process/stage_input'      addParams( params.modules["stage_input"] )
include { HLATYPE_COMPARE } from '../process/hlatype_compare'  addParams( params.modules["hlatype_compare"] )
include { PHASE_HLA       } from '../process/phase_hla'        addParams( params.modules["phase_hla"] )

ch_genome_fasta = file(params.genome_fasta)

include { create_hlacn_channel } from '../lib/functions.nf'
workflow HLACN {
    take:
        ch_input // tuple path(manifest), path(input_folder)
        ch_hisat2_hlatypes // tuple specimen_id, path(hisat2_hlatypes)

    main:

    STAGE_INPUT (
        ch_input
    )
    STAGE_INPUT.out.tsv
      .splitCsv(sep: '\t', header: true)
      .map { create_hlacn_channel(it) }
      .map { it -> 
             [ it[0].specimen_id, it ] 
       }
      .combine(ch_hisat2_hlatypes, by: 0)
      .map { specimen_id, data, hisat2_hlatypes ->
             data + hisat2_hlatypes 
      }
      .set { ch_hlacn_input } // tuple val(meta), val[normal_bam, normal_bai, tumor_bam, tumor_bai, sequenzaModelRData, epic_hlatypes], file(hisat2_hlatypes)

    HLATYPE_COMPARE (
        ch_hlacn_input
    )

    PHASE_HLA (
        HLATYPE_COMPARE.out.hlatypes, // tuple val(meta), val[normal_bam, normal_bai, tumor_bam, tumor_bai, sequenzaModelRData, epic_hlatypes], file(hlatypes)
        ch_genome_fasta
    )
}