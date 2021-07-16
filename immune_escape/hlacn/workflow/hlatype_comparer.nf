// Stage inputs

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/


// Check input path parameters to see if they exist
//checkPathParamList = [
//    params.input
    // params.fasta, params.transcript_fasta, params.additional_fasta,
//]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
//if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { HLA_COMPARE } from '../process/hlatype_compare'  addParams( params.modules["hlatype_compare"] )

/*
========================================================================================
    RUN WORKFLOW
========================================================================================
*/ 
ch_hlahat_hlatypes = HISAT_GENOTYPE.out.hla_types
//ch_imgt_alignments = file(params.imgt_alignments)
//ch_ref_hla_genotypes = file(params.ref_hla_genotypes)
//ch_ref_hla_bed = file(params.ref_hla_bed)

workflow HLATYPE_COMPARE {
    ch_input_hlacn,
    ch_hlahat_hlatypes
}