// RNAseq batch correction 

params.pact_xena = [:]
/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { PACT_XENA } from '../process/pact_xena'  addParams( params.modules["pact_xena"] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/ 

workflow BATCH_CORRECT {
    take:
    ch_input
    ch_pact_emat
    ch_batch_correct

    main:
    PACT_XENA (
        ch_input,
        ch_pact_emat,
        ch_batch_correct
    )

    emit:
    bc_hgnc_tpm = PACT_XENA.out.bc_emat_updated // tuple val(meta), path(bc_log_updated), path(bc_tpm_updated)
}