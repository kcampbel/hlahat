// Update PACT RNAseq expression matrix

include { SAMPLE_TO_PACT_EMAT } from '../process/pact_exprs'  addParams( params.modules["update_pact_emat"] )
include { HELLOWORLD       } from '../workflow/helloworld'       addParams( params.modules["pact_emat"] )

workflow UPDATE_PACT_EMAT {
    take:
    ch_input
    ch_pact_emat
    ch_metadata_tsv

    main:
//    HELLOWORLD (
//        ch_input,
//        ch_pact_emat,
//        ch_metadata_tsv
//    )

    SAMPLE_TO_PACT_EMAT (
        ch_input,
        ch_pact_emat,
        ch_metadata_tsv
    )
    
    emit:
    pact_emat_updated = SAMPLE_TO_PACT_EMAT.out.pact_emat_updated
}