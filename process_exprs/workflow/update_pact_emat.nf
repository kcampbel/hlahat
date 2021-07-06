// Update PACT RNAseq expression matrix

include { SAMPLE_TO_PACT_EMAT } from '../process/pact_exprs'  addParams( params.modules["update_pact_emat"] )
//include { TIMESTAMP_EMAT } from '../process/timestamp_emat'  addParams( params.modules["update_pact_emat"] )
include { HELLOWORLD       } from '../workflow/helloworld'       addParams( params.modules["pact_emat"] )

workflow UPDATE_PACT_EMAT {
    take:
    ch_input
    ch_pact_emat
    ch_metadata_tsv

    main:
    HELLOWORLD (
        ch_input,
        ch_pact_emat,
        ch_metadata_tsv
    )

    SAMPLE_TO_PACT_EMAT (
        ch_input,
        ch_pact_emat,
        ch_metadata_tsv
    )
    
//    update_timestamp.py\
//     --exprs ${pact_hgnc_tpm_updated}\
//     --exprs_log ${pact_emat_log_updated}\
//     --metadata_tsv ${metadata_tsv_updated}\
//     -o ${pact_emat_s3

//    TIMESTAMP_EMAT (
//        SAMPLE_TO_PACT_EMAT.out.pact_hgnc_updated //tuple val(meta), path(meta_tsv_updated), path(pact_emat_log_updated, path(pact_hgnc_tpm_updated)
//        outdir = params.pact_emat_outdir
//    )

//    TIMESTAMP_EMAT (
//        SAMPLE_TO_PACT_EMAT.out.pact_gid_updated tuple val(meta), path(meta_tsv_updated), path(pact_emat_log_updated),  path(pact_gid_tpm_updated)
//        outdir = params.pact_emat_outdir
//    )

   // emit:
    //pact_emat = SAMPLE_TO_PACT_EMAT.out.pact_emat_updated
}