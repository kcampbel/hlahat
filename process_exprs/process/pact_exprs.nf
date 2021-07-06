/*
    Update PACT expression set, log and metadata
*/

include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]
params.pact_emat_s3 = false

process SAMPLE_TO_PACT_EMAT {
 //   label 'process_medium'
    tag "${meta}"
//    publishDir "${params.outdir}",
//        mode: params.publish_dir_mode,
//        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'merged_fastq', meta:meta, publish_by_meta:['id']) }
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }

    input:
    tuple val(meta), path(gene_counts), path(manifest)
    tuple path(pact_gid_tpm), path(pact_hgnc_tpm), path(tcga_gtex_map), path(pact_emat_log)
    path metadata_tsv 
//    tuple val(meta), path(gene_counts), path(manifest)
//    tuple path(pact_gid_tpm), path(pact_hgnc_tpm), path(tcga_gtex_map), path(pact_emat_log)
//    path metadata_tsv 

    output:
    tuple val(meta), path(metadata_tsv_updated), path(pact_emat_log_updated), path(pact_hgnc_tpm_updated), emit: pact_hgnc_emat_updated
    //tuple val(meta), path(metadata_tsv_updated), path(pact_emat_log_updated), path(pact_gid_tpm_updated), emit: pact_gid_emat_updated
    path outdir

    script:
    outdir                = "./pact_emat_updated"
    metadata_tsv_updated  = outdir + '/' + metadata_tsv
    pact_emat_log_updated = outdir + '/' + pact_emat_log
    pact_hgnc_tpm_updated = outdir + '/' + pact_hgnc_tpm
    pact_gid_tpm_updated  = outdir + '/' + pact_gid_tpm
    //log.info "$metadata_tsv_updated"
    if( params.pact_emat_s3 )
        """
        update_exprs.py ${gene_counts}\
        --exprs ${pact_hgnc_tpm}\
        --manifest ${manifest}\
        --exprs_log ${pact_emat_log}\
        --metadata_tsv ${metadata_tsv}\
        --tcga_gtex_map ${tcga_gtex_map}\
        -o ${outdir}

        timestamp_exprs.py\
        --exprs ${pact_hgnc_tpm_updated}\
        --exprs_log ${pact_emat_log_updated}\
        --metadata_tsv ${metadata_tsv_updated}\
        -o ${params.pact_emat_s3}

        update_exprs.py ${gene_counts}\
        --exprs ${pact_gid_tpm}\
        --manifest ${manifest}\
        --exprs_log ${pact_emat_log}\
        -o ${outdir}

        timestamp_exprs.py\
        --exprs ${pact_gid_tpm_updated}\
        --exprs_log ${pact_emat_log_updated}\
        -o ${params.pact_emat_s3}
        """
    else
        """
        update_exprs.py ${gene_counts}\
        --exprs ${pact_hgnc_tpm}\
        --manifest ${manifest}\
        --exprs_log ${pact_emat_log}\
        --metadata_tsv ${metadata_tsv}\
        --tcga_gtex_map ${tcga_gtex_map}\
        -o ${outdir}

        update_exprs.py ${gene_counts}\
        --exprs ${pact_gid_tpm}\
        --manifest ${manifest}\
        --exprs_log ${pact_emat_log}\
        -o ${outdir}
        """
}
