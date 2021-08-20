/*
    Batch correction 
*/

include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]
params.bc_emat_s3 = false

process BATCH_CORRECT {
 //   label 'process_medium'
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }
    conda params.conda_basedir + params.conda_envt

    input:
    tuple val(meta), path(manifest), val(primary_site), path(gene_counts), path(multiqc)
    tuple val(meta), path(metadata_tsv), path(pact_hgnc_tpm), path(blacklist)
    tuple path(xena_hgnc_tpm), path(bc_emat_log)
        
    output:
    tuple val(meta), path(bc_log_updated), path(bc_tpm_updated), emit: bc_hgnc_tpm_updated
    path ".command*"

    script:
    outpath = './bc_emat'
    bc_tpm_updated = outpath + '/' + 'PactTcgaGtex_' + primary_site + '_hgnc_bctpm.tsv.gz'
    bc_log_updated = outpath + '/' + bc_emat_log.toString()

    if( params.bc_emat_s3 )
        """
        batch_correct.py\
        ${xena_hgnc_tpm} ${pact_hgnc_tpm}\
        --batch_name study\
        --primary_site ${primary_site}\
        --metadata_tsv ${metadata_tsv}\
        --exprs_log ${bc_emat_log}\
        -o ${bc_tpm_updated}
        
        timestamp_exprs.py\
        --exprs ${bc_tpm_updated}\
        --exprs_log ${bc_log_updated}\
        -o ${params.bc_emat_s3}
        """
    else
        """
        batch_correct.py\
        ${xena_hgnc_tpm} ${pact_hgnc_tpm}\
        --batch_name study\
        --primary_site ${primary_site}\
        --metadata_tsv ${metadata_tsv}\
        --exprs_log ${bc_emat_log}\
        -o ${bc_tpm_updated}
        """
}
