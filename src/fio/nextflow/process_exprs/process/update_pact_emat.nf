/*
    Update PACT expression set, log and metadata
*/

include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]
params.pact_emat_s3 = false

process UPDATE_PACT_EMAT {
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
    tuple path(pact_gid_counts), path(pact_gid_tpm), path(pact_hgnc_tpm), path(tcga_gtex_map), path(pact_emat_log), path(blacklist)
    path metadata_tsv 

    output:
    tuple val(meta), path(metadata_tsv_updated), path(pact_hgnc_tpm_updated), path(blacklist_updated), emit: pact_hgnc_tpm_updated
    tuple val(meta), path(metadata_tsv_updated), path(pact_gid_counts_updated), emit: pact_gid_counts_updated
    path pact_gid_tpm_updated
    path pact_emat_log_updated
    path ".command*"

    script:
    outpath                  = "./pact_emat"
    metadata_tsv_updated     = outpath + '/' + metadata_tsv.toString()
    pact_emat_log_updated    = outpath + '/' + pact_emat_log.toString()
    pact_hgnc_tpm_updated    = outpath + '/' + pact_hgnc_tpm.toString()
    pact_gid_tpm_updated     = outpath + '/' + pact_gid_tpm.toString()
    pact_gid_counts_updated  = outpath + '/' + pact_gid_counts.toString()
    blacklist_updated        = outpath + '/' + blacklist.toString()

    if( params.pact_emat_s3 )
        """
        update_exprs.py ${gene_counts}\
        --exprs ${pact_hgnc_tpm}\
        --manifest ${manifest}\
        --exprs_log ${pact_emat_log}\
        --metadata_tsv ${metadata_tsv}\
        --tcga_gtex_map ${tcga_gtex_map}\
        --blacklist ${blacklist}\
        --multiqc ${multiqc}\
        -o ${outpath}

        timestamp_exprs.py\
        --exprs ${pact_hgnc_tpm_updated}\
        --exprs_log ${pact_emat_log_updated}\
        --metadata_tsv ${metadata_tsv_updated}\
        --blacklist ${blacklist_updated}\
        -o ${params.pact_emat_s3}

        update_exprs.py ${gene_counts}\
        --exprs ${pact_gid_tpm}\
        --manifest ${manifest}\
        --exprs_log ${pact_emat_log_updated}\
        -o ${outpath}

        timestamp_exprs.py\
        --exprs ${pact_gid_tpm_updated}\
        --exprs_log ${pact_emat_log_updated}\
        -o ${params.pact_emat_s3}

        update_exprs.py ${gene_counts}\
        --exprs ${pact_gid_counts}\
        --manifest ${manifest}\
        --exprs_log ${pact_emat_log_updated}\
        -o ${outpath}

        timestamp_exprs.py\
        --exprs ${pact_gid_counts_updated}\
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
        --blacklist ${blacklist}\
        --multiqc ${multiqc}\
        -o ${outpath}

        update_exprs.py ${gene_counts}\
        --exprs ${pact_gid_tpm}\
        --manifest ${manifest}\
        --exprs_log ${pact_emat_log_updated}\
        -o ${outpath}

        update_exprs.py ${gene_counts}\
        --exprs ${pact_gid_counts}\
        --manifest ${manifest}\
        --exprs_log ${pact_emat_log_updated}\
        -o ${outpath}

        """
}
