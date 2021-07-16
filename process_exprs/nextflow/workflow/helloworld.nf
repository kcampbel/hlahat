
workflow HELLOWORLD {
    take:
    ch_input
    ch_pact_emat
    ch_metadata_tsv

    main:
    helloworld (
    ch_input,
    ch_pact_emat,
    ch_metadata_tsv
    )
}

process helloworld {
    echo true

    input:
    tuple val(meta), path(gene_counts), path(manifest)
    tuple path(pact_gid_tpm), path(pact_hgnc_tpm), path(tcga_gtex_map), path(pact_emat_log)
    //path pact
    path metadata_tsv 
    
    //log.info "$meta, $gene_counts, $manifest"
    script:
    //log.info "$meta"
    //log.info "$pact_gid_tpm"
    """
    """
}