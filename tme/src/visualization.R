dt_ext <- c('Buttons', 'ColReorder', 'FixedHeader')
dt_opts <- list(
    dom = 'flrtBip', buttons = I(c('colvis', 'copy', 'csv', 'excel', 'print')), # Buttons
    colReorder = TRUE, # ColReorder
    fixedHeader = TRUE#, # FixedHeader
    #fixedColumns = list(leftColumns = 3), scrollX = TRUE# FixedHeader
    #deferRender = TRUE, scrollY = 200, scroller = TRUE # Scroller
)

link_gene_string <- function(String, source = link_cbioportal, gois = NULL, collapser = ', '){
    ### Add html links to gene names 
    genes <- strsplit(String,'[ ]?,[ ]?')[[1]]
    if(!is.null(gois)){
        genes = genes[genes %in% gois]
    }
    outs <- paste(sapply(
        genes,
        source), 
        collapse = collapser)
    return(outs)
}

link_cbioportal <- function(gene){
    ### Helper function for link_gene_string
    require(kableExtra)
    link <- paste0('https://www.cbioportal.org/results/cancerTypesSummary?Action=Submit&RPPA_SCORE_THRESHOLD=2.0&Z_SCORE_THRESHOLD=2.0&cancer_study_list=laml_tcga_pan_can_atlas_2018%2Cacc_tcga_pan_can_atlas_2018%2Cblca_tcga_pan_can_atlas_2018%2Clgg_tcga_pan_can_atlas_2018%2Cbrca_tcga_pan_can_atlas_2018%2Ccesc_tcga_pan_can_atlas_2018%2Cchol_tcga_pan_can_atlas_2018%2Ccoadread_tcga_pan_can_atlas_2018%2Cdlbc_tcga_pan_can_atlas_2018%2Cesca_tcga_pan_can_atlas_2018%2Cgbm_tcga_pan_can_atlas_2018%2Chnsc_tcga_pan_can_atlas_2018%2Ckich_tcga_pan_can_atlas_2018%2Ckirc_tcga_pan_can_atlas_2018%2Ckirp_tcga_pan_can_atlas_2018%2Clihc_tcga_pan_can_atlas_2018%2Cluad_tcga_pan_can_atlas_2018%2Clusc_tcga_pan_can_atlas_2018%2Cmeso_tcga_pan_can_atlas_2018%2Cov_tcga_pan_can_atlas_2018%2Cpaad_tcga_pan_can_atlas_2018%2Cpcpg_tcga_pan_can_atlas_2018%2Cprad_tcga_pan_can_atlas_2018%2Csarc_tcga_pan_can_atlas_2018%2Cskcm_tcga_pan_can_atlas_2018%2Cstad_tcga_pan_can_atlas_2018%2Ctgct_tcga_pan_can_atlas_2018%2Cthym_tcga_pan_can_atlas_2018%2Cthca_tcga_pan_can_atlas_2018%2Cucs_tcga_pan_can_atlas_2018%2Cucec_tcga_pan_can_atlas_2018%2Cuvm_tcga_pan_can_atlas_2018&case_set_id=all&data_priority=0&gene_list=',
                   gene, '&geneset_list=%20&tab_index=tab_visualize')
    return(text_spec(gene, link = link, tooltip = 'Open gene in cBioPortal'))
}

dt_quantile_ramp <- function(df, cols, caption="None",  gene_name_col = 'Gene'){
    df$Gene <- sapply(df[[gene_name_col]], link_gene_string)
    if(caption == "None"){
        dt <- datatable(df,
                        extensions = dt_ext,
                        options = dt_opts,
                        rownames=FALSE,
                        escape = FALSE
                        )
    } else {
        title <- htmltools::tags$caption(
          caption,
          style="font-size:16px"#; font-weight:bold"
        )
        dt <- datatable(df, 
                        caption=title,
                        extensions = dt_ext,
                        options = dt_opts,
                        rownames=FALSE,
                        escape = FALSE
                        )
    }
    colRamp <- colorRamp(c("lightblue","white","red"))
    for(column in cols){
        x <- na.omit(df[[column]])
        brks <- quantile(x, probs = seq(.05, .95, .005))
        RGB <- colRamp(c(0, (brks-min(x))/(max(x)-min(x))))
        clrs <- apply(RGB, 1, function(rgb){
          sprintf("rgb(%s)", toString(round(rgb,0)))
        })
        dt <- dt %>% 
          formatStyle(column, backgroundColor = styleInterval(brks, clrs))
    }
    return(dt)
}

boxplot_goi <- function(df){
    pal_3grp <- c("red","#113C63","grey50","grey81")
    ax <- list(title = '')
    fig <- plot_ly(df,
                    x = ~gene, y = ~bcTPM, color = ~Source, 
                    type='box', colors=pal_3grp, width=1400) %>%
      layout(boxmode = "group", xaxis=ax)
}

heatmap_annot <- function(sid, emat, subtype_df='None'){
    # Sample annot
    sids <- colnames(emat)
    annCol <- data.frame(
        Sample = factor(if_else(sids == sid, sid, 'Other')),
        row.names=colnames(emat)
    )

    # Subtype annot
    if(subtype_df != 'None'){
        df <- subtype_df[rownames(annCol),] %>%
            dplyr::select(Subtype) 
        annCol <- merge(annCol, df, by="row.names",all.x=TRUE) %>%
            column_to_rownames('Row.names')
        annCol <- annCol[colnames(emat), ]
    }
    return(annCol)
}

