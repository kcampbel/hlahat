mutation_cn <- function(mutTable, cellularity, ploidy, avg.depth.ratio, vaf = 'obsVaf') {
    # Mutation copy number caller
    # Args:
    #     mutTable(df): mutations dataframe
    #     cellularity(float): cellularity of the solution
    #     ploidy(float): ploidy of the solution
    #     avg.depth.ratio(float): average depth ratio of solution
    #     vaf(str): mutTable VAF column name
    # 
    # Returns:
    #     df: allelic copy number calls from model fit
    calAllele <- function(mufreq, depth.ratio) {
        allele = mufreq.bayes(mufreq = mufreq,
                           depth.ratio = depth.ratio,
                           cellularity = cellularity, ploidy = ploidy,
                           avg.depth.ratio = avg.depth.ratio)
        return(allele)
    }

    alleles = matrix(NA, ncol=8, nrow=dim(mutTable)[1])
    colnames(alleles) = c("CNn", "CNt", "Mt",  "LPP", "genotype", "adjVaf", "cellularity", "ploidy")
    for (i in 1:dim(mutTable)[1]) {
        depth.ratio = mutTable[i,]$depth.ratio.seg
        mutfreq = mutTable[i, vaf]
        allele = calAllele(mutfreq, depth.ratio)
        if (dim(allele)[1] > 0) {
            alleles[i,1:4] = unlist(allele)
            genotype = paste(c(rep("A", allele$CNt-allele$Mt), rep("B", allele$Mt)), collapse="")
            alleles[i,5] = genotype
            CNt = allele$CNt
            obsVaf = mutfreq 
            adjVaf = obsVaf*(1+2*(1-cellularity)/(cellularity*CNt)) 
            alleles[i,6] = adjVaf
        }
        alleles[i,7] = cellularity
        alleles[i,8] = ploidy
    }
    sumFile = cbind(mutTable, alleles)

    return(sumFile)
}

segment_cn <- function(seg.tab, avg.depth.ratio, cellularity, ploidy, female,
    CNt.min=1, CNt.max=20){   
    # Segment copy number caller
    # from modified from sequenza.extract
    # Args:
    #    seg.tab{df}: segments data frame
    #    avg.depth.ratio(float): average depth ratio of solution
    #    cellularity(float): cellularity of the solution
    #    ploidy(float): ploidy of the solution
    #    female(bool): female?
    #    CNt.max(int): Maximum tumor copy number
    #    CNt.min(int): Min tumor copy number
    #
    # Returns:
    #     df: segment copy number
    seg.len <- (seg.tab$end.pos - seg.tab$start.pos) / 1e6
    weight.ratio <- seg.tab$start.pos - seg.tab$end.pos

    XY = c(X="chrX", Y="chrY")
    if (female){
        segs.is.xy <- seg.tab$chromosome == XY["Y"]
    } else{
        segs.is.xy <- seg.tab$chromosome %in% XY
    }
    cn.alleles  <- baf.bayes(Bf = seg.tab$Bf[!segs.is.xy], 
        CNt.min = CNt.min, CNt.max = CNt.max,
        depth.ratio = seg.tab$depth.ratio[!segs.is.xy],
        cellularity = cellularity, ploidy = ploidy,
        avg.depth.ratio = avg.depth.ratio,
        sd.ratio = seg.tab$sd.ratio[!segs.is.xy],
        weight.ratio = seg.len[!segs.is.xy],
        sd.Bf = seg.tab$sd.BAF[!segs.is.xy],
        weight.Bf = 1, ratio.priority = FALSE, CNn = 2)
    seg.res <- cbind(seg.tab[!segs.is.xy, ], cn.alleles)

    if (!female){
        if (sum(segs.is.xy) >= 1) {
            cn.alleles  <- baf.bayes(Bf = NA, CNt.max = CNt.max,
                depth.ratio = seg.tab$depth.ratio[segs.is.xy],
                cellularity = cellularity, ploidy = ploidy,
                avg.depth.ratio = avg.depth.ratio,
                sd.ratio = seg.tab$sd.ratio[segs.is.xy],
                weight.ratio = seg.len[segs.is.xy], sd.Bf = NA,
                weight.Bf = NA, ratio.priority = FALSE, CNn = 1)
            seg.xy <- cbind(seg.tab[segs.is.xy, ], cn.alleles)
            seg.res <- rbind(seg.res, seg.xy)
        }
    }
    return(seg.res)
}


mutSeg_cn <- function(id, vars, sequenzaExtractRData, cellularity, ploidy, 
    vaf_col = 'obsVaf', snp = FALSE, method = 'segDR'){
    # Copy number caller for mutations and their segments
    # 
    # Args:
    #     id(str): Sample ID
    #     vars(str): mutation dataframe to fit to sequenza model
    #     sequenzaExtractRData(R): R object
    #     cellularity(float): cellularity of the solution
    #     ploidy(float): ploidy of the solution
    #     vaf_col(str): vars VAF column name
    #     snp(bool): TRUE = process SNPs, FALSE = process SNVs
    #     method(str): if SNP=TRUE, 'segDR' for using segment depth ratio, 'segCNt'to use segment 
    #                  CNt
    # Returns:
    #     df: copy number calls for variants from varsFile and their assigned segments 
    load(sequenzaExtractRData)
    seqz.extract <- get(paste0(id, "_sequenza_extract"))
    seg.tab <- do.call(rbind, seqz.extract$segments)
    avg.depth.ratio <- seqz.extract$avg.depth.ratio

    # Calculate segment copy number
    # Female = TRUE to retain segments on the X. Males have no X mutations in varsFile, 
    # so these segments will be dropped in distancetoNearest below.
    segCn <- segment_cn(seg.tab, avg.depth.ratio, cellularity, ploidy, female=TRUE) 
    segmentRanges = makeGRangesFromDataFrame(segCn, seqnames.field='chromosome', 
        start.field="start.pos", end.field="end.pos", keep.extra.columns=TRUE)

    # read mutations and find nearest segment
    seqnames <- grep('chr|chromosome', colnames(vars), value=TRUE)
    mutRanges = makeGRangesFromDataFrame(vars, seqnames.field=seqnames, 
        start.field="start", end.field="end", keep.extra.columns=TRUE)

    maxDis2Seg = 10000
    hits = distanceToNearest(mutRanges, segmentRanges)
    hits = hits[which(mcols(hits)$distance < maxDis2Seg)]

    # Merge segment and mutation data
    seg.mcols <- mcols(segmentRanges)
    seg.colnames <- sapply(names(seg.mcols), function(x) paste0(x, '.seg'))
    seg.df = as.data.frame(
        matrix(nrow=length(mutRanges), ncol=ncol(seg.mcols), 
        dimnames=list(seq(1, length(mutRanges)), seg.colnames)
        )
    )
    seg.df[queryHits(hits), ] <- as.data.frame(mcols(segmentRanges[subjectHits(hits)]))
    mutTable = cbind(mcols(mutRanges), seg.df)
     
    # Call mutation copy number 
    if(snp) {
        if(method == 'segCNt'){
            print('Using segCNt and snp depth.ratio')
        }
        if(method == 'segDR'){
            print('Using segment depth.ratio')
        }
        snp.cn = data.frame()
        for(ii in 1:nrow(mutTable)){
            row = mutTable[ii, ]
            if(method == 'segCNt'){
                tmp = mufreq.bayes(mufreq = row[[vaf_col]],
                                    depth.ratio = row$depth.ratio,
                                    cellularity = cellularity, ploidy = ploidy,
                                    avg.depth.ratio = avg.depth.ratio, 
                                    CNt.min = row$CNt.seg, CNt.max = row$CNt.seg)
            } 
            if(method == 'segDR') {
                tmp = mufreq.bayes(mufreq = row[[vaf_col]],
                                    depth.ratio = row$depth.ratio.seg,
                                    cellularity = cellularity, ploidy = ploidy,
                                    avg.depth.ratio = avg.depth.ratio,
                                    CNt.min = 1, CNt.max = 7
                                    )
            }
            snp.cn = rbind(snp.cn, tmp)
        }
        allele.df = cbind(mutTable, snp.cn)
    } else {
        allele.df <- mutation_cn(mutTable, cellularity, ploidy, avg.depth.ratio, vaf_col)
    }
    cols <- colnames(allele.df) %in% seg.colnames
    out <- cbind(allele.df[!cols], allele.df[cols])
    return(out)
}
