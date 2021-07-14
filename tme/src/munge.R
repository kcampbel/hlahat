symbol_check <- function(query, target, ref){
    # Gene symbol checker
    # Useful for replacing geneset symbols with those in the gene expression matrix.
    # Searches ref for query and returns the replacement symbols found in target.
    # Ref is the hgnc database, which must have columns 'approved_symbol', 'alias_symbol', 'previous_symbol'
    # 
    # Args:
    #     query(vector): gene symbols to process
    #     target(vector): gene symbols from expression matrix
    #     ref(dataframe): hgnc database
    #
    # Returns:
    #     vector of replaced gene symbols
    pull_symbol <- function(df, query_col, output_col){
      df %>%
        filter(get(query_col)) %>%
        pull(get(output_col)) %>%
        unique()
    }

    # Filter queries already in target
    not_in_target = query[!query %in% target]

    # Search hgnc database for query genes. Preference is given to approved, alias, and previous symbols in 
    # that order. If the selected gene is not found in the target, return NA.
    map_l = list()
    for(gene in not_in_target){
        df <- hgnc %>% 
          filter(approved_symbol %in% gene | alias_symbol %in% gene | previous_symbol %in% gene) %>%
          mutate(query = gene,
                 approved = approved_symbol %in% target,
                 alias = alias_symbol %in% target,
                 previous = previous_symbol %in% target)
        for(ii in c('approved', 'alias', 'previous')){
            if(any(df[[ii]])){
                hit = pull_symbol(df, ii, paste0(ii, '_symbol'))
                if(length(hit) != 0){
                  break
		}
            } else {
                hit <- NA
	    }
	}
        if(length(hit) != 1){
          cat(paste('WARNING:', gene, 'matches:', paste(hit, collapse=', '), '\n'))
        } else if(is.na(hit)){
          cat(paste('WARNING:', gene, 'gene not found in either ref or target', '\n'))
        } else {
          #cat(paste(gene, '->', hit, '\n'))
          query[which(query == gene)] <- hit
          map_l[[gene]] = hit
        }
    }
    return(query)
}

zscore<- function(x){
  z<- (x - mean(x)) / sd(x)
  return(z)
}

percentile <- function(x, value) ecdf(x)(value)

findoverlaps_metadata <- function(query,subject){
  hits <- GenomicRanges::findOverlaps(query, subject)
  subject_mcols <- GenomicRanges::mcols(subject)
  subj.df = as.data.frame(
    matrix(nrow=length(query), ncol=ncol(subject_mcols), 
           dimnames=list(seq(1, length(query)), names(subject_mcols))
    )
  )
  
  subj.df[S4Vectors::queryHits(hits), ] <- as.data.frame(
    GenomicRanges::mcols(subject[S4Vectors::subjectHits(hits)])
    )
  out <- cbind(GenomicRanges::mcols(query), subj.df)
  return(as.data.frame(out))
}

