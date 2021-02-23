library(Biobase)
library(purrr)
library(dplyr)
library(stringr)
library(tibble)
# Genes
pattern = '*rsem.*.txt.gz'
in_path = '/home/csmith/git/RNAseqDB/data/expected_count'
files <- list.files(in_path, pattern, full.names=TRUE, recursive=TRUE)
if(length(files) == 0){
  print(paste0(pattern, ' files found'))
}
tissue2tcga <- read.table('rnaseqdb/tissue_map.tsv', sep='\t', header=TRUE)

counts_l <- list()
meta <- data.frame()
for(ii in files){
  print(ii)
  # Expression matrix
  dat <- read.table(ii, sep = '\t', header=1)
  expressed <- rowSums(dat[,3:ncol(dat)]) != 0
  counts_l[[ii]] <- dat[expressed,] %>%
    dplyr::select(-Entrez_Gene_Id)
  sampleId <- colnames(dat[3:ncol(dat)])

  # Metadata
  tmp <- data.frame(tissue=NA, algorithm=NA, count=NA, db=NA)
  fields <- c('tissue', 'algorithm', 'count', 'db')
  tmp[, fields] = str_split(basename(ii), '-|\\.', simplify=TRUE)[1:4]

  if(tmp$db == 'gtex'){
    tmp$tcga <- tissue2tcga %>% filter(tissue == 'bladder') %>% pull(tcga) %>% unique()
  } else {
    tmp$tcga <- tmp$tissue
  }

  if(grepl('*-t.txt.gz', ii)){
    tmp$source <- 'tumor'
  } else {
    tmp$source <- 'normal'
  }

  tmp <- cbind(sampleId, tmp)
  rownames(tmp) <- tmp$sampleId
  tmp <- subset(tmp, select=-sampleId)
  meta <- rbind(meta, tmp)
}

# Merge list
gene_id <- 'Hugo_Symbol'
exprs <- counts_l %>%
  reduce(full_join, by = gene_id) %>%
  replace(is.na(.), 0) %>%
  data.frame(row.names='Hugo_Symbol', check.names=FALSE)

# Eset
phenoData <- new("AnnotatedDataFrame", data = meta)
eset <- Biobase::ExpressionSet(as.matrix(exprs), phenoData)#, featureData)
saveRDS(eset, 'rnaseqdb/rnaseqdb_rsem_expectedcount.RDS')

