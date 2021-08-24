#!/usr/bin/env Rscript
#options(echo=TRUE)
message(Sys.time(), "::Starting fitSequenza")
library(optparse)
option_list = list(
  make_option(c("-i", "--id"), type="character", default=NULL,
              help="Patient id", metavar="character"),
  make_option(c("-a", "--altMode"), type="character", default=NULL,
              help="Run alt solutions", metavar="character"),
  make_option(c("-o", "--outDir"), type="character", default=NULL,
              help="Output directory", metavar="character"),
  make_option(c("-r", "--sequenzaModelRData"), type="character", default=NULL,
              help="Path to sequenzaModel.RData", metavar="character"),
  make_option(c("-v", "--varsFile"), type="character", default=NULL,
              help="Path variants tsv", metavar="character"),
  make_option(c("-t", "--sequenzaTools"), type="character", default=NULL,
              help="Path to sequenzaTools.R", metavar="character")#,
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# strip trailing / from all parameters in opt
for(i in 1:length(opt)){
  if(substr(opt[i], nchar(opt[i]), nchar(opt[i])) == '/'){
    opt[i] <- substr(opt[i], 1, nchar(opt[i])-1)
  }
}
message(print(opt))
altMode                 <- opt$altMode
id                      <- opt$id
outDir                  <- opt$outDir
sequenzaModelRData      <- opt$sequenzaModelRData
varsFile                <- opt$varsFile
sequenzaTools           <- opt$sequenzaTools

if(length(names(opt)) == 1){
        print_help(opt_parser)
        stop(call.=FALSE)
}

option_names <- unlist(lapply(option_list, function(x) x@dest))
for(ii in option_names){
    if(is.null(opt[[ii]])){
        print_help(opt_parser)
        stop(paste("Missing argument", ii, ".\n", call.=FALSE))
    }
}

# load R libraries
library(sequenza)
library(GenomicRanges)
source(sequenzaTools)
library(sessioninfo)

# load sequenza data model
message(Sys.time(), "::Loading ", sequenzaModelRData)
load(sequenzaModelRData)
seqz_root <- dirname(sequenzaModelRData)

message(Sys.time(), "::Loading ", varsFile)
vars = read.delim(file=varsFile, header=TRUE, sep="\t")

nsols <- nrow(allSols)
if(altMode == 'y'){
    solutions <- nsols
} else {
    solutions <- 1
}

for(s in 1:solutions){
    cellularity = allSols$cellularity[s]
    ploidy = allSols$ploidy[s]
    message(Sys.time(), "::processing solution ", paste(s, 'of', nsols), 
        " cellularity: ", cellularity, " ploidy: ", ploidy)
    sequenzaExtractRData <- file.path(
        seqz_root, 
        paste0('alt', s),
        paste0(id, '_sequenza_extract.RData')
        )
    if(!file.exists(sequenzaExtractRData)){
        message(Sys.time(), "::", sequenzaExtractRData, ' does not exist. Skipping.')
        next
    }
    outDirAlt = file.path(outDir, paste0("alt", s))
    if (!dir.exists(outDirAlt)) {
        dir.create(outDirAlt, recursive=TRUE, showWarnings=FALSE)
    }
    alleleTsvAlt = file.path(outDirAlt, paste(id, "_alleles.tsv", sep=""))
    cn <- mutSeg_cn(id, vars, sequenzaExtractRData, cellularity, ploidy, snp=TRUE, vaf_col='obsVafNorm_flip', method='segDR')
    message(Sys.time(), "::Writing ", alleleTsvAlt)
    write.table(cn, file=alleleTsvAlt, quote=FALSE, sep="\t", row.names=FALSE)
}
session_info()

message(Sys.time(), "::finished R session.")
