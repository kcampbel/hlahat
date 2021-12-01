library(dplyr)
library(tidyr)
library(data.table)

##### Read in arguments #####
args <- commandArgs(TRUE)
# bed_file <- args[1]
pileups_file <- args[1]
output_prefix <- args[2]

# read_bed <- fread(bed_file, sep = "\t", col.names = c("CHROM", "START", "END", "INFO"))
read_pileups <- fread(pileups_file, col.names = c('CHROM','POS','REF','DEPTH','BASES','BQ_ASCII'))

summBases <- function(test, ref_base) {
  counts <- rep(0, 6)
  names(counts) <- c("A", "T", "C", "G", "N", "Other")
  if(is.na(test)) {
    return(counts)
  } else {
    i = 1
    while (i < nchar(test)+1){
      this_char <- toupper(substring(test, i, i))
      if (this_char %in% c('.', ',')) {
        counts[ref_base] <- counts[ref_base] + 1
        i = i + 1
      } else if (this_char %in% c('A','T','C','G','N')) {
        counts[this_char] <- counts[this_char] + 1
        i = i + 1
      } else if (this_char %in% c("-", "+")){
        n_nt <- as.numeric(substring(test, i+1, i+1))
        # indel <- toupper(substring(test, i+2, i+1+n_nt))
        counts["Other"] <- counts["Other"]+ 1
        i = i + 2 + n_nt
      } else if (this_char == "^") {
        i = i + 2
      } else if (this_char %in% c("$", "*", ">", "<")) {
        i = i + 1
      } else {
        i = i + 1
      }
    }
    return(counts)
  }
}
summ_pileups <- read_pileups %>% group_by(CHROM, POS) %>% 
  mutate(summBase = list(summBases(BASES, REF))) %>% unnest_wider(summBase) %>% mutate(BASES = NULL, BQ_ASCII = NULL)

# final_pileups <- inner_join(read_bed, summ_pileups, by = c("CHROM" = "CHROM", "END" = "POS"))
write.table(summ_pileups, file = paste0(output_prefix, ".summarize_pileups.txt"), sep = "\t", quote = F, row.names = F)