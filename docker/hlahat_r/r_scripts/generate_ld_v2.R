library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(docopt)
# library(rmarkdown)

##### Read in arguments #####
"Usage:
  generate_ld.R [options]

Options:
--hla_find_types=PATH
--tumor_depth_cutoff=<Threshold>    [default: 15]
--tumor_pileups=PATH
--tumor_flagstat=PATH
--tumor_idxstats=PATH
--normal_depth_cutoff=<Threshold>   [default: 15]
--normal_pileups=PATH               [default: NA]
--normal_flagstat=PATH              [default: NA]
--normal_idxstats=PATH              [default: NA]
--hla_snps=PATH
--rna                               [default: FALSE]
--exon_only                         [default: TRUE]
--hla_exonbed=PATH
--tumor_purity=Numeric
--out_prefix=OUTPUT                 [default: out]
" -> doc

opts <- docopt(doc)

hla_find_types <- "~/hlahat_final/working/tm_dev/SKCM-3N-A9WD-NB.find_hlatypes.tsv"#opts$hla_find_types
tumor_depth_cutoff <- 15#as.numeric(opts$tumor_depth_cutoff)
tumor_pileups <- "~/hlahat_final/working/tm_dev/SKCM-3N-A9WD-TM.summarize_pileups.txt"#opts$tumor_pileups
tumor_flagstat <- "~/hlahat_final/working/tm_dev/SKCM-3N-A9WD-TM.flagstat.txt"#opts$tumor_flagstat
tumor_idxstats <- "~/hlahat_final/working/tm_dev/SKCM-3N-A9WD-TM.idxstats.txt"#opts$tumor_idxstats
normal_depth_cutoff <- 15#as.numeric(opts$normal_depth_cutoff)
normal_pileups <- "~/hlahat_final/working/tm_dev/SKCM-3N-A9WD-NB.summarize_pileups.txt"#opts$normal_pileups
normal_flagstat <- "~/hlahat_final/working/tm_dev/SKCM-3N-A9WD-NB.flagstat.txt"#opts$normal_flagstat
normal_idxstats <- "~/hlahat_final/working/tm_dev/SKCM-3N-A9WD-NB.idxstats.txt"#opts$normal_idxstats
hla_snps <- "~/hlahat_final/working/tm_dev/SKCM-3N-A9WD-NB.custom_hla.SNPs.bed"#opts$hla_snps
rna <- FALSE#opts$rna
exon_only <- TRUE#opts$exon_only
hla_exonbed <- "~/hlahat_final/working/tm_dev/SKCM-3N-A9WD-NB.custom_hla.cds_exons.bed"#opts$hla_custom_fasta
tumor_purity <- 100#opts$tumor_purity
out_prefix <- "~/hlahat_final/working/tm_dev/testout"#opts$out_prefix

# Default values
if(is.na(tumor_depth_cutoff)) {
  tumor_depth_cutoff <- 15
}
if(is.na(normal_depth_cutoff)) {
  normal_depth_cutoff <- 15
}
if(is.na(out_prefix)) {
  out_prefix <- "out"
}
if(rna){
  exon_only <- TRUE
}

hla_gene_order <- c("HFE","F","V","G","H","K","A","L","E","C","B","MICA","MICB","DRA","DRB1","DQA1","DQB1","DOB","TAP2","TAP1","DMB","DMA","DOA","DPA1","DPB1","DPB2")

##### Read in files #####
# Reference information
types <- fread(hla_find_types, sep = '\t')
snps <- fread(hla_snps, sep = '\t', col.names = c('allele_ref','st','sp','info','score','strand'))
snps <- snps %>% mutate(Gene = gsub("^(\\w+)\\*\\d+.+", "\\1", allele_ref)) %>%
  separate(info, c("dir", "allele_ref", "allele_alt", "nt_ref", "nt_alt"), sep = "\\||>")
vars_prep <- snps %>%
  group_by(Gene, allele_ref, dir) %>%
  mutate(snp_no = 1:n(), st = NULL, score = NULL, strand = NULL) %>%
  separate(dir, c("comp_ref", "comp_alt"), sep = "v")
vars <- vars_prep %>%
  left_join(vars_prep, by = c("Gene", "snp_no", 
                              "allele_ref" = "allele_alt", "nt_ref" = "nt_alt", "comp_ref" = "comp_alt",
                              "allele_alt" = "allele_ref", "nt_alt" = "nt_ref", "comp_alt" = "comp_ref"),
            suffix = c("_ref", "_alt"))
if(exon_only == TRUE) {
  exons <- fread(hla_exonbed, sep = '\t', col.names = c('allele_ref','st','sp','info','score','strand'))
  exon_pos <- exons %>% group_by(allele_ref) %>% mutate(exon_no = 1:n()) %>% 
    group_by(allele_ref, exon_no, info) %>% summarise(sp_ref = st:sp)
  vars <- vars %>% left_join(exon_pos)
} else {
  vars <- vars %>% mutate(exon_no = NA, info = NA)
}
vars <- vars %>%
  mutate(heterozygous = ifelse((comp_ref == 1 & comp_alt == 2) | (comp_ref == 2 & comp_alt == 1), TRUE, FALSE),
         exon = ifelse(is.na(exon_no), FALSE, TRUE))
#
##### Tumor metrics #####
pu_t <- fread(tumor_pileups, col.names = c("allele_ref", "sp_ref", "nt_ref", "depth", "A", "T", "C", "G", "N", "Other"))
flag_t <- fread(tumor_flagstat, sep = '\t', head = F, col.names = c("line"))
idx_t <- fread(tumor_idxstats, sep = '\t', col.names = c('allele_ref','length','allele_mapped','allele_unmapped'))
# Format flagstat
flagmet_t <- flag_t %>%
  mutate(metric = gsub("^\\d+ \\+ \\d+ (\\w+.*)\\(*.*", "\\1", line),
         metric = gsub(" \\(\\d+.+\\)\\b", "", metric),
         value = as.numeric(gsub("(^\\d+) .+", "\\1", line)),
         line = NULL) %>%
  deframe %>% as.list
# Tumor full summary
tum <- idx_t %>% mutate(Gene = gsub("^(\\w+)\\*\\d+.+", "\\1", allele_ref)) %>%
  mutate(total_mapped = flagmet_t$mapped) %>%
  left_join(pu_t) %>%
  full_join(vars) %>%
  group_by(allele_ref, sp_ref) %>%
  mutate(count_ref = ifelse(is.na(nt_ref), NA, get(nt_ref)),
         count_alt = ifelse(is.na(nt_alt), NA, get(nt_alt))) %>%
  group_by(Gene, heterozygous, snp_no) %>%
  mutate(snp_cov = ifelse(heterozygous == FALSE, NA, sum(count_ref))) %>%
  ungroup()
tum_filt <- tum %>% filter(heterozygous == TRUE & exon %in% c(TRUE, exon_only) & snp_cov >= tumor_depth_cutoff)

##### Normal (if provided) metrics #####
if(!is.na(normal_pileups)){
  pu_n <- fread(normal_pileups, col.names = c("allele_ref", "sp_ref", "nt_ref", "depth", "A", "T", "C", "G", "N", "Other"))
  flag_n <- fread(normal_flagstat, sep = '\t', head = F, col.names = c("line"))
  idx_n <- fread(normal_idxstats, sep = '\t', col.names = c('allele_ref','length','allele_mapped','allele_unmapped'))
  # Format flagstat
  flagmet_n <- flag_n %>%
    mutate(metric = gsub("^\\d+ \\+ \\d+ (\\w+.*)\\(*.*", "\\1", line),
           metric = gsub(" \\(\\d+.+\\)\\b", "", metric),
           value = as.numeric(gsub("(^\\d+) .+", "\\1", line)),
           line = NULL) %>%
    deframe %>% as.list
  # Normal full summary
  norm <- idx_n %>% mutate(Gene = gsub("^(\\w+)\\*\\d+.+", "\\1", allele_ref)) %>%
    mutate(total_mapped = flagmet_n$mapped) %>%
    left_join(pu_n) %>%
    full_join(vars) %>%
    group_by(allele_ref, sp_ref) %>%
    mutate(count_ref = ifelse(is.na(nt_ref), NA, get(nt_ref)),
           count_alt = ifelse(is.na(nt_alt), NA, get(nt_alt))) %>%
    group_by(Gene, heterozygous, snp_no) %>%
    mutate(snp_cov = ifelse(heterozygous == FALSE, NA, sum(count_ref))) %>%
    ungroup()
  norm_filt <- norm %>% filter(heterozygous == TRUE & exon %in% c(TRUE, exon_only) & snp_cov >= normal_depth_cutoff)
  
  # Combine Tumor-Normal metrics
  tumnorm <- tum %>% full_join(norm, by = c("allele_ref", "length", "Gene", "sp_ref", "nt_ref", "comp_ref", "comp_alt", "allele_alt", "nt_alt", "snp_no", "exon_no", "info", "heterozygous", "exon"), suffix = c(".tumor", ".normal"))   
  tumnorm_filts <- tum_filt %>% full_join(norm_filt, by = c("allele_ref", "length", "Gene", "sp_ref", "nt_ref", "comp_ref", "comp_alt", "allele_alt", "nt_alt", "snp_no", "exon_no", "info", "heterozygous", "exon"), suffix = c(".tumor", ".normal"))   
  
  # write.table(tumnorm, file = paste0(out_prefix, ".annotated_pileups.tsv"), sep = "\t", row.names = F, quote = F)
  # Coverage metrics by SNP, allele, and locus
}
# else {
# }

tumnorm_hets %>% 
  group_by(Gene, allele_ref, allele_alt, comp_ref, comp_alt, snp_no, heterozygous) %>%
  

##### Compare SNP sequencing coverage #####
snp_t <- pu_t %>%
  separate(INFO, c("a1_nuc", "mut", "a2_gen", "a2_nuc", "a2_pos"), sep = "\\|") %>%
  separate(mut, c("ref", "alt"), sep = ">") %>%
  gather(type, nt, ref, alt) %>% gather(tm, count, A, T, C, G, N, Other) %>%
  filter(nt == tm) %>% gather(tm, value, nt, count) %>% 
  unite(label, type, tm, sep = "_") %>% mutate(label = factor(label, levels = c("ref_nt", "alt_nt", "ref_count", "alt_count"))) %>% 
  spread(label, value) %>% 
  mutate(START = NULL, REF = NULL, 
         a1_pos = as.numeric(a1_pos), a2_pos = as.numeric(a2_pos),
         ref_count = as.numeric(ref_count), alt_count = as.numeric(alt_count))
snpcomp_t <- snp_t %>%
  inner_join(snp_t, by = c("a1_gen" = "a2_gen", "a1_nuc" = "a2_nuc", "a1_pos" = "a2_pos", "ref_nt" = "alt_nt",
                          "a2_gen" = "a1_gen", "a2_nuc" = "a1_nuc", "a2_pos" = "a1_pos", "alt_nt" = "ref_nt"),
            suffix = c(".a1", ".a2"))

##### If Normal provided, normalize to normal values #####
if(normal_pileups != "NA") {
  pu_n <- fread(normal_pileups, col.names = c("a1_gen","START", "a1_pos", "INFO", "REF", "a1_depth", "A", "T", "C", "G", "N", "Other"))
  #
  flag_n <- fread(normal_flagstat, sep = '\t', head = F, col.names = c("line"))
  flagmet_n <- flag_n %>%
    mutate(metric = gsub("^\\d+ \\+ \\d+ (\\w+.*)\\(*.*", "\\1", line),
           metric = gsub(" \\(\\d+.+\\)\\b", "", metric),
           value = as.numeric(gsub("(^\\d+) .+", "\\1", line)),
           line = NULL) %>%
    deframe %>% as.list
  snp_n <- pu_n %>%
    separate(INFO, c("a1_nuc", "mut", "a2_gen", "a2_nuc", "a2_pos"), sep = "\\|") %>%
    separate(mut, c("ref", "alt"), sep = ">") %>%
    gather(type, nt, ref, alt) %>% gather(tm, count, A, T, C, G, N, Other) %>%
    filter(nt == tm) %>% gather(tm, value, nt, count) %>% 
    unite(label, type, tm, sep = "_") %>% mutate(label = factor(label, levels = c("ref_nt", "alt_nt", "ref_count", "alt_count"))) %>% 
    spread(label, value) %>% 
    mutate(START = NULL, REF = NULL, 
           a1_pos = as.numeric(a1_pos), a2_pos = as.numeric(a2_pos),
           ref_count = as.numeric(ref_count), alt_count = as.numeric(alt_count))
  snpcomp_n <- snp_n %>%
    inner_join(snp_n, by = c("a1_gen" = "a2_gen", "a1_nuc" = "a2_nuc", "a1_pos" = "a2_pos", "ref_nt" = "alt_nt",
                             "a2_gen" = "a1_gen", "a2_nuc" = "a1_nuc", "a2_pos" = "a1_pos", "alt_nt" = "ref_nt"),
               suffix = c(".a1_matchednormal", ".a2_matchednormal"))
  snpcomp_tn <- snpcomp_t %>% inner_join(snpcomp_n) %>%
    mutate(Gene = gsub("(^\\w+\\d*)\\*\\d+.+", "\\1", a1_gen),
           ld = log2(ref_count.a1) - log2(ref_count.a1_matchednormal) - log2(ref_count.a2) + log2(ref_count.a2_matchednormal))
  snp_ratios <- snpcomp_tn %>%
    filter(as.numeric(ref_count.a1) + as.numeric(ref_count.a2) >= tumor_depth_cutoff &
             as.numeric(ref_count.a1_matchednormal) + as.numeric(ref_count.a2_matchednormal) >= normal_depth_cutoff &
             as.numeric(ref_count.a1_matchednormal) > 0 & as.numeric(ref_count.a2_matchednormal)>0)
  allele_ratios <- snp_ratios %>% group_by(Gene, a1_gen, a2_gen) %>%
    summarise(value = quantile(ld, c(0, 0.25, 0.5, 0.75, 1)), 
              metric = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"),
              n_snps = n()) %>%
    mutate(metric = factor(metric, levels = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"))) %>%
    spread(metric, value)
  locus_ratio <- allele_ratios %>% filter(ld_median > 0) %>% ungroup %>%
    summarise(value = quantile(ld_median, c(0, 0.25, 0.5, 0.75, 1)), 
              metric = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"),
              n_snps = sum(n_snps)) %>%
    mutate(metric = factor(metric, levels = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"))) %>%
    spread(metric, value)
} else {
  ##### Otherwise just normalize to tumor only values #####
  snpcomp_tn <- snpcomp_t %>%
    mutate(Gene = gsub("(^\\w+\\d*)\\*\\d+.+", "\\1", a1_gen),
           ld = log2(ref_count.a1) - log2(ref_count.a2))
  snp_ratios <- snpcomp_tn %>%
    filter(as.numeric(ref_count.a1) + as.numeric(ref_count.a2) >= tumor_depth_cutoff)
  allele_ratios <- snp_ratios %>% group_by(Gene, a1_gen, a2_gen) %>%
    summarise(value = quantile(ld, c(0, 0.25, 0.5, 0.75, 1)),
              metric = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"),
              n_snps = n()) %>%
    mutate(metric = factor(metric, levels = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"))) %>%
    spread(metric, value)
  locus_ratio <- allele_ratios %>% filter(ld_median > 0) %>% ungroup %>%
    summarise(value = quantile(ld_median, c(0, 0.25, 0.5, 0.75, 1)),
              metric = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"),
              n_snps = sum(n_snps)) %>%
    mutate(metric = factor(metric, levels = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"))) %>%
    spread(metric, value)
}

# Table outputs #
write.table(snpcomp_tn, file = paste0(out_prefix, "_allSnps.tsv"), sep = "\t", row.names = F, quote = F)
write.table(snp_ratios, file = paste0(out_prefix, "_snpLD.tsv"), sep = "\t", row.names = F, quote = F)
write.table(allele_ratios, file = paste0(out_prefix, "_alleleLD.tsv"), sep = "\t", row.names = F, quote = F)
write.table(locus_ratio, file = paste0(out_prefix, "_locusLD.tsv"), sep = "\t", row.names = F, quote = F)







##### END #####

