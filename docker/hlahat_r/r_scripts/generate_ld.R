library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(docopt)

##### Read in arguments #####
"Usage:
  generate_ld.R [options]

Options:
--tumor_depth_cutoff=<Threshold>    [default: 15]
--tumor_pileups=PATH
--tumor_flagstat=PATH
--normal_depth_cutoff=<Threshold>   [default: 15]
--normal_pileups=PATH               [default: NA]
--normal_flagstat=PATH              [default: NA]
--out_prefix=OUTPUT                 [default: out]
" -> doc

opts <- docopt(doc)

tumor_depth_cutoff <- as.numeric(opts$tumor_depth_cutoff)
tumor_pileups <- opts$tumor_pileups
tumor_flagstat <- opts$tumor_flagstat
normal_depth_cutoff <- as.numeric(opts$normal_depth_cutoff)
normal_pileups <- opts$normal_pileups
normal_flagstat <- opts$normal_flagstat
out_prefix <- opts$out_prefix

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

##### Read in files #####
pu_t <- fread(tumor_pileups, col.names = c("a1_gen","START", "a1_pos", "INFO", "REF", "a1_depth", "A", "T", "C", "G", "N", "Other"))
flag_t <- fread(tumor_flagstat, sep = '\t', head = F, col.names = c("line"))

##### Format flagstats #####
flagmet_t <- flag_t %>%
  mutate(metric = gsub("^\\d+ \\+ \\d+ (\\w+.*)\\(*.*", "\\1", line),
         metric = gsub(" \\(\\d+.+\\)\\b", "", metric),
         value = as.numeric(gsub("(^\\d+) .+", "\\1", line)),
         line = NULL) %>%
  deframe %>% as.list

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

