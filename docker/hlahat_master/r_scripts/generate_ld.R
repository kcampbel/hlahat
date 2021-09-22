library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(docopt)

##### Read in arguments #####
"Usage:
  generate_ald.R [options]

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

##### Read in files #####
pu_t <- fread(tumor_pileups)
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
  separate(INFO, c("closest_nuc", "mut", "alt_gen", "alt_nuc", "alt_pos"), sep = "\\|") %>%
  separate(mut, c("REF", "ALT"), sep = ">") %>%
  gather(type, nt, REF, ALT) %>% gather(tm, count, A, T, C, G, N, Other) %>%
  filter(nt == tm) %>% gather(tm, value, nt, count) %>%
  unite(label, type, tm, sep = "_") %>% mutate(label = factor(label, levels = c("REF_nt", "ALT_nt", "REF_count", "ALT_count"))) %>%
  spread(label, value) %>% 
  mutate(START = NULL, 
         END = as.numeric(END), alt_pos = as.numeric(alt_pos)) %>%
  mutate(totalMapped = flagmet_t$`in total (QC-passed reads + QC-failed reads)`,
         REF_count.MappedNormalized = as.numeric(REF_count)/(totalMapped/1e6))
snpcomp_t <- snp_t %>%
  inner_join(snp_t, by = c("CHROM" = "alt_gen", "closest_nuc" = "alt_nuc", "END" = "alt_pos", "REF_nt" = "ALT_nt",
                          "alt_gen" = "CHROM", "alt_nuc" = "closest_nuc", "alt_pos" = "END", "ALT_nt" = "REF_nt", "totalMapped"),
            suffix = c(".ref", ".alt"))

##### If Normal provided, normalize to normal values #####
if(normal_pileups != "NA") {
  pu_n <- fread(normal_pileups)
  snp_n <- pu_n %>%
    separate(INFO, c("closest_nuc", "mut", "alt_gen", "alt_nuc", "alt_pos"), sep = "\\|") %>%
    separate(mut, c("REF", "ALT"), sep = ">") %>%
    gather(type, nt, REF, ALT) %>% gather(tm, count, A, T, C, G, N, Other) %>%
    filter(nt == tm) %>% gather(tm, value, nt, count) %>%
    unite(label, type, tm, sep = "_") %>% spread(label, value)
  #
  flag_n <- fread(normal_flagstat, sep = '\t', head = F, col.names = c("line"))
  flagmet_n <- flag_n %>%
    mutate(metric = gsub("^\\d+ \\+ \\d+ (\\w+.*)\\(*.*", "\\1", line),
           metric = gsub(" \\(\\d+.+\\)\\b", "", metric),
           value = as.numeric(gsub("(^\\d+) .+", "\\1", line)),
           line = NULL) %>%
    deframe %>% as.list
  snp_n <- pu_n %>%
    separate(INFO, c("closest_nuc", "mut", "alt_gen", "alt_nuc", "alt_pos"), sep = "\\|") %>%
    separate(mut, c("REF", "ALT"), sep = ">") %>%
    gather(type, nt, REF, ALT) %>% gather(tm, count, A, T, C, G, N, Other) %>%
    filter(nt == tm) %>% gather(tm, value, nt, count) %>%
    unite(label, type, tm, sep = "_") %>%
    mutate(label = factor(label, levels = c("REF_nt", "ALT_nt", "REF_count", "ALT_count"))) %>%
    spread(label, value) %>%
    mutate(START = NULL,
           END = as.numeric(END), alt_pos = as.numeric(alt_pos)) %>%
    mutate(totalMapped_NORMAL = flagmet_n$`in total (QC-passed reads + QC-failed reads)`,
           REF_count.MappedNormalized = as.numeric(REF_count)/(totalMapped_NORMAL/1e6))
  snpcomp_n <- snp_n %>%
    inner_join(snp_n, by = c("CHROM" = "alt_gen", "closest_nuc" = "alt_nuc", "END" = "alt_pos", "REF_nt" = "ALT_nt",
                            "alt_gen" = "CHROM", "alt_nuc" = "closest_nuc", "alt_pos" = "END", "ALT_nt" = "REF_nt",
                            "totalMapped_NORMAL"),
              suffix = c(".ref_NORMAL", ".alt_NORMAL"))
  snpcomp_tn <- snpcomp_t %>% inner_join(snpcomp_n) %>%
    mutate(Gene = gsub("(^\\w+\\d*)\\*\\d+.+", "\\1", CHROM),
           ld = log2(REF_count.MappedNormalized.ref/REF_count.MappedNormalized.ref_NORMAL) - log2(REF_count.MappedNormalized.alt/REF_count.MappedNormalized.alt_NORMAL))
  snp_ratios <- snpcomp_tn %>%
    filter(as.numeric(REF_count.ref) + as.numeric(REF_count.alt) >= tumor_depth_cutoff &
             as.numeric(REF_count.ref_NORMAL) + as.numeric(REF_count.alt_NORMAL) >= normal_depth_cutoff &
             as.numeric(REF_count.ref_NORMAL) > 0 & as.numeric(REF_count.alt_NORMAL)>0) %>%
    group_by(Gene, CHROM, END, alt_gen, alt_pos, ld) %>% select()
  allele_ratios <- snp_ratios %>% group_by(Gene, CHROM, alt_gen) %>%
    summarise(#min_ld = min(ld, na.rm = T), median_ld = median(ld, na.rm = T), max_ld = max(ld, na.rm = T),
      value = quantile(ld, c(0, 0.25, 0.5, 0.75, 1)), metric = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"),
              n_snps = n()) %>%
    mutate(metric = factor(metric, levels = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"))) %>%
    spread(metric, value)
  locus_ratio <- allele_ratios %>% filter(ld_median > 0) %>% ungroup %>%
    summarise(#min_ld = min(median_ld), median_ld = median(median_ld), max_ld = max(median_ld),
      value = quantile(ld_median, c(0, 0.25, 0.5, 0.75, 1)), metric = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"),
              n_snps = sum(n_snps)) %>%
    mutate(metric = factor(metric, levels = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"))) %>%
    spread(metric, value)
} else {
  ##### Otherwise just normalize to tumor only values #####
  snpcomp_tn <- snpcomp_t %>%
    mutate(Gene = gsub("(^\\w+\\d*)\\*\\d+.+", "\\1", CHROM),
           ld = log2(REF_count.MappedNormalized.ref) - log2(REF_count.MappedNormalized.alt))
  snp_ratios <- snpcomp_tn %>%
    filter(as.numeric(REF_count.ref) + as.numeric(REF_count.alt) >= tumor_depth_cutoff) %>%
    group_by(Gene, CHROM, END, alt_gen, alt_pos, ld) %>% select()
  allele_ratios <- snp_ratios %>% group_by(Gene, CHROM, alt_gen) %>%
    summarise(#min_ld = min(ld, na.rm = T), median_ld = median(ld, na.rm = T), max_ld = max(ld, na.rm = T),
      value = quantile(ld, c(0, 0.25, 0.5, 0.75, 1)), metric = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"),
              n_snps = n()) %>%
    mutate(metric = factor(metric, levels = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"))) %>%
    spread(metric, value)
  locus_ratio <- allele_ratios %>% filter(ld_median > 0) %>% ungroup %>%
    summarise(#min_ld = min(median_ld), median_ld = median(median_ld), max_ld = max(median_ld),
      value = quantile(ld_median, c(0, 0.25, 0.5, 0.75, 1)), metric = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"),
              n_snps = sum(n_snps)) %>%
    mutate(metric = factor(metric, levels = c("ld_min", "ld_25", "ld_median", "ld_75", "ld_max"))) %>%
    spread(metric, value)
}

# Table outputs #
write.table(snp_ratios, file = paste0(out_prefix, "_snpLD.tsv"), sep = "\t", row.names = F, quote = F)
write.table(allele_ratios, file = paste0(out_prefix, "_alleleLD.tsv"), sep = "\t", row.names = F, quote = F)
write.table(locus_ratio, file = paste0(out_prefix, "_locusLD.tsv"), sep = "\t", row.names = F, quote = F)







##### END #####

