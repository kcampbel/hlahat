library(seqinr)
library(Biostrings)
library(biomaRt)
library(tidyverse)
library(msa)


# Reference genotypes #
load("~/Desktop/hs37d5_ref_HLA.Rda")
ref_types <- data.frame(type = unlist(lapply(ref_haplotypes, function(i) { i$allele })),
                        id = unlist(lapply(ref_haplotypes, function(i) { i$id }))) %>%
  mutate(gene_alias = gsub("(^\\w+\\d*)\\*\\d+:.+", "\\1", type)) %>%
  mutate(gene = case_when(grepl("^HFE|MIC|TAP", gene_alias) == FALSE ~ paste0("HLA-", gene_alias), TRUE ~ gene_alias))

##### Genomic reference #####
# Load reference sequence information #
# ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl", version = "Ensembl Genes 104")
load("~/Desktop/ensembl104.Rda")
# Genes
find_gene <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
                            filters = c("hgnc_symbol", "chromosome_name"),
                            values = list(hgnc_symbol = ref_types$gene, chromosome_name = "6"),
                            mart = ensembl)
# Transcripts, since ensembl gives transcript sequences
find_transcript <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", 
                                                 "transcript_start", "transcript_end", "strand",
                                                 # 'transcript_mane_select','ccds',
                                                 # 'transcript_appris','transcript_tsl',
                                                 'transcript_is_canonical'),
                                  filters = c("ensembl_gene_id", "chromosome_name"),
                                  values = list(ensembl_gene_id = find_gene$ensembl_gene_id, chromosome_name = "6"),
                                  mart = ensembl) %>% 
  filter(transcript_is_canonical == 1)

# Transcript flank starts at beginning of transcript (including UTR)
upstream_transcript_flank <- biomaRt::getSequence(id = find_transcript$ensembl_transcript_id,
                                                  type = "ensembl_transcript_id",
                                                  seqType = "transcript_flank",
                                                  upstream = 5000,
                                                  mart = ensembl)
colnames(upstream_transcript_flank)[which(colnames(upstream_transcript_flank) == "transcript_flank")] <- "upstream_transcript_flank"
downstream_transcript_flank <- biomaRt::getSequence(id = find_transcript$ensembl_transcript_id,
                                                    type = "ensembl_transcript_id",
                                                    seqType = "transcript_flank",
                                                    downstream = 5000,
                                                    mart = ensembl)
colnames(downstream_transcript_flank)[which(colnames(downstream_transcript_flank) == "transcript_flank")] <- "downstream_transcript_flank"
flank_info <- find_gene %>% full_join(find_transcript) %>% 
  full_join(upstream_transcript_flank) %>% 
  full_join(downstream_transcript_flank)

# Transcript sequences +/- flanks
transcript_seq <- biomaRt::getSequence(id = find_transcript$ensembl_transcript_id,
                                       type = "ensembl_transcript_id",
                                       seqType = "transcript_exon_intron",
                                       mart = ensembl)
transcript_seq <- find_gene %>% full_join(find_transcript) %>% 
  full_join(transcript_seq) %>% mutate(length = nchar(transcript_exon_intron))
transcripts <- left_join(transcript_seq, flank_info) %>%
  mutate(context_seq = paste0(upstream_transcript_flank, transcript_exon_intron, downstream_transcript_flank))

##### Get genomic sequence of reference HLA types #####
# Get Genomic MSF for all genes
gens <- paste0("~/hisat_test/tm/hisat2-hisat2_v2.2.0_beta/hisatgenotype_db/HLA/msf/", ref_types$gene_alias, "_gen.msf")
gen_msf <- lapply(gens, function(msf) { read.alignment(msf, format = "msf") })
names(gen_msf) <- ref_types$gene_alias

ref_seqs <- apply(ref_types, 1, function(i){
  if( any(i['type'] %in% gen_msf[[ i['gene_alias'] ]]$nam) ) {
    gen_used <- i['type']
    gen_seq <- gsub("-", "", toupper(gen_msf[[ i['gene_alias'] ]]$seq[[ which(gen_msf[[ i['gene_alias'] ]]$nam == gen_used) ]]))
  } else if (any(gsub(":\\d+$", "", i['type']) %in% gen_msf[[ i['gene_alias'] ]]$nam)) {
    gen_used <- gsub(":\\d+$", "", i['type'])
    gen_seq <- gsub("-", "", toupper(gen_msf[[ i['gene_alias'] ]]$seq[[ which(gen_msf[[ i['gene_alias'] ]]$nam == gen_used) ]]))
  } else {
    gen_used <- NULL
    gen_seq <- NULL
  }
  data.frame(gen_used, gen_seq)
}) %>% enframe %>%
  mutate(gene_alias = ref_types$gene_alias) %>% full_join(ref_types) %>% unnest(value)

##### Extract positions of genomic reference sequences of reference HLA types from transcript sequences #####
ref_seqsXtranscripts <- ref_seqs %>%
  left_join(transcripts, by = c("gene" = "hgnc_symbol"))
find.gen_seq <- apply(ref_seqsXtranscripts, 1, function(i){
  pa <- pairwiseAlignment(pattern = i['gen_seq'], subject = i['context_seq'], type = "local")
  subject_positions <- data.frame(pa@subject@range)
  return(subject_positions)
})

ref_locsXtranscript_info <- find.gen_seq %>% enframe  %>%
  mutate(gene_alias = ref_types$gene_alias) %>% full_join(ref_types) %>% unnest(value) %>%
  left_join(find_gene, by = c("gene" = "hgnc_symbol")) %>%
  left_join(find_transcript)

ref_genomicpos <- ref_locsXtranscript_info %>%
  # Make 0-based for bed file
  mutate(ref_start = ifelse(strand == 1, transcript_start - 5000 -1 + start - 1, transcript_end + 5000 - end),
         ref_end = ifelse(strand == 1, transcript_start - 5000 + end - 1, transcript_end + 5000 - start + 1))
# ref_genomicpos %>% mutate(chr = "chr6", score = 0, strand = ifelse(strand == 1, "+", "-")) %>%
#   select(chr, ref_start, ref_end, type, score, strand) %>%
#   write.table(., file = "~/hlahat_final/hlahat/docker/hlahat_r/r_scripts/reftype_gDNA.bed", row.names = F, quote = F, sep = '\t', col.names = F)

##### Get genomic coordinates of exons ####
nucs <- paste0("~/hisat_test/tm/hisat2-hisat2_v2.2.0_beta/hisatgenotype_db/HLA/msf/", ref_types$gene_alias, "_nuc.msf")
nuc_msf <- lapply(nucs, function(msf) { read.alignment(msf, format = "msf") })
names(nuc_msf) <- ref_types$gene_alias
find_cds <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", "strand",
                                            "coding"),
                             filters = c("ensembl_transcript_id", "chromosome_name"),
                             values = list(ensembl_transcript_id = find_transcript$ensembl_transcript_id, chromosome_name = "6"),
                             mart = ensembl)
find_exons <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "strand",
                                                 "ensembl_exon_id",
                                                 "exon_chrom_start", "exon_chrom_end", "gene_exon",
                                            "genomic_coding_start", "genomic_coding_end"),
                                  filters = c("ensembl_transcript_id", "chromosome_name"),
                                  values = list(ensembl_transcript_id = find_transcript$ensembl_transcript_id, chromosome_name = "6"),
                                  mart = ensembl)
exons_getcds <- find_exons %>%
  group_by(ensembl_exon_id) %>%
  mutate(gene_exon = ifelse(strand == 1, gene_exon,
                            as.character(complement(reverseComplement(DNAString(gene_exon)))))) %>% ungroup %>%
  mutate(cds_sequence = ifelse(genomic_coding_start == "", NA, 
                               substr(gene_exon, 
                                      as.numeric(genomic_coding_start) - as.numeric(exon_chrom_start) + 1,
                                      as.numeric(genomic_coding_end) - as.numeric(exon_chrom_start) + 1)),
         cds_sequence_len = ifelse(genomic_coding_start == "", NA, nchar(cds_sequence)),
         gene_exon_len = nchar(gene_exon),
         genomic_coding_len = ifelse(genomic_coding_start == "", NA, as.numeric(genomic_coding_end) - as.numeric(genomic_coding_start) + 1))
exoncds <- exons_getcds %>% 
  arrange(exon_chrom_start) %>% filter(!is.na(cds_sequence)) %>%
  group_by(ensembl_gene_id, ensembl_transcript_id, strand) %>%
  summarise(exoncds = ifelse(as.numeric(strand) == 1, 
                         paste0(cds_sequence, collapse = ""),
                         as.character(complement(reverseComplement(DNAString(paste0(cds_sequence, collapse = ""))))))) %>%
  unique() %>% 
  mutate(strand = as.integer(strand))
exonXcoding <- exoncds  %>%
  left_join(find_cds) %>%
  mutate(status = exoncds == coding)
# 5 genes missing: "HLA-DPB2" "HLA-H"    "HLA-K"    "HLA-L"    "HLA-V" 
missing_codings <- setdiff(find_gene$hgnc_symbol, exonXcoding$external_gene_name)

#
ref_nuc <- apply(ref_types, 1, function(i){
  if( any(i['type'] %in% nuc_msf[[ i['gene_alias'] ]]$nam) ) {
    gen_used <- i['type']
    nuc_seq <- gsub("-", "", toupper(nuc_msf[[ i['gene_alias'] ]]$seq[[ which(nuc_msf[[ i['gene_alias'] ]]$nam == gen_used) ]]))
  } else if (any(gsub(":\\d+$", "", i['type']) %in% nuc_msf[[ i['gene_alias'] ]]$nam)) {
    gen_used <- gsub(":\\d+$", "", i['type'])
    nuc_seq <- gsub("-", "", toupper(nuc_msf[[ i['gene_alias'] ]]$seq[[ which(nuc_msf[[ i['gene_alias'] ]]$nam == gen_used) ]]))
  } else {
    gen_used <- NULL
    nuc_seq <- NULL
  }
  data.frame(gen_used, nuc_seq)
}) %>% enframe %>%
  mutate(gene_alias = ref_types$gene_alias) %>% full_join(ref_types) %>% unnest(value)

# Some ref_seqs don't match coding/cds
# [1] "HLA-F"    "HLA-V"    "HLA-H"    "HLA-K"    "HLA-L"    "MICA"     "HLA-DOB"  "TAP2"     "HLA-DPB2"

wrong_tx <- ref_nuc %>% 
  left_join(ref_seqs) %>%
  left_join(find_cds, by = c("gene" = "external_gene_name")) %>%
  left_join(exoncds) %>%
  mutate(status_nucXcoding = nuc_seq == coding,
         status_nucXexoncds = nuc_seq == exoncds,
         len_nuc = nchar(nuc_seq),
         len_gen = nchar(gen_seq),
         len_coding = nchar(coding)) %>% 
  left_join(find_transcript) %>%
  filter(status_nucXcoding == FALSE)

##### Find other transcripts that match refseq ####
other_transcripts <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id",
                                                   "transcript_start", "transcript_end", "strand",
                                          "coding"),
                           filters = c("external_gene_name", "chromosome_name"),
                           values = list(external_gene_name = wrong_tx$gene, chromosome_name = "6"),
                           mart = ensembl)
other_transcripts %>% left_join(ref_nuc, by = c("external_gene_name" = "gene")) %>%
  mutate(status_codingXnuc = coding == nuc_seq) %>%
  filter(status_codingXnuc == TRUE)
# Use HLA-F ENST00000334668
f_exons <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "strand",
                                            "ensembl_exon_id",
                                            "exon_chrom_start", "exon_chrom_end", "gene_exon",
                                            "genomic_coding_start", "genomic_coding_end"),
                             filters = c("ensembl_transcript_id", "chromosome_name"),
                             values = list(ensembl_transcript_id = "ENST00000334668", chromosome_name = "6"),
                             mart = ensembl)
f_exons_getcds <- f_exons %>%
  group_by(ensembl_exon_id) %>%
  mutate(gene_exon = ifelse(strand == 1, gene_exon,
                            as.character(complement(reverseComplement(DNAString(gene_exon)))))) %>% ungroup %>%
  mutate(cds_sequence = ifelse(genomic_coding_start == "", NA, 
                               substr(gene_exon, 
                                      as.numeric(genomic_coding_start) - as.numeric(exon_chrom_start) + 1,
                                      as.numeric(genomic_coding_end) - as.numeric(exon_chrom_start) + 1)),
         cds_sequence_len = ifelse(genomic_coding_start == "", NA, nchar(cds_sequence)),
         gene_exon_len = nchar(gene_exon),
         genomic_coding_len = ifelse(genomic_coding_start == "", NA, as.numeric(genomic_coding_end) - as.numeric(genomic_coding_start) + 1))
f_exoncds <- f_exons_getcds %>% 
  arrange(exon_chrom_start) %>% filter(!is.na(cds_sequence)) %>%
  group_by(ensembl_gene_id, ensembl_transcript_id, strand) %>%
  summarise(exoncds = ifelse(as.numeric(strand) == 1, 
                             paste0(cds_sequence, collapse = ""),
                             as.character(complement(reverseComplement(DNAString(paste0(cds_sequence, collapse = ""))))))) %>%
  unique() %>% 
  mutate(strand = as.integer(strand))

ref_nuc %>% left_join(other_transcripts, by = c("gene" = "external_gene_name")) %>%
  # left_join(ref_seqs) %>%
  # left_join(find_cds, by = c("gene" = "external_gene_name")) %>%
  inner_join(f_exoncds) %>%
  mutate(status_nucXcoding = nuc_seq == coding,
         status_nucXexoncds = nuc_seq == exoncds,
         len_nuc = nchar(nuc_seq),
         len_coding = nchar(coding)) %>% View

#### Read through into UTR #####
# MICA ref_nuc is longer than coding sequence
# Check to see whether coding sequence is in ref_nuc
library(DECIPHER)
mica_pa <- pairwiseAlignment(pattern = filter(wrong_tx, gene == "MICA")$coding, subject = filter(wrong_tx, gene == "MICA")$nuc_seq, type = "local-global")
c(aligned(pattern(mica_pa)), aligned(subject(mica_pa))) %>% BrowseSeqs(.)
# MICA nuc seq reads into what's labeled as the 5'UTR or non-CDS part of last exon

unread_utr <- wrong_tx %>%
  filter(gene_alias %in% c("TAP2", "MICA", "DOB"))

unread_getcds <- find_exons %>% filter(ensembl_gene_id %in% unread_utr$ensembl_gene_id) %>%
  group_by(ensembl_exon_id) %>%
  mutate(gene_exon = ifelse(strand == 1, gene_exon,
                            as.character(complement(reverseComplement(DNAString(gene_exon)))))) %>% ungroup
unread_exoncds <- unread_getcds %>% 
  arrange(exon_chrom_start) %>% filter(!is.na(gene_exon)) %>%
  group_by(ensembl_gene_id, ensembl_transcript_id, strand) %>%
  summarise(all_exons = ifelse(as.numeric(strand) == 1, 
                             paste0(gene_exon, collapse = ""),
                             as.character(complement(reverseComplement(DNAString(paste0(gene_exon, collapse = ""))))))) %>%
  unique() %>% 
  mutate(strand = as.integer(strand))

all_exonseqs <- wrong_tx %>% 
  inner_join(unread_exoncds)

nuc_in_exonseq <- apply(all_exonseqs, 1, function(i){
  pa <- pairwiseAlignment(pattern = as.character(i['nuc_seq']), subject = i['all_exons'], type = "overlap")
  data.frame(pa@subject@range)
}) %>%
  enframe %>% mutate(ensembl_transcript_id = all_exonseqs$ensembl_transcript_id) %>% unnest(value) %>% mutate(name = NULL)

unread_fixgenpos <- unread_getcds %>% arrange(as.numeric(strand) * as.numeric(exon_chrom_start)) %>%
  mutate(len_gene_exon = nchar(gene_exon)) %>%
  group_by(ensembl_transcript_id) %>%
  mutate(coding_pos_end = cumsum(len_gene_exon),
         coding_pos_start = coding_pos_end - len_gene_exon + 1) %>%
  left_join(nuc_in_exonseq) %>% 
  # filter(ensembl_transcript_id != "ENST00000449934") %>%
  mutate(genomic_pos_start = ifelse(start > coding_pos_end, NA, ifelse(start >= coding_pos_start & start <= coding_pos_end,
                                                                       start, coding_pos_start)),
         genomic_pos_end = ifelse(start > coding_pos_end, NA, ifelse(end >= coding_pos_start & end <= coding_pos_end,
                                                                     end, coding_pos_end))) %>% 
  mutate(cds_chrom_start = ifelse(strand == 1,
                                    as.numeric(exon_chrom_start) + (as.numeric(genomic_pos_start) - as.numeric(coding_pos_start) + 1) - 1,
                                    as.numeric(exon_chrom_end) - (as.numeric(genomic_pos_end) - as.numeric(coding_pos_start) + 1) + 1),
         cds_chrom_end = ifelse(strand == 1,
                                  as.numeric(exon_chrom_start) + (as.numeric(genomic_pos_end) - as.numeric(coding_pos_start) + 1) - 1,
                                  as.numeric(exon_chrom_end) - (as.numeric(genomic_pos_start) - as.numeric(coding_pos_start) + 1) + 1))
unread_exons_getcds <- unread_fixgenpos %>%
  group_by(ensembl_exon_id) %>%
  mutate(exoncds_seq = ifelse(genomic_pos_start == "", NA,
                              ifelse(strand == 1,
                                     substr(gene_exon,
                                            as.numeric(genomic_pos_start) - as.numeric(coding_pos_start) + 1,
                                            as.numeric(genomic_pos_end) - as.numeric(coding_pos_start) + 1),
                                     substr(gene_exon,
                                            len_gene_exon - (as.numeric(genomic_pos_end) - as.numeric(coding_pos_start) + 1) + 1,
                                            len_gene_exon - (as.numeric(genomic_pos_start) - as.numeric(coding_pos_start) + 1) + 1)
                               )
                              )
         )
unread_exoncds <- unread_exons_getcds %>% 
  arrange(exon_chrom_start) %>% filter(!is.na(exoncds_seq)) %>%
  group_by(ensembl_gene_id, ensembl_transcript_id, strand) %>%
  summarise(exoncds = ifelse(as.numeric(strand) == 1, 
                             paste0(exoncds_seq, collapse = ""),
                             as.character(complement(reverseComplement(DNAString(paste0(exoncds_seq, collapse = ""))))))) %>%
  unique() %>% 
  mutate(strand = as.integer(strand))
unread_exoncds %>% 
  left_join(other_transcripts) %>% 
  left_join(ref_nuc, by = c("external_gene_name" = "gene")) %>%
  mutate(status_nucXcoding = as.character(nuc_seq) == as.character(coding),
         status_nucXexoncds = as.character(nuc_seq) == as.character(exoncds),
         len_nuc = nchar(nuc_seq),
         len_exoncds = nchar(exoncds)) %>% View

##### Final 5 ######
missing_codings
exons_of_missing <- find_exons %>% mutate(strand = as.integer(strand)) %>% left_join(find_transcript) %>% left_join(find_gene) %>% 
  filter(hgnc_symbol %in% missing_codings) %>% arrange(as.numeric(strand) * as.numeric(exon_chrom_start))


tm <- ref_nuc %>% filter(gene %in% missing_codings) %>% mutate(len_nuc_seq = nchar(nuc_seq))

exons_of_missing %>% group_by(ensembl_gene_id, hgnc_symbol) %>% summarise(exon_seq = paste0(gene_exon, collapse = "")) %>% 
  mutate(len_exon_seq = nchar(exon_seq)) %>% 
  left_join(tm, by = c("hgnc_symbol" = "gene")) %>%
  mutate(nuc_in_exon = grepl(nuc_seq, exon_seq),
         exon_in_nuc = grepl(exon_seq, nuc_seq))
# HLA-DPB2 nuc is exactly the exon sequences

# TODO: HLA-V, HLA-H, HLA-K, and HLA-L are currently missing from the CDS mapping
library(GenomicRanges)
remaining_seqs <- ref_seqs %>% filter(gene_alias %in% c('V','H','K','L')) %>% left_join(ref_nuc) %>%
  left_join(ref_genomicpos)

found_exons5 <- apply(remaining_seqs, 1, function(i){
  pa <- pairwiseAlignment(pattern = i['nuc_seq'], subject = i['gen_seq'], type = "global")
  full_gen <- GRanges("chr6", IRanges(c(1), c(nchar(i['gen_seq']))))
  nuc_start <- pa@subject@range@start
  nuc_end <- pa@subject@range@start + pa@subject@range@width - 1
  cds_region <- GRanges("chr6", IRanges(nuc_start, nuc_end))
  introns <- GRanges("chr6", IRanges(nuc_start + deletion(pa)@unlistData@start - 1, 
                                     nuc_start + (deletion(pa)@unlistData@start + deletion(pa)@unlistData@width - 1) - 1))
  exon_regions <- setdiff(intersect(full_gen, cds_region), introns)
  data.frame(exon_regions) %>%
    mutate(nuc_end = cumsum(width),
           nuc_start = nuc_end - width + 1,
           genomic_start = start + as.numeric(i['ref_start']) - 1,
           genomic_end = end + as.numeric(i['ref_start']), type = i['type'], strand = as.character(i['strand']),
           exon_no = paste0("exon_", 1:n()))
}) %>% enframe %>% unnest(value)
View(found_exons5)


#### Final exons ####
# 17 original
found_exons1 <- ref_nuc %>% 
  left_join(find_cds, by = c("gene" = "external_gene_name")) %>%
  left_join(exoncds) %>%
  mutate(status_nucXcoding = nuc_seq == coding,
         status_nucXexoncds = nuc_seq == exoncds) %>% 
  left_join(find_transcript) %>%
  filter(status_nucXcoding == TRUE) %>% mutate(strand = as.character(strand)) %>%
  left_join(find_exons)
# HLA-F
found_exons2 <- ref_nuc %>% 
  left_join(other_transcripts, by = c("gene" = "external_gene_name")) %>%
  inner_join(f_exoncds) %>%
  mutate(status_nucXcoding = nuc_seq == coding,
         status_nucXexoncds = nuc_seq == exoncds) %>% 
  left_join(find_transcript) %>%
  filter(status_nucXcoding == TRUE) %>% 
  left_join(f_exons)
# HLA-DOB, MICA, TAP2
found_exons3 <- ref_nuc %>%
  left_join(other_transcripts, by = c("gene" = "external_gene_name")) %>%
  inner_join(unread_exoncds) %>%
  mutate(status_nucXcoding = nuc_seq == coding,
         status_nucXexoncds = nuc_seq == exoncds) %>% 
  left_join(find_transcript) %>% mutate(strand = as.character(strand)) %>%
  left_join(unread_exons_getcds)
# HLA-DPB2
found_exons4 <- ref_nuc %>%
  left_join(other_transcripts, by = c("gene" = "external_gene_name")) %>%
  filter(gene == "HLA-DPB2") %>% mutate(strand = as.character(strand)) %>%
  inner_join(find_exons)
  
bed1 <- found_exons1 %>%
  filter(genomic_coding_start != "") %>%
  mutate(chr = "chr6", score = 0, strand = ifelse(strand == 1, "+", "-"),
         st = as.numeric(genomic_coding_start)-1, sp = as.numeric(genomic_coding_end),
         name = paste0(type, "|", ensembl_exon_id)) %>%
  select(chr, st, sp, name, score, strand)
bed2 <- found_exons2 %>% 
  filter(genomic_coding_start != "") %>%
  mutate(chr = "chr6", score = 0, strand = ifelse(strand == 1, "+", "-"),
         st = as.numeric(genomic_coding_start)-1, sp = as.numeric(genomic_coding_end),
         name = paste0(type, "|", ensembl_exon_id)) %>%
  select(chr, st, sp, name, score, strand)
bed3 <- found_exons3 %>%
  filter(!is.na(cds_chrom_start)) %>%
  mutate(chr = "chr6", score = 0, strand = ifelse(strand == 1, "+", "-"),
         st = as.numeric(cds_chrom_start)-1, sp = as.numeric(cds_chrom_end),
         name = paste0(type, "|", ensembl_exon_id, "_mod")) %>%
  select(chr, st, sp, name, score, strand)
bed4 <- found_exons4 %>% 
  mutate(chr = "chr6", score = 0, strand = ifelse(strand == 1, "+", "-"),
         st = as.numeric(exon_chrom_start)-1, sp = as.numeric(exon_chrom_end),
         name = paste0(type, "|", ensembl_exon_id)) %>%
  select(chr, st, sp, name, score, strand)
bed5 <- found_exons5 %>%
  mutate(chr = "chr6", score = 0, strand = ifelse(strand == 1, "+", "-"),
         st = as.numeric(genomic_start), sp = as.numeric(genomic_end),
         name = paste0(type, "|", exon_no)) %>%
  select(chr, st, sp, name, score, strand)

# list(bed1, bed2, bed3, bed4, bed5) %>% enframe(name = "idx") %>%
#   unnest(value) %>% mutate(idx = NULL) %>%
#   write.table(., file = "~/hlahat_final/hlahat/docker/hlahat_r/r_scripts/reftype_cds.bed", row.names = F, quote = F, sep = '\t', col.names = F)
# THIS FILE HAS BEEN MANUALLY REVIEWED. DO NOT OVERWRITE


##### END #####