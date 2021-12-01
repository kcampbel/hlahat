library(dplyr)
library(tidyr)
library(Biostrings)
library(seqinr)
library(msa)
library(data.table)
library(tibble)

##### Read in arguments #####
args <- commandArgs(TRUE)
name <- args[1]
hlatypes_file <- args[2]
n_fields <- args[3] # Number of fields to resolve to
gen_msf_list <- args[4] # comma delimited list of gen_msf files
nuc_msf_list <- args[5] # comma delimited list of nuc_msf files

##### Read in reference files
list_gen <- unlist(strsplit(as.character(gen_msf_list), split = ','))
gen_msf <- lapply(list_gen, read.alignment, format = "msf")
names(gen_msf) <- gsub(".*/([A-Z]+\\d*)_gen.msf", "\\1", list_gen)
#
list_nuc <- unlist(strsplit(as.character(gen_msf_list), split = ','))
nuc_msf <- lapply(list_nuc, read.alignment, format = "msf")
names(nuc_msf) <- gsub(".*/([A-Z]+\\d*)_nuc.msf", "\\1", list_nuc)

reftype_gdna <- fread("/code/GRCh38reftype_gDNA.bed", col.names = c("chr", "st", "sp", "info", "score", "strand"))
reftype_gdna[, gene := gsub("(^\\w+\\d*)\\*\\d+:.+", "\\1", info)]
reftype_cds <- fread("/code/GRCh38reftype_cds.bed", col.names = c("chr", "st", "sp", "info", "score", "strand"))
reftype_cds[, gene := gsub("(^\\w+\\d*)\\*\\d+:.+", "\\1", info)]
#
reftype_cds_prep <- reftype_cds %>% mutate(width = sp - st) %>%
  arrange(st * ifelse(strand == "+", 1, -1)) %>%
  group_by(gene) %>% 
  mutate(nuc_pos.sp = cumsum(width), nuc_pos.st = nuc_pos.sp - width + 1)
reftype_gdnaXcds_map <- reftype_cds_prep %>% 
  left_join(reftype_gdna, by = c("chr", "gene", "score", "strand"), suffix = c(".cds", ".gdna")) %>%
  mutate(gen_pos.st = ifelse(strand == "+", (st.cds + 1) - (st.gdna + 1) + 1, sp.gdna - sp.cds + 1),
         gen_pos.sp = gen_pos.st + width - 1, gen_pos.st)
reftype_genXnuc_map <- reftype_gdnaXcds_map %>% group_by(gene, info.cds, strand, st.gdna, sp.gdna) %>%
  summarise(nuc_pos.ref = seq(nuc_pos.st, nuc_pos.sp, 1), 
            gen_pos.ref = seq(gen_pos.st, gen_pos.sp, 1)) %>%
  mutate(grch38_pos = ifelse(strand == "+", st.gdna + gen_pos.ref, sp.gdna - gen_pos.ref + 1)) %>%
  ungroup() %>% dplyr::select(gene, strand, matches("pos"))

ref_gen <- lapply(reftype_gdna$info, function(type){
  gene <- gsub("(^\\w+\\d*)\\*\\d+:.+", "\\1", type)
  alleleN = "ref"
  if( type %in% gen_msf[[gene]]$nam){
    closest_gen <- type
    gen_seq <- toupper(gen_msf[[gene]]$seq[[which(gen_msf[[gene]]$nam == closest_gen)]])
  } else if (gsub(":\\d+$", "", type) %in% gen_msf[[gene]]$nam) {
    closest_gen <- gsub(":\\d+$", "", type)
    gen_seq <- toupper(gen_msf[[gene]]$seq[[which(gen_msf[[gene]]$nam == closest_gen)]])
  } else {
    closest_gen <- NULL
    gen_seq <- NULL
  }
  return(data.frame(gene, closest_gen, gen_seq, alleleN))
}) %>% enframe %>% unnest(value) %>% 
  mutate(name = NULL) %>% data.table()

ref_nuc <- lapply(reftype_gdna$info, function(type){
  gene <- gsub("(^\\w+\\d*)\\*\\d+:.+", "\\1", type)
  alleleN = "ref"
  if( type %in% nuc_msf[[gene]]$nam){
    closest_nuc <- type
    nuc_seq <- toupper(nuc_msf[[gene]]$seq[[which(nuc_msf[[gene]]$nam == closest_nuc)]])
  } else if (gsub(":\\d+$", "", type) %in% nuc_msf[[gene]]$nam) {
    closest_nuc <- gsub(":\\d+$", "", type)
    nuc_seq <- toupper(nuc_msf[[gene]]$seq[[which(nuc_msf[[gene]]$nam == closest_nuc)]])
  } else {
    closest_nuc <- NULL
    nuc_seq <- NULL
  }
  return(data.frame(gene, closest_nuc, nuc_seq, alleleN))
}) %>% enframe %>% unnest(value) %>% 
  left_join(ref_gen) %>% mutate(name = NULL, gen_seq = NULL) %>% 
  data.table()

##### Read in haplotypes #####
summarize_hlatypes <- function(hlatypes_file, name) {
  read <- read.delim(hlatypes_file, sep = '\t', header = FALSE) %>% unlist %>% grep("ranked", ., value = T)
  if(length(read)>0){
    summ_hlatypes <- data.frame(name = name,
                                ranks = gsub("(\\d+) ranked .+", "\\1", read) %>% as.numeric,
                                alleles = gsub(".+ ranked (.+) \\(abundance.+", "\\1", read) %>% as.character,
                                gene = gsub(".+ ranked (\\w+\\d*)\\*\\d+.+ \\(abundance.+", "\\1", read) %>% as.character,
                                perc_abundance = gsub(".+ ranked .+ \\(abundance: (\\d+\\.*\\d*)\\%\\)", "\\1", read) %>% as.numeric)
  } else {
    read <- fread(hlatypes_file, header = FALSE, col.names = c("input"))
    summ_hlatypes <- read %>%
      mutate(name = name,
             alleles = gsub("HLA-", "", input),
             gene = gsub("(\\w+\\d*)\\*\\d+.+", "\\1", alleles),
             perc_abundance = 100) %>%
      group_by(gene) %>%
      mutate(ranks = 1:n()) %>% ungroup
  }
  return(summ_hlatypes)
}
read_hlatypes <- summarize_hlatypes(hlatypes_file, name)
write.table(read_hlatypes, paste0(name, ".all_hlatypes.tsv"), row.names = F, quote = F, sep = '\t')
# Use percent abundance/rank to retrieve top 2 HLA types
if(n_fields == 3) {
  top_hlatypes <- read_hlatypes %>% 
    mutate(allele = ifelse(grepl("\\w+.*\\*\\d+:\\d+:\\d+:\\d+.*", alleles), gsub("(\\w+.*\\*\\d+:\\d+:\\d+):\\d+.*", "\\1", alleles), alleles)) %>%
    group_by(gene, allele) %>% 
    filter(perc_abundance == max(perc_abundance, na.rm = T)) %>% ungroup %>% 
    filter(perc_abundance > 5) %>%
    group_by(gene) %>% mutate(new_rank = 1:n()) %>% 
    filter(new_rank %in% c(1,2))
} else {
  top_hlatypes <- read_hlatypes %>% 
    mutate(allele = ifelse(grepl("\\w+.*\\*\\d+:\\d+:\\d+:\\d+.*", alleles) , gsub("(\\w+.*\\*\\d+:\\d+):\\d+:\\d+.*", "\\1", alleles),
                           ifelse(grepl("\\w+.*\\*\\d+:\\d+:\\d+.*", alleles) , gsub("(\\w+.*\\*\\d+:\\d+.*):\\d+.*", "\\1", alleles), alleles))) %>%
    group_by(gene, allele) %>% 
    filter(perc_abundance == max(perc_abundance, na.rm = T)) %>% ungroup %>% 
    filter(perc_abundance > 5) %>%
    group_by(gene) %>% mutate(new_rank = 1:n()) %>%
    filter(new_rank %in% c(1,2))
}
write.table(top_hlatypes, paste0(name, ".top_hlatypes.tsv"), row.names = F, quote = F, sep = '\t')

# Generate reference files
# make_reference <- function(name){
# Use allele with highest level of resolution
hlatypes_list <- top_hlatypes$alleles
hlatypes_called <- data.table(hla_call = hlatypes_list)
hlatypes_called[, gene := gsub("(^.*)\\*\\d.+", "\\1", hla_call), by = 'hla_call']
hlatypes_called[, n_fields := length(unlist(strsplit(hla_call, ":"))), by = 'hla_call']
hlatypes_called <- hlatypes_called[, list(fields = strsplit(gsub("^.*\\*", "", hla_call), split = ":")), by = c('gene', 'hla_call', 'n_fields')] %>% 
  unnest(fields) %>% data.table
hlatypes_called[, field_n := 1:.N, by = c('hla_call','n_fields')]
hlatypes_called <- spread(hlatypes_called, field_n, fields)

# Find closest genomic and CDS sequences for HLA type
find_hlatypes <- apply(hlatypes_called, 1, function(call){
  hla_call <- call[['hla_call']]
  closest_gen <- NULL
  closest_nuc <- NULL
  
  # Find closest genomic DNA sequence
  if(hla_call %in% gen_msf[[call[['gene']]]]$nam) { # If genomic DNA sequence for allele is provided
    closest_gen <- hla_call
    gen_match = "EXACT"
  } else if(any(grepl(paste0(hla_call, ":"), gen_msf[[call[['gene']]]]$nam, fixed = T))) { # If field is one less another that already exists
    closest_gen <- grep(paste0(hla_call, ":"), gen_msf[[call[['gene']]]]$nam, value = T, fixed = T)[1]
    gen_match = "EXACT"
  } else {
    if(call[['n_fields']] == 1){ # If allele called only had one field
      closest_gen <- gen_msf[[call[['gene']]]]$nam[1]
      gen_match = "FIRST AVAILABLE"
    } else { # Otherwise try and find closest based upon 1 less field annotation
      less_field <- as.numeric(call[['n_fields']]) - 1
      while(is.null(closest_gen)) {
        new_call <- paste0(unlist(strsplit(hla_call, ":"))[1:less_field], collapse = ":")
        grep_new_call <- grep(new_call, gen_msf[[call[['gene']]]]$nam, value = T, fixed = T)
        if(length(grep_new_call)>0) {
          closest_gen <- grep_new_call[1]
          gen_match = "FIELD MATCH"
        } else {
          if(less_field == 1){
            closest_gen <- gen_msf[[call[['gene']]]]$nam[1]
            gen_match = "FIRST AVAILABLE"
          } else {
            less_field = less_field - 1
          }
        }
      }
    }
  }
  
  # Find closest CDS sequence
  if(hla_call %in% nuc_msf[[call[['gene']]]]$nam) { # If CDS sequence for allele is provided
    closest_nuc <- hla_call
    nuc_match = "EXACT"
  } else {
    if(call[['n_fields']] == 1){ # If allele called only had one field
      closest_nuc <- nuc_msf[[call[['gene']]]]$nam[1]
      nuc_match = "FIRST AVAILABLE"
    } else if(any(grepl(paste0(hla_call, ":"), nuc_msf[[call[['gene']]]]$nam, fixed = T))) { # If field is one less another that already exists
      closest_nuc <- grep(paste0(hla_call, ":"), nuc_msf[[call[['gene']]]]$nam, value = T, fixed = T)[1]
      nuc_match = "EXACT"
    } else { # Otherwise try and find closest based upon 1 less field annotation
      less_field <- as.numeric(call[['n_fields']]) - 1
      while(is.null(closest_nuc)) {
        new_call <- paste0(unlist(strsplit(hla_call, ":"))[1:less_field], collapse = ":")
        grep_new_call <- grep(new_call, nuc_msf[[call[['gene']]]]$nam, value = T, fixed = T)
        if(length(grep_new_call)>0) {
          closest_nuc <- grep_new_call[1]
          nuc_match = "FIELD MATCH"
        } else {
          if(less_field == 1){
            closest_nuc <- nuc_msf[[call[['gene']]]]$nam[1]
            nuc_match = "FIRST AVAILABLE"
          } else {
            less_field = less_field - 1
          }
        }
      }
    }
  }
  return(data.table(gene = call[['gene']],
                    hla_call = hla_call,
                    closest_gen = closest_gen,
                    closest_nuc = closest_nuc,
                    gen_match = gen_match,
                    nuc_match = nuc_match,
                    gen_seq = toupper(gen_msf[[call[['gene']]]]$seq[which(gen_msf[[call[['gene']]]]$nam == closest_gen)]),
                    nuc_seq = toupper(nuc_msf[[call[['gene']]]]$seq[which(nuc_msf[[call[['gene']]]]$nam == closest_nuc)]),
                    closest_gen.nuc = toupper(nuc_msf[[call[['gene']]]]$seq[which(nuc_msf[[call[['gene']]]]$nam == closest_gen)])))
}) %>% do.call(rbind, .)
find_hlatypes[, alleleN := 1:.N, by = c("gene")]
write.table(find_hlatypes[, c('gene','hla_call','alleleN','closest_gen','closest_nuc','gen_match','nuc_match')],
file = paste0(name, ".find_hlatypes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

##### STEP 1. Create Custom HLA Reference #####
# Report which HLA types are duplicates
hla_gen_fasta <- find_hlatypes %>% group_by(closest_gen) %>%
  summarise(hla_calls = list(hla_call),
            gen_seq = toupper(gsub("-", "", unique(gen_seq))))
write.fasta(sequences = as.list(hla_gen_fasta$gen_seq), # sequences
            names = hla_gen_fasta$closest_gen, # names
            paste0(name, ".custom_hla.fasta"))

##### STEP 2. Create Exon BED file for custom reference #####
allgen <- find_hlatypes %>% dplyr::select(gene, closest_gen, gen_seq) %>% 
  left_join(ref_gen, by = c("gene"), suffix = c(".custom", ".ref")) %>% mutate(alleleN = NULL) %>%
  gather(key, gen_seq, matches("gen_seq"))
allgen_unnest <- allgen %>%
  group_by(gene, closest_gen.custom, closest_gen.ref, key) %>%
  mutate(gen_seq = strsplit(as.character(gen_seq), split = "")) %>% unnest(gen_seq) %>%
  mutate(msf_pos = 1:n()) %>% 
  filter(gen_seq != "-") %>% mutate(gen_pos = 1:n())
allgen_compseq <- allgen_unnest %>% 
  gather(key2, value, gen_seq, gen_pos) %>%
  unite(key, key, key2, sep = ".") %>% spread(key, value)
# Map exon locations to custom gDNA for custom exon bed
exonbed_prep <- reftype_gdnaXcds_map %>% 
  mutate(gen_pos.st = as.character(gen_pos.st), gen_pos.sp = as.character(gen_pos.sp)) %>%
  left_join(allgen_compseq, by = c("gene", "gen_pos.st" = "gen_seq.ref.gen_pos")) %>%
  left_join(allgen_compseq, by = c("gene", "gen_pos.sp" = "gen_seq.ref.gen_pos", "closest_gen.custom", "closest_gen.ref"), suffix = c(".st", ".sp")) %>%
  filter(!is.na(gen_seq.custom.gen_pos.st) & !is.na(gen_seq.custom.gen_pos.sp))
exonbed.custom <- data.frame(chrom = exonbed_prep$closest_gen.custom,
                                  st = as.numeric(exonbed_prep$gen_seq.custom.gen_pos.st)-1,
                                  sp = as.numeric(exonbed_prep$gen_seq.custom.gen_pos.sp),
                                  info = paste(exonbed_prep$closest_gen.custom, exonbed_prep$info.cds, sep = "_"),
                                  score = 0,
                                  strand = "+")
# Create custom exon bed with GRCh38 coordinates
exonbed.ref <- data.frame(chrom = "chr6",
                                  st = as.numeric(exonbed_prep$st.cds),
                                  sp = as.numeric(exonbed_prep$sp.cds),
                                  info = paste(exonbed_prep$closest_gen.custom, exonbed_prep$info.cds, sep = "_"),
                                  score = 0,
                                  strand = exonbed_prep$strand)

##### STEP 3. Create SNP bed files #####
## Scenario 1: All alleles have gDNA available
scen1 <- find_hlatypes %>% group_by(gene) %>% 
  filter(all(gen_match == "EXACT")) %>% mutate(alleleN = as.character(alleleN))
if(nrow(scen1)>0){
  scen1_prep <- scen1 %>% dplyr::select(gene, closest_gen, gen_seq, alleleN) %>% 
    full_join(ref_gen) %>% filter(gene %in% scen1$gene)
  scen1_unnest <- scen1_prep %>%
    group_by(gene, alleleN) %>% mutate(gen_seq = strsplit(as.character(gen_seq), split = "")) %>% unnest(gen_seq) %>%
    mutate(msf_pos = 1:n()) %>% 
    filter(gen_seq != "-") %>% mutate(gen_pos = 1:n())
  scen1_compseq <- scen1_unnest %>% 
    gather(key, value, closest_gen, gen_seq, gen_pos) %>%
    unite(key, key, alleleN, sep = ".") %>% spread(key, value) %>% 
    inner_join(reftype_gdna) %>%
    mutate(ref_pos.ref = ifelse(strand == "+", st + as.numeric(gen_pos.ref), sp - as.numeric(gen_pos.ref) + 1))
  ## Obtain variants 1 vs. ref
  scen1_1ref <- scen1_compseq %>% filter(gen_seq.1 != gen_seq.ref)
  if(nrow(scen1_1ref)>0){
    scen1_1ref.custom <- data.frame(chrom = scen1_1ref$closest_gen.1, st = as.numeric(scen1_1ref$gen_pos.1) - 1, sp = as.numeric(scen1_1ref$gen_pos.1),
                                    info = paste("1vref", paste0(scen1_1ref$closest_gen.1, ">", scen1_1ref$closest_gen.ref), paste0(scen1_1ref$gen_seq.1, ">", scen1_1ref$gen_seq.ref), sep = "|"),
                                    score = 0, strand = "+")
    scen1_1ref.grch38 <- data.frame(chrom = scen1_1ref$chr, st = as.numeric(scen1_1ref$ref_pos.ref) - 1, sp = as.numeric(scen1_1ref$ref_pos.ref),
                                    info = paste("refv1", paste0(scen1_1ref$closest_gen.ref, ">", scen1_1ref$closest_gen.1), paste0(scen1_1ref$gen_seq.ref, ">", scen1_1ref$gen_seq.1), sep = "|"),
                                    score = 0, strand = scen1_1ref$strand)
  }
  #
  if("gen_seq.2" %in% colnames(scen1_compseq)) {
    ## Obtain variants 2 vs. ref
    scen1_2ref <- scen1_compseq %>% filter(gen_seq.2 != gen_seq.ref)
    if(nrow(scen1_2ref)>0){
      scen1_2ref.custom <- data.frame(chrom = scen1_2ref$closest_gen.2, st = as.numeric(scen1_2ref$gen_pos.2) - 1, sp = as.numeric(scen1_2ref$gen_pos.2),
                                      info = paste("2vref", paste0(scen1_2ref$closest_gen.2, ">", scen1_2ref$closest_gen.ref), paste0(scen1_2ref$gen_seq.2, ">", scen1_2ref$gen_seq.ref), sep = "|"),
                                      score = 0, strand = "+")
      scen1_2ref.grch38 <- data.frame(chrom = scen1_2ref$chr, st = as.numeric(scen1_2ref$ref_pos.ref) - 1, sp = as.numeric(scen1_2ref$ref_pos.ref),
                                      info = paste("refv2", paste0(scen1_2ref$closest_gen.ref, ">", scen1_2ref$closest_gen.2), paste0(scen1_2ref$gen_seq.ref, ">", scen1_2ref$gen_seq.2), sep = "|"),
                                      score = 0, strand = scen1_2ref$strand)
    }
    ## Obtain variants 1 vs. 2
    scen1_12 <- scen1_compseq %>% filter(gen_seq.1 != gen_seq.2)
    if(nrow(scen1_12)>0){
      scen1_12.custom <- data.frame(chrom = c(scen1_12$closest_gen.1, scen1_12$closest_gen.2), st = c(as.numeric(scen1_12$gen_pos.1) - 1, as.numeric(scen1_12$gen_pos.2) - 1), sp = c(as.numeric(scen1_12$gen_pos.1), as.numeric(scen1_12$gen_pos.2)),
                                    info = c(paste("1v2", paste0(scen1_12$closest_gen.1, ">", scen1_12$closest_gen.2), paste0(scen1_12$gen_seq.1, ">", scen1_12$gen_seq.2), sep = "|"),
                                             paste("2v1", paste0(scen1_12$closest_gen.2, ">", scen1_12$closest_gen.1), paste0(scen1_12$gen_seq.2, ">", scen1_12$gen_seq.1), sep = "|")),
                                    score = 0, strand = "+")
      scen1_12.grch38 <- data.frame(chrom = scen1_12$chr, st = as.numeric(scen1_12$ref_pos.ref) - 1, sp = as.numeric(scen1_12$ref_pos.ref),
                                    info = paste("1v2", paste0(scen1_12$closest_gen.1, ">", scen1_12$closest_gen.2), paste0(scen1_12$gen_seq.1, ">", scen1_12$gen_seq.2), sep = "|"),
                                    score = 0, strand = scen1_12$strand)
    }
  }
}

## Scenario 2: One or both alleles don't have gDNA available
scen2 <- find_hlatypes %>% group_by(gene) %>% 
  filter(any(gen_match != "EXACT")) %>% mutate(alleleN = as.character(alleleN))
if(nrow(scen2) > 0) {
  scen2_prep.gen <- scen2 %>% dplyr::select(gene, closest_gen, gen_seq, alleleN) %>% 
    full_join(ref_gen) %>% filter(gene %in% scen2$gene)
  scen2_unnest.gen <- scen2_prep.gen %>%
    group_by(gene, alleleN) %>% mutate(gen_seq = strsplit(as.character(gen_seq), split = "")) %>% unnest(gen_seq) %>%
    mutate(gen.msf_pos = 1:n()) %>% 
    filter(gen_seq != "-") %>% mutate(gen_pos = 1:n())
  scen2_compseq.gen <- scen2_unnest.gen %>% 
    gather(key, value, closest_gen, gen_seq, gen_pos) %>%
    unite(key, key, alleleN, sep = ".") %>% spread(key, value) %>% mutate(gen_pos.ref = as.numeric(gen_pos.ref))
  #
  scen2_prep <- scen2 %>% dplyr::select(gene, closest_gen, closest_nuc, nuc_seq, closest_gen.nuc, alleleN) %>%
    gather(key, nuc_seq, nuc_seq, closest_gen.nuc) %>%
    unite(alleleN, alleleN, key, sep = ".") %>%
    full_join(ref_nuc) %>% filter(gene %in% scen2$gene)
  scen2_unnest <- scen2_prep %>%
    group_by(gene, alleleN) %>% mutate(nuc_seq = strsplit(as.character(nuc_seq), split = "")) %>% unnest(nuc_seq) %>%
    mutate(msf_pos = 1:n()) %>% 
    filter(nuc_seq != "-") %>% mutate(nuc_pos = 1:n())
  scen2_compseq <- scen2_unnest %>% 
    gather(key, value, closest_gen, closest_nuc, nuc_seq, nuc_pos) %>%
    unite(key, key, alleleN, sep = ".") %>% spread(key, value) %>% mutate(nuc_pos.ref = as.numeric(nuc_pos.ref)) %>%
    left_join(reftype_genXnuc_map) %>%
    left_join(scen2_compseq.gen, by = c("gene", "closest_gen.1.closest_gen.nuc" = "closest_gen.1", "closest_gen.2.closest_gen.nuc" = "closest_gen.2", "closest_gen.ref",
                                        "gen_pos.ref"))
  ## Obtain variants 1 vs. ref
  scen2_1ref <- scen2_compseq %>% filter(nuc_seq.1.nuc_seq != nuc_seq.ref)
  if(nrow(scen2_1ref)>0){
    scen2_1ref.custom <- data.frame(chrom = scen2_1ref$closest_gen.1.closest_gen.nuc, st = as.numeric(scen2_1ref$gen_pos.1) - 1, sp = as.numeric(scen2_1ref$gen_pos.1),
                                    info = paste("1vref", paste0(scen2_1ref$closest_nuc.1.nuc_seq, ">", scen2_1ref$closest_nuc.ref), paste0(scen2_1ref$nuc_seq.1.nuc_seq, ">", scen2_1ref$nuc_seq.ref), sep = "|"),
                                    score = 0, strand = "+")
    scen2_1ref.grch38 <- data.frame(chrom = "chr6", st = as.numeric(scen2_1ref$grch38_pos) - 1, sp = as.numeric(scen2_1ref$grch38_pos),
                                    info = paste("refv1", paste0(scen2_1ref$closest_nuc.ref, ">", scen2_1ref$closest_nuc.1.nuc_seq), paste0(scen2_1ref$nuc_seq.ref, ">", scen2_1ref$nuc_seq.1.nuc_seq), sep = "|"),
                                    score = 0, strand = scen2_1ref$strand)
  }
  # If there are differences between closest_nuc 1 and closest_gen 1
  scen2_1closestgen <- scen2_compseq %>% filter(nuc_seq.1.nuc_seq != nuc_seq.1.closest_gen.nuc)
  if(nrow(scen2_1closestgen)>0){
    scen2_1closestgen.custom <- data.frame(chrom = scen2_1closestgen$closest_gen.1.closest_gen.nuc, st = as.numeric(scen2_1closestgen$gen_pos.1) - 1, sp = as.numeric(scen2_1closestgen$gen_pos.1),
                                    info = paste("1", paste0(scen2_1closestgen$closest_gen.1.closest_gen.nuc, ">", scen2_1closestgen$closest_nuc.1.nuc_seq), paste0(scen2_1closestgen$nuc_seq.1.closest_gen.nuc, ">", scen2_1closestgen$nuc_seq.1.nuc_seq), sep = "|"),
                                    score = 0, strand = "+")
  }
  if("gen_seq.2" %in% colnames(scen2_compseq)) {
    ## Obtain variants 2 vs. ref
    scen2_2ref <- scen2_compseq %>% filter(nuc_seq.2.nuc_seq != nuc_seq.ref)
    if(nrow(scen2_2ref)>0){
      scen2_2ref.custom <- data.frame(chrom = scen2_2ref$closest_gen.2.closest_gen.nuc, st = as.numeric(scen2_2ref$gen_pos.2) - 1, sp = as.numeric(scen2_2ref$gen_pos.2),
                                      info = paste("2vref", paste0(scen2_2ref$closest_nuc.2.nuc_seq, ">", scen2_2ref$closest_nuc.ref), paste0(scen2_2ref$nuc_seq.2.nuc_seq, ">", scen2_2ref$nuc_seq.ref), sep = "|"),
                                      score = 0, strand = "+")
      scen2_2ref.grch38 <- data.frame(chrom = "chr6", st = as.numeric(scen2_2ref$grch38_pos) - 1, sp = as.numeric(scen2_2ref$grch38_pos),
                                      info = paste("refv2", paste0(scen2_2ref$closest_nuc.ref, ">", scen2_2ref$closest_nuc.2.nuc_seq), paste0(scen2_2ref$nuc_seq.ref, ">", scen2_2ref$nuc_seq.2.nuc_seq), sep = "|"),
                                      score = 0, strand = scen2_2ref$strand)
    }
    # If there are differences between closest_nuc 2 and closest_gen 2
    scen2_2closestgen <- scen2_compseq %>% filter(nuc_seq.2.nuc_seq != nuc_seq.2.closest_gen.nuc)
    if(nrow(scen2_2closestgen)>0){
      scen2_2closestgen.custom <- data.frame(chrom = scen2_2closestgen$closest_gen.2.closest_gen.nuc, st = as.numeric(scen2_2closestgen$gen_pos.2) - 1, sp = as.numeric(scen2_2closestgen$gen_pos.2),
                                             info = paste("2", paste0(scen2_2closestgen$closest_gen.2.closest_gen.nuc, ">", scen2_2closestgen$closest_nuc.2.nuc_seq), paste0(scen2_2closestgen$nuc_seq.2.closest_gen.nuc, ">", scen2_2closestgen$nuc_seq.2.nuc_seq), sep = "|"),
                                             score = 0, strand = "+")
    }
    # Obtain variants 1 v 2
    scen2_12 <- scen2_compseq %>% filter(nuc_seq.1.nuc_seq != nuc_seq.2.nuc_seq)
    if(nrow(scen2_12)>0){
      scen2_12.custom <- data.frame(chrom = c(scen2_12$closest_gen.1.closest_gen.nuc, scen2_12$closest_gen.2.closest_gen.nuc), st = c(as.numeric(scen2_12$gen_pos.1) - 1, as.numeric(scen2_12$gen_pos.2) - 1), sp = c(as.numeric(scen2_12$gen_pos.1), as.numeric(scen2_12$gen_pos.2)),
                                    info = c(paste("1v2", paste0(scen2_12$closest_nuc.1.nuc_seq, ">", scen2_12$closest_nuc.2.nuc_seq), paste0(scen2_12$nuc_seq.1.nuc_seq, ">", scen2_12$nuc_seq.2.nuc_seq), sep = "|"),
                                             paste("2v1", paste0(scen2_12$closest_nuc.2.nuc_seq, ">", scen2_12$closest_nuc.1.nuc_seq), paste0(scen2_12$nuc_seq.2.nuc_seq, ">", scen2_12$nuc_seq.1.nuc_seq), sep = "|")),
                                    score = 0, strand = "+")
      scen2_12.grch38 <- data.frame(chrom = "chr6", st = as.numeric(scen2_12$grch38_pos) - 1, sp = as.numeric(scen2_12$grch38_pos),
                                    info = paste("1v2", paste0(scen2_12$closest_nuc.1.nuc_seq, ">", scen2_12$closest_nuc.2.nuc_seq), paste0(scen2_12$nuc_seq.1.nuc_seq, ">", scen2_12$nuc_seq.2.nuc_seq), sep = "|"),
                                    score = 0, strand = scen2_12$strand)
    }
  }
}

##### WRITE OUTPUTS #####
# Exon bed file: Custom genome
exonbed.custom %>% arrange(chrom, st, sp) %>%
  write.table(., file = paste0(name, ".custom_hla.cds_exons.bed"), sep = '\t', row.names = F, quote = F, col.names = F)
# Exon bed file: GRCh38
exonbed.ref %>% arrange(chrom, st, sp) %>%
  write.table(., file = paste0(name, ".GRCh38.cds_exons.bed"), sep = '\t', row.names = F, quote = F, col.names = F)

# SNPs: Custom genome
ls(pattern = "\\.custom\\b") %>%
  grep("bed", ., invert = T, value = T) %>%
  lapply(., get) %>% enframe %>% unnest(value) %>% mutate(name = NULL) %>% 
  arrange(chrom, st, sp) %>% 
  write.table(., file = paste0(name, ".custom_hla.SNPs.bed"), sep = '\t', row.names = F, quote = F, col.names = F)
# SNPs: GRCh38
ls(pattern = "\\.grch38\\b") %>%
  grep("bed", ., invert = T, value = T) %>%
  lapply(., get) %>% enframe %>% unnest(value) %>% mutate(name = NULL) %>% 
  arrange(chrom, st, sp) %>% filter(!is.na(st)) %>%
  write.table(., file = paste0(name, ".GRCh38.SNPs.bed"), sep = '\t', row.names = F, quote = F, col.names = F)

# SNPs: Custom genome, exon only
custom_gen_pos_inexons <- allgen_compseq %>% 
  mutate(gen_seq.ref.gen_pos = as.numeric(gen_seq.ref.gen_pos), gen_seq.custom.gen_pos = as.numeric(gen_seq.custom.gen_pos)) %>% 
  inner_join(reftype_genXnuc_map, by = c("gene", "gen_seq.ref.gen_pos" = "gen_pos.ref")) %>% ungroup %>%
  dplyr::select(closest_gen.custom, gen_seq.custom.gen_pos)
ls(pattern = "\\.custom\\b") %>%
  grep("bed", ., invert = T, value = T) %>%
  lapply(., get) %>% enframe %>% unnest(value) %>% mutate(name = NULL) %>%
  arrange(chrom, st, sp) %>% 
  inner_join(custom_gen_pos_inexons, by = c("chrom" = "closest_gen.custom", "sp" = "gen_seq.custom.gen_pos")) %>%
  write.table(., file = paste0(name, ".custom_hla.SNPs_inexons.bed"), sep = '\t', row.names = F, quote = F, col.names = F)

# SNPs: GRCh38, exon only
ls(pattern = "\\.grch38\\b") %>%
  grep("bed", ., invert = T, value = T) %>%
  lapply(., get) %>% enframe %>% unnest(value) %>% mutate(name = NULL) %>%
  arrange(chrom, st, sp) %>% filter(!is.na(st)) %>%
  filter(sp %in% reftype_genXnuc_map$grch38_pos) %>%
  write.table(., file = paste0(name, ".GRCh38.SNPs_inexons.bed"), sep = '\t', row.names = F, quote = F, col.names = F)

##### END #####

