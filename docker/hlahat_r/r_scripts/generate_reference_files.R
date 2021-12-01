library(dplyr)
library(tidyr)
library(Biostrings)
library(seqinr)
library(msa)
library(data.table)

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
names(nuc_msf) <- gsub(".*/([A-Z]+\\d*)_gen.msf", "\\1", list_nuc)

##### Read in haplotypes #####
summarize_hlatypes <- function(hlatypes_file, name) {
  read <- read.delim(hlatypes_file, sep = '\t', header = FALSE) %>% unlist %>% grep("ranked", ., value = T)
  summ_hlatypes <- data.frame(name = name,
                              ranks = gsub("(\\d+) ranked .+", "\\1", read) %>% as.numeric,
                              alleles = gsub(".+ ranked (.+) \\(abundance.+", "\\1", read) %>% as.character,
                              gene = gsub(".+ ranked (\\w+\\d*)\\*\\d+.+ \\(abundance.+", "\\1", read) %>% as.character,
                              perc_abundance = gsub(".+ ranked .+ \\(abundance: (\\d+\\.*\\d*)\\%\\)", "\\1", read) %>% as.numeric)
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
make_reference <- function(name){
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
                      gen_seq = gen_msf[[call[['gene']]]]$seq[which(gen_msf[[call[['gene']]]]$nam == closest_gen)],
                      nuc_seq = nuc_msf[[call[['gene']]]]$seq[which(nuc_msf[[call[['gene']]]]$nam == closest_nuc)],
                      closest_gen.nuc = nuc_msf[[call[['gene']]]]$seq[which(nuc_msf[[call[['gene']]]]$nam == closest_gen)]))
  }) %>% do.call(rbind, .)
  find_hlatypes[, alleleN := 1:.N, by = c("gene")]
  table(find_hlatypes$gen_match)
  table(find_hlatypes$nuc_match)
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
  
  # For genes with 2 alleles called, get differing position calls
  unnest_gen_seq <- find_hlatypes[, list(gen_seq = strsplit(as.character(gen_seq), "")), by = c("gene", "closest_gen", "closest_nuc", "alleleN")] %>%
    unnest(gen_seq) %>% 
    mutate(gen_seq = toupper(gen_seq)) %>%
    group_by(gene, closest_gen, closest_nuc, alleleN) %>% mutate(msf_pos = 1:n()) 
  get_ref_pos <- unnest_gen_seq %>%
    filter(gen_seq != "-") %>% mutate(ref_pos = 1:n())
  comp_gen_seq <- get_ref_pos[, c('gene','alleleN','gen_seq','msf_pos','ref_pos')] %>%
    gather(tm, value, gen_seq, ref_pos) %>%
    unite(alleleStat, alleleN, tm, sep = "-") %>%
    spread(alleleStat, value)
  # Annotate differences between alleles
  get_ref_var <- comp_gen_seq %>% filter(`1-gen_seq` != `2-gen_seq`) %>%
    group_by(gene, msf_pos) %>%
    mutate(alleleN = list(c(1, 2)), REF = list(c(`1-gen_seq`,`2-gen_seq`)), ALT = list(c(`2-gen_seq`, `1-gen_seq`)), ref_pos = list(c(`1-ref_pos`, `2-ref_pos`)), alt_pos = list(c(`2-ref_pos`,`1-ref_pos`))) %>%
    unnest(c(alleleN, REF, ALT, ref_pos, alt_pos)) %>% mutate(ref_pos = as.numeric(ref_pos), alt_pos = as.numeric(alt_pos))
  pair_alleles <- find_hlatypes %>% group_by(gene) %>% 
    summarise(closest_gen = list(c(unique(closest_gen))), closest_nuc = list(c(unique(closest_nuc)))) %>% 
    ungroup %>% mutate(alt_gen = closest_gen, alt_nuc = closest_nuc) %>% 
    unnest(c(closest_gen, closest_nuc)) %>% unnest(c(alt_gen, alt_nuc)) %>% 
    filter(closest_gen != alt_gen)
  map_allele_pos_diff <- get_ref_pos %>% 
    left_join(pair_alleles) %>%
    inner_join(get_ref_var) %>% ungroup %>%
    mutate(CHROM = closest_gen, 
           START = ref_pos-1,
           END = ref_pos,
           INFO = paste0(closest_nuc, "|", REF, ">", ALT, "|", alt_gen, "|", alt_nuc, "|", alt_pos)
    )
  write.table(map_allele_pos_diff[, c('CHROM','START','END','INFO')], 
              paste0(name, ".custom_hla.allelic_differences.bed"), 
              sep = '\t', quote = F, row.names = F, col.names = F)
}
make_reference(name)
