library(seqinr)
library(Biostrings)
library(biomaRt)
library(tidyverse)
library(msa)

#### USER SPECIFIC: SET WD OF MSF FILES
msf_dir <- "~/hisat_test/tm/hisat2-hisat2_v2.2.0_beta/hisatgenotype_db/HLA/msf/"

#### Hard coded  #####
hlas <- c('HFE','F','V','G','H',
          'K','A','L','E','C',
          'B','MICA','MICB','DRA','DRB1',
          'DQA1','DQB1','DOB','TAP2','TAP1',
          'DMB','DMA','DOA','DPA1','DPB1',
          'DPB2')
hla_genes <- gsub("(\\bD\\w.+\\d*)", "HLA-\\1", gsub("(\\b\\w\\b)", "HLA-\\1", hlas))
ref_hla_types <- c('HFE*001:01:01','F*01:03:01:01','V*01:01:01:01','G*01:01:01:05','H*02:04:01',
                 'K*01:01:01:01','A*03:01:01:01','L*01:01:01:03','E*01:03:02:01','C*07:02:01:03',
                 'B*07:02:01:01','MICA*008:04:01','MICB*004:01:01','DRA*01:02:03','DRB1*15:01:01:01',
                 'DQA1*01:02:01:01','DQB1*06:02:01:01','DOB*01:01:01:01','TAP2*01:01:03:01','TAP1*01:01:01:05',
                 'DMB*01:03:01:03','DMA*01:01:01:03','DOA*01:01:02:02','DPA1*01:03:01:02','DPB1*04:01:01:01',
                 'DPB2*03:01:01:01')
ref_hla <- data.frame(id = hlas, gene = hla_genes, type = ref_hla_types)

###### Genomic reference sequences #####
# Ensembl mart
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl", version = "Ensembl Genes 102", host = "https://nov2020.archive.ensembl.org")
# HLA genes
find_gene <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
                            filters = c("hgnc_symbol", "chromosome_name"),
                            values = list(hgnc_symbol = ref_hla$gene, chromosome_name = "6"),
                            mart = ensembl)
# HLA transcripts
find_transcript <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_start", "transcript_end", "strand"),
                                  filters = c("ensembl_gene_id", "chromosome_name"),
                                  values = list(ensembl_gene_id = find_gene$ensembl_gene_id, chromosome_name = "6"),
                                  mart = ensembl)

# Transcript flank starts at beginning of transcript (including UTR)
upstream_transcript_flank <- biomaRt::getSequence(id = find_transcript$ensembl_transcript_id,
                                                  type = "ensembl_transcript_id",
                                                  seqType = "transcript_flank",
                                                  upstream = 1000,
                                                  mart = ensembl)
colnames(upstream_transcript_flank)[which(colnames(upstream_transcript_flank) == "transcript_flank")] <- "upstream_transcript_flank"
downstream_transcript_flank <- biomaRt::getSequence(id = find_transcript$ensembl_transcript_id,
                                                    type = "ensembl_transcript_id",
                                                    seqType = "transcript_flank",
                                                    downstream = 1000,
                                                    mart = ensembl)
colnames(downstream_transcript_flank)[which(colnames(downstream_transcript_flank) == "transcript_flank")] <- "downstream_transcript_flank"
flank_info <- find_gene %>% full_join(find_transcript) %>% 
  full_join(upstream_transcript_flank) %>% 
  full_join(downstream_transcript_flank)

# Transcript sequence
transcript_seq <- biomaRt::getSequence(id = find_transcript$ensembl_transcript_id,
                                       type = "ensembl_transcript_id",
                                       seqType = "transcript_exon_intron",
                                       mart = ensembl)
transcript_seq <- find_gene %>% full_join(find_transcript) %>% 
  full_join(transcript_seq) %>% mutate(length = nchar(transcript_exon_intron))
# Get reference genome HLA types using pairwiseAlignment of reference_fastas to 
transcripts <- left_join(transcript_seq, flank_info) %>%
  mutate(context_seq = paste0(upstream_transcript_flank, transcript_exon_intron, downstream_transcript_flank))

##### Map Type to Reference Genome to get IMGT coordinates ####
gens <- paste0(msf_dir, ref_hla$id, "_gen.msf")
gen_msf <- lapply(gens, function(msf) { read.alignment(msf, format = "msf") })
names(gen_msf) <- ref_hla$id

apply(ref_hla, 1, function(hla){
  # Find MSF in genomic DNA
  find_msf <- grep(hla['type'], gen_msf[[hla['id']]]$nam, fixed = TRUE)
  # If not found, take away one field to find closest genomic DNA sequence
  if(length(find_msf)==0) {
    take_away_field <- gsub("(.+):\\d+\\b", "\\1", hla['type'])
    find_msf <- grep(take_away_field, gen_msf[[hla['id']]]$nam, fixed = TRUE)
  }
})





##### OLD #####


#####Get reference HLA type sequences #####
# Check if gen sequence provided
gens <- paste0("~/immutype/ref/hla/msf/", sortedHlaGenes, "_gen.msf")
gen_msf <- lapply(gens, function(msf) { read.alignment(msf, format = "msf") })
names(gen_msf) <- sortedHlaGenes
nucs <- paste0("~/immutype/ref/hla/msf/", sortedHlaGenes, "_nuc.msf")
nuc_msf <- lapply(nucs, function(msf) { read.alignment(msf, format = "msf") })
names(nuc_msf) <- sortedHlaGenes

find_type <- function(gene, allele){
  ref <- reference_hlatypes %>% filter(hgnc_symbol == gsub("(\\bD\\w.+\\d*)", "HLA-\\1", gsub("(\\b\\w\\b)", "HLA-\\1", gene)))
  gen <- gen_msf[[gene]]
  nuc <- nuc_msf[[gene]]
  # 1. If allele matches reference genotype, no changes need to be made to the reference genome
  if(allele %in% ref$allele) {
    # warning(paste(allele, "matches the reference genome. No changes will be made."))
    # return(paste0(allele, "REF"))
    return("REF")
    
  # 2. Check if gDNA is available for allele
  } else if(allele %in% gen$nam) {
    # return(paste0(allele, "gDNA"))
    return("gDNA")
  # 3. Check if allele has a field matched gDNA available
  } else if( any( grepl(allele, gen$nam, fixed = TRUE) ) ) {
    gDNA_status <- "FIELD"
    gDNA_hlatype <- grep(paste0(allele, ":"), gen$nam, fixed = TRUE, value = T)[1]
    # return(paste0(allele, "Field-matched gDNA"))
    return(grep(paste0(allele, ":"), gen$nam, fixed = TRUE, value = T)[1])
  
  # 4. Check if cDNA for allele is available
  } else if(allele %in% nuc$nam) {
    # return(paste0(allele, "cDNA"))
    return("cDNA")
    
  # 5. Lose fields to try and find closest match. Hopefully not necessary
  } else {
    gDNA_hlatype <- "tmp"
    lose_field <- "tmp"
    while( gDNA_hlatype == "tmp" ) {
      previous_field <- lose_field
      lose_field <- gsub("(\\w+\\d*\\*\\d.+):\\d+.*\\b", "\\1", allele)
      find_field <- grep(paste0(lose_field, ":"), union(gen$nam, nuc$nam), fixed = TRUE, value = T)
      if (length(find_field)>0) {
        gDNA_status <- "FIELD MATCH"
        gDNA_hlatype <- find_field[1]
      } else if (previous_field == lose_field) {
        gDNA_status <- "FIRST gDNA AVAILABLE"
        gDNA_hlatype <- gen$nam[1]
      } else {
        warning("Neither condition")
      }
    }
    warning(paste(allele, "was not found.", gDNA_hlatype, "will be used for gDNA."))
  }
}
chosen$status <- mapply(find_type, gene = chosen$gene, allele = chosen$haplotype, SIMPLIFY = FALSE) %>%
  as.character()

##### Prepare reform inputs #####
tmp <- chosen %>% mutate(row = 1:n()) %>% filter(status == "cDNA")
lapply(tmp$row, function(ix) {
gene = chosen[ix,'gene']
gene
allele = chosen[ix, 'haplotype']
allele
ref <- reference_hlatypes %>% filter(hgnc_symbol == gsub("(\\bD\\w.+\\d*)", "HLA-\\1", gsub("(\\b\\w\\b)", "HLA-\\1", gene)))
strand <- unique(filter(transcript_seq, hgnc_symbol == gsub("(\\bD\\w.+\\d*)", "HLA-\\1", gsub("(\\b\\w\\b)", "HLA-\\1", gene)))$strand)
strand
gen <- gen_msf[[gene]]
nuc <- nuc_msf[[gene]]
this_nuc <- nuc$seq[which(nuc$nam == allele)] %>% toupper
ref_nuc <- nuc$seq[which(nuc$nam == ref$allele)] %>% toupper
ref_gen <- gen$seq[which(gen$nam == ref$allele)] %>% toupper
pa_nuc <- pairwiseAlignment(pattern = gsub("-", "", this_nuc), subject = gsub("-", "", ref_gen), type = "local-global")
### Get location of CDS in gDNA
range_loc <- pa_nuc@subject@range
forGtf_CDS <- data.frame(seqname = "6 ",
                         source = gsub("\\*", "_", gsub(":", "_", allele)),
                         feature = "CDS",
                         start = ifelse(strand == 1, range_loc@start, nchar(this_gen) - (range_loc@start + range_loc@width - 1) + 1),
                         end = ifelse(strand == 1, range_loc@start + range_loc@width - 1, nchar(this_gen) - range_loc@start +1),
                         score = ".",
                         strand = ifelse(strand == 1, "+", "-"),
                         frame = "0",
                         attribute = paste0('gene_id "', gsub("\\*", "_", gsub(":", "_", allele)),
                                            '"; transcript_id "', gsub("\\*", "_", gsub(":", "_", allele)),
                                            '";'))
### Get exon locations
cds_in_gDNA <- IRanges(start = range_loc@start, 
                       end = range_loc@start + range_loc@width - 1)
# Intron positions
introns <- as.data.frame(deletion(pa_nuc)) %>% 
  mutate(cumwidth = cumsum(width),
         newstart = start + lag(cumwidth, default = 0),
         newend = newstart + width - 1)
introns_in_gDNA <- IRanges(start = introns$newstart + range_loc@start,
                           end = introns$newend + range_loc@start)
# Exon locations
exon_locs <- setdiff.Vector(cds_in_gDNA, introns_in_gDNA)
if(strand == 1) {
  forGtf_exon <- data.frame(seqname = "6",
                            source = gsub("\\*", "_", gsub(":", "_", allele)),
                            feature = "exon",
                            start = exon_locs@start,
                            end = exon_locs@start + exon_locs@width - 1,
                            score = ".",
                            strand = ifelse(strand == 1, "+", "-"),
                            frame = "0",
                            attribute = paste0('gene_id "', gsub("\\*", "_", gsub(":", "_", allele)),
                                               '"; transcript_id "', gsub("\\*", "_", gsub(":", "_", allele)),
                                               '"; exon_id "', gsub("\\*", "_", gsub(":", "_", allele)), ".", 1:length(exon_locs), '";'))
} else {
  forGtf_exon <- data.frame(seqname = "6",
                            source = gsub("\\*", "_", gsub(":", "_", allele)),
                            feature = "exon",
                            start = nchar(this_gen) - (exon_locs@start + exon_locs@width - 1) + 1,
                            end = nchar(this_gen) - exon_locs@start + 1,
                            score = ".",
                            strand = "-",
                            frame = "0",
                            attribute = paste0('gene_id "', gsub("\\*", "_", gsub(":", "_", allele)),
                                               '"; transcript_id "', gsub("\\*", "_", gsub(":", "_", allele)),
                                               '"; exon_id "', gsub("\\*", "_", gsub(":", "_", allele)), ".", 1:length(exon_locs), '";'))
  
}
### Get Up/Downstream sequences for allele
trans <- ref %>% inner_join(summ_IMGT_to_ref) %>% inner_join(transcripts)
trans_pa <- apply(trans, 1, function(tx){
  if(strand == 1){
    gen_context_pa <- pairwiseAlignment(pattern = gsub("-", "", ref_gen), # Need to use location of ref allele gDNA
                                        subject = as.character(tx['context_seq']),
                                        type = "global-local")
  } else {
    gen_context_pa <- pairwiseAlignment(pattern = as.character(reverseComplement(DNAString(gsub("-", "", ref_gen)))), # Need to use location of ref allele gDNA
                                        subject = as.character(reverseComplement(DNAString(tx['context_seq']))),
                                        type = "global-local")
  }
  return(gen_context_pa)
})

trans$ref_context_start <- sapply(trans_pa, function(i){ i@subject@range@start })
trans$ref_context_end <- sapply(trans_pa, function(i){ i@subject@range@start }) + sapply(trans_pa, function(i){ i@subject@range@width })
flanks <- apply(trans, 1, function(tx){
  start_in_context <- as.numeric(tx['ref_context_start'])
  end_in_context <- as.numeric(tx['ref_context_end'])
  if(strand == 1) {
    upstream_100 <- substr(tx['context_seq'], start_in_context-100, start_in_context-1)
    downstream_100 <- substr(tx['context_seq'], end_in_context, end_in_context+100-1)
  } else {
    upstream_100 <- substr(as.character(reverseComplement(DNAString(tx['context_seq']))), start_in_context-100, start_in_context-1)
    downstream_100 <- substr(as.character(reverseComplement(DNAString(tx['context_seq']))), end_in_context, end_in_context+100-1)
  }
  data.frame(up_flank = upstream_100, down_flank = downstream_100)
}) %>% do.call(rbind, .) %>% filter(nchar(up_flank) == 100 & nchar(down_flank)==100) %>% unique

# Replace CDS differences using the mapping between pa_nuc and 
this_msa <- msa(c(as.character(pa_nuc@subject), as.character(pa_nuc@pattern)), type = "dna")
# print(this_msa, show = "complete")
correct_sequence <- data.frame(gen = as.character(pa_nuc@subject),
           nuc = as.character(pa_nuc@pattern),
           cons = msaConsensusSequence(this_msa)) %>%
  gather(key, value) %>% group_by(key) %>% mutate(value = strsplit(value, split = "")) %>%
  unnest(value) %>% mutate(i = 1:n()) %>% spread(key, value) %>%
  mutate(pick_value = ifelse(cons != "?", cons, ifelse(nuc == "-", gen, nuc)))

beg <- substr(as.character(pa_nuc@subject@unaligned), 1, cds_in_gDNA@start-1)
mid <- paste0(correct_sequence$pick_value, collapse = "")
end <- substr(as.character(pa_nuc@subject@unaligned), cds_in_gDNA@start + cds_in_gDNA@width, nchar(as.character(pa_nuc@subject@unaligned)))
# Write in FASTA
write.fasta(ifelse(strand == 1,
                   list(paste0(beg, mid, end)),
                   list(as.character(reverseComplement(DNAString(paste0(beg, mid, end)))))),
            gsub("\\*", "_", gsub(":", "_", allele)),
            paste0(gsub("\\*", "_", gsub(":", "_", allele)), ".fasta"))
# Write GTF
write.table(full_join(forGtf_CDS, forGtf_exon), file = paste0(gsub("\\*", "_", gsub(":", "_", allele)), ".gtf"),
            sep = "\t", quote = F, row.names = F, col.names = F)
# Write upstream fasta
write.fasta(as.character(flanks$up_flank), 
            gsub("\\*", "_", gsub(":", "_", allele)),
            paste0("upstream_", gsub("\\*", "_", gsub(":", "_", allele)), ".fasta"))
# Write downstream fasta
write.fasta(as.character(flanks$down_flank), 
            gsub("\\*", "_", gsub(":", "_", allele)),
            paste0("downstream_", gsub("\\*", "_", gsub(":", "_", allele)), ".fasta"))
})


## If gDNA, Field matched, or matches REF ##
tmp <- chosen %>% mutate(row = 1:n()) %>% filter(status != "cDNA")
lapply(tmp$row, function(ix) {
  gene = chosen[ix,'gene']
  gene
  allele = ifelse(chosen[ix, 'status'] %in% c('gDNA','REF'), chosen[ix, 'haplotype'], chosen[ix, 'status'])
  allele
  ref <- reference_hlatypes %>% filter(hgnc_symbol == gsub("(\\bD\\w.+\\d*)", "HLA-\\1", gsub("(\\b\\w\\b)", "HLA-\\1", gene)))
  strand <- unique(filter(transcript_seq, hgnc_symbol == gsub("(\\bD\\w.+\\d*)", "HLA-\\1", gsub("(\\b\\w\\b)", "HLA-\\1", gene)))$strand)
  strand
  gen <- gen_msf[[gene]]
  nuc <- nuc_msf[[gene]]
  
  this_gen <- gen$seq[which(gen$nam == allele)] %>% gsub("-", "", .) %>% toupper
  this_nuc <- nuc$seq[which(nuc$nam == allele)] %>% gsub("-", "", .) %>% toupper
  ref_nuc <- nuc$seq[which(nuc$nam == ref$allele)] %>% gsub("-", "", .) %>% toupper
  ref_gen <- gen$seq[which(gen$nam == ref$allele)] %>% gsub("-", "", .) %>% toupper
  pa_nuc <- pairwiseAlignment(pattern = this_nuc, subject = this_gen, type = "local-global")
  ### Get location of CDS in gDNA
  range_loc <- pa_nuc@subject@range
  forGtf_CDS <- data.frame(seqname = "6",
                           source = gsub("\\*", "_", gsub(":", "_", allele)),
                           feature = "CDS",
                           start = ifelse(strand == 1, range_loc@start, nchar(this_gen) - (range_loc@start + range_loc@width - 1) + 1),
                           end = ifelse(strand == 1, range_loc@start + range_loc@width - 1, nchar(this_gen) - range_loc@start +1),
                           score = ".",
                           strand = ifelse(strand == 1, "+", "-"),
                           frame = "0",
                           attribute = paste0('gene_id "', gsub("\\*", "_", gsub(":", "_", allele)),
                                              '"; transcript_id "', gsub("\\*", "_", gsub(":", "_", allele)),
                                              '";'))
  ### Get exon locations
  cds_in_gDNA <- IRanges(start = range_loc@start, 
                         end = range_loc@start + range_loc@width - 1)
  # Intron positions
  introns <- as.data.frame(deletion(pa_nuc)) %>% 
    mutate(cumwidth = cumsum(width),
           newstart = start + lag(cumwidth, default = 0),
           newend = newstart + width - 1)
  introns_in_gDNA <- IRanges(start = introns$newstart + range_loc@start,
                             end = introns$newend + range_loc@start)
  # Exon locations
  exon_locs <- setdiff.Vector(cds_in_gDNA, introns_in_gDNA)
  if(strand == 1) {
    forGtf_exon <- data.frame(seqname = "6",
                              source = gsub("\\*", "_", gsub(":", "_", allele)),
                              feature = "exon",
                              start = exon_locs@start,
                              end = exon_locs@start + exon_locs@width - 1,
                              score = ".",
                              strand = ifelse(strand == 1, "+", "-"),
                              frame = "0",
                              attribute = paste0('gene_id "', gsub("\\*", "_", gsub(":", "_", allele)),
                                                 '"; transcript_id "', gsub("\\*", "_", gsub(":", "_", allele)),
                                                 '"; exon_id "', gsub("\\*", "_", gsub(":", "_", allele)), ".", 1:length(exon_locs), '";'))
  } else {
    forGtf_exon <- data.frame(seqname = "6",
                              source = gsub("\\*", "_", gsub(":", "_", allele)),
                              feature = "exon",
                              start = nchar(this_gen) - (exon_locs@start + exon_locs@width - 1) + 1,
                              end = nchar(this_gen) - exon_locs@start + 1,
                              score = ".",
                              strand = "-",
                              frame = "0",
                              attribute = paste0('gene_id "', gsub("\\*", "_", gsub(":", "_", allele)),
                                                 '"; transcript_id "', gsub("\\*", "_", gsub(":", "_", allele)),
                                                 '"; exon_id "', gsub("\\*", "_", gsub(":", "_", allele)), ".", 1:length(exon_locs), '";'))
    
  }
  
  ### Get Up/Downstream sequences for allele
  trans <- ref %>% inner_join(summ_IMGT_to_ref) %>% inner_join(transcripts)
  trans_pa <- apply(trans, 1, function(tx){
    if(strand == 1){
      gen_context_pa <- pairwiseAlignment(pattern = ref_gen, # Need to use location of ref allele gDNA
                                          subject = as.character(tx['context_seq']),
                                          type = "global-local")
    } else {
      gen_context_pa <- pairwiseAlignment(pattern = as.character(reverseComplement(DNAString(ref_gen))), # Need to use location of ref allele gDNA
                                          subject = as.character(reverseComplement(DNAString(tx['context_seq']))),
                                          type = "global-local")
    }
    return(gen_context_pa)
  })
  
  trans$ref_context_start <- sapply(trans_pa, function(i){ i@subject@range@start })
  trans$ref_context_end <- sapply(trans_pa, function(i){ i@subject@range@start }) + sapply(trans_pa, function(i){ i@subject@range@width })
  flanks <- apply(trans, 1, function(tx){
    start_in_context <- as.numeric(tx['ref_context_start'])
    end_in_context <- as.numeric(tx['ref_context_end'])
    if(strand == 1) {
      upstream_100 <- substr(tx['context_seq'], start_in_context-100, start_in_context-1)
      downstream_100 <- substr(tx['context_seq'], end_in_context, end_in_context+100-1)
    } else {
      upstream_100 <- substr(as.character(reverseComplement(DNAString(tx['context_seq']))), start_in_context-100, start_in_context-1)
      downstream_100 <- substr(as.character(reverseComplement(DNAString(tx['context_seq']))), end_in_context, end_in_context+100-1)
    }
    data.frame(up_flank = upstream_100, down_flank = downstream_100)
  }) %>% do.call(rbind, .) %>% filter(nchar(up_flank) == 100 & nchar(down_flank)==100) %>% unique
  
  # Write in FASTA
  write.fasta(ifelse(strand == 1, 
                     list(this_gen), 
                     list(as.character(reverseComplement(DNAString(this_gen))))), 
              gsub("\\*", "_", gsub(":", "_", allele)),
              paste0(gsub("\\*", "_", gsub(":", "_", chosen[ix, 'haplotype'])), ".fasta"))
#   # Write GTF
#   write.table(full_join(forGtf_CDS, forGtf_exon), 
#               file = paste0(gsub("\\*", "_", gsub(":", "_", chosen[ix, 'haplotype'])), ".gtf"),
#               sep = "\t", quote = F, row.names = F, col.names = F)
#   # Write upstream fasta
#   write.fasta(as.character(flanks$up_flank), 
#               gsub("\\*", "_", gsub(":", "_", allele)),
#               paste0("upstream_", gsub("\\*", "_", gsub(":", "_", chosen[ix, 'haplotype'])), ".fasta"))
#   # Write downstream fasta
#   write.fasta(as.character(flanks$down_flank), 
#               gsub("\\*", "_", gsub(":", "_", allele)),
#               paste0("downstream_", gsub("\\*", "_", gsub(":", "_", chosen[ix, 'haplotype'])), ".fasta"))
})


#####
a <- which(gen_msf$A$nam == "A*02:01:01:01")
b <- which(gen_msf$A$nam == "A*01:01:01:01")
this_msa <- msa(c(toupper(as.character(gen_msf$A$seq[a])), toupper(as.character(gen_msf$A$seq[b]))), type = "dna")
print(this_msa, show = "complete")
correct_sequence <- data.frame(a = toupper(as.character(gen_msf$A$seq[a])),
                               b = toupper(as.character(gen_msf$A$seq[b])),
                               cons = msaConsensusSequence(this_msa)) %>%
  gather(key, value) %>% group_by(key) %>% mutate(value = strsplit(value, split = "")) %>%
  unnest(value) %>% mutate(i = 1:n()) %>% spread(key, value)
correct_sequence %>% 
  filter(! (a == "-" & b == "-") & cons == "?") %>% nrow



##### END #####
