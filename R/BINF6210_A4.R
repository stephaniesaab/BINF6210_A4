# Assignment 4: BINF6210: Software for Bioinformatics
# Date: December 5, 2025
# Author: Stephanie Saab
# Taxonomy of interest: Mus musculus
# Research Questions: Gene classification fo mitochondrial DNA in mice

# Software Tools - Anopheles Classification Activity ----

##***************************
## Software Tools - Anopheles Classification Activity with Random Forest
## Script: anopheles_scrambled.R
##
## Student: Stephanie Saab
## Author: Karl Cottenie
## Updated: October 23rd, 2025 by Stephanie Saab
## Description: This script fetches Anopheles COI gene sequences from GenBank,
## and builds a Random Forest classifier to predict species based on their 
## dinucleotide frequencies. Data are split into training (80%) and validation
## (20%) datasets to assess the model's generalization.
##***************************

# PART 1: LOAD PACKAGES & SETUP ----

# Loop checks if the packages are already installed, if not it installs them, if they are it loads them.
required_packages <- c("DECIPHER", "BiocManager", "conflicted", "rentrez", "seqinr", "Biostrings", "randomForest", 
                       "tidyverse", "dplyr", "vegan", "countrycode", "viridis", 
                       "ggplot2", "igraph")
#Install pwalign for version of R if not available
BiocManager::install("pwalign")
for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package)
  } else {
    library(package, character.only = TRUE)
  }
}

# Adjust for conflicts
conflict_prefer("filter", "dplyr")

# PART 2: DOWNLOAD & IMPORT DATA ----
# Ensure working directory is set to source file location (".../BINF6210_A4/R")

#Provenance (record in README or variables)
# Get the GenBank search hits based on the search terms of interest (database = nucleotide, organism = mus musculus, gene = COI and cytb)
df_mus_coi_search <- entrez_search(db = "nuccore", term = "mus musculus [ORGN] AND COI [gene]", use_history = T)
df_mus_coi_search # Check that it resulted in 136 hits, 20 IDs and a web_history object
df_mus_cytb_search <- entrez_search(db = "nuccore", term = "mus musculus [ORGN] AND CYTB [gene]", use_history = T)
df_mus_cytb_search # Check that it resulted in 2455 hits, 20 IDs and a web_history object

#Downloaded from: https://www.ncbi.nlm.nih.gov/nuccore/
#Search terms: (Mus musculus[Organism] AND COI[Gene]) and (Mus musculus[Organism] AND cytb[Gene])
# Read in the sequence data as a DNA StringSet from data GenBank database
st_coi <- readDNAStringSet("../Data/COI_mouse.fasta")
st_coi # Check that it is a DNAStringSet object of length 136
st_cytb <- readDNAStringSet("../Data/CYTB_mouse.fasta")
st_cytb # Check that it is a DNAStringSet object of length 2455

#Download reference sequences for COI and CYTB in mice from NCBI
# Link for COI: https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=17708#:~:text=Official%20Symbol%20mt%2DCo1provided,.p6%20(GCF_000001635.26)
# Link for CYTB: https://www.ncbi.nlm.nih.gov/datasets/gene/17711/
#Write sequences as fasta file in data folder
coi_acc <- "NC_005089.1" #Accession for COI
coi_seq <- entrez_fetch(db = "nuccore", 
                        id = coi_acc, 
                        rettype = "fasta", 
                        retmode = "text",
                        seq_start = 5328,
                        seq_stop = 6872)
write(coi_seq, file = "../Data/COI_ref.fasta")
cytb_acc <- "NC_005089.1"
cytb_seq <- entrez_fetch(db = "nuccore", 
                         id = cytb_acc, 
                         rettype = "fasta", 
                         retmode = "text",
                         seq_start = 14145,
                         seq_stop = 15288)
write(cytb_seq, file = "../Data/CYTB_ref.fasta") #Write seque ce as fasta file in data folder

# PART 3: PRE-PROCESSING & QC ----

#Clean the sequences to remove non-base characters (non ACGTN)
clean_dna <- function(dna_str) {
  #Make all bases uppercase
  seq <- toupper(as.character(dna_str))

  #Replace all ambiguity codes with N, and remove gaps
  seq <- gsub("-", "", seq)
  seq <- gsub("[^ACGTN]", "N", seq)
  
  # Make DNAStringSet and preserve original names
  dna_clean <- DNAStringSet(seq)
  names(dna_clean) <- names(dna_str)
  
  return(dna_clean)
}

##Clean the fasta, remove headers and formatting characters
clean_fasta <- function(fasta) {
  
  #Read file
  raw <- readDNAStringSet(fasta)
  
  #Remove header and non-base characters (formatting characters)
  cleaned_seq <- lapply(as.character(raw), function(x) {
      seq <- clean_dna(x)
    
    #Remove empty of NA sequences
    if (nchar(seq) == 0 || is.na(seq)) 
      return(NULL)
    
    return(DNAString(seq))
  })
  
  #Remove nulls
  cleaned_seq <- cleaned_seq[!sapply(cleaned_seq, is.null)]
  
  #Make DNAStringSet
  cleaned_seq <- DNAStringSet(cleaned_seq)
  
  #Keep original names
  names(cleaned_seq) <- names(raw)[!sapply(cleaned_seq, is.null)]
  return(cleaned_seq)
}

#Make clean reference sequences as DNAStringSet objects
coi_ref <- clean_fasta("../Data/COI_ref.fasta")
cytb_ref <- clean_fasta("../Data/CYTB_ref.fasta")

#Make clean genome sequences
coi_clean <- clean_dna(st_coi)
cytb_clean <- clean_dna(st_cytb)

#Shorten sequences to only the gene of interest, only do this for sequences that are very long to avoid going through all sequences and speed the process
#Lengths based on the lengths of the genes
cytb_long <- DNAStringSet(cytb_clean[ Biostrings::width(cytb_clean) > 15000])
coi_long <- DNAStringSet(coi_clean[ Biostrings::width(coi_clean) > 15000])
cytb_short <- DNAStringSet(cytb_clean[ Biostrings::width(cytb_clean) < 15000])
coi_short <- DNAStringSet(coi_clean[ Biostrings::width(coi_clean) < 15000])

#Check there are some full mtDNA genomes
length(cytb_long) #Should return 382
length(coi_long) #Should return 15
length(cytb_short) #Should return 2073
length(coi_short) #Should return 136

extract_gene <- function(sequences, 
                         start_pos,
                         end_pos, 
                         ref,
                         min_score = 200) {
  #sequences: DNAStringSet of genomes
  #Start_pos, end_pos: gene coordinates for extraction
  #ref_seq: DNAString of reference gene
  #min_score: for alignment, default is 200
  
  dna_mat <- pwalign::nucleotideSubstitutionMatrix(match = 2, 
                                                   mismatch = -1, 
                                                   baseOnly = FALSE)
  extracted <- lapply(sequences, function(genome) {
    #Extract window first to make it faster
    start_use = max(1, start_pos)
    end_use = min(end_pos, length(genome))
    genome <- subseq(genome, start = start_use, end = end_use)
  
    #Alignment-based extraction
    aln <- pwalign::pairwiseAlignment(pattern = ref, 
                                      subject = genome,
                                      type = "local", 
                                      substitutionMatrix = dna_mat)
    if (score(aln) < min_score) return(NULL)
  
    subseq(genome, start(subject(aln)), end(subject(aln)))
  })
  #Remove nulls from failed alignments
  extracted <- extracted[!sapply(extracted, is.null)]
  
  #Convert to DNAStringSet and keep names
  extracted <- DNAStringSet(extracted)
  names(extracted) <- names(sequences)[!sapply(extracted, is.null)]
  return(extracted)
}

#CYTB Coordinate in Mus musculus reference mtDNA based on NCBI reference sequence
cytb_start <- 14000
cytb_end   <- 15400
coi_start <- 5200
coi_end <- 7000
cytb_extracted <- extract_gene(cytb_long,
                               cytb_start,
                               cytb_end, 
                               cytb_ref[[1]],)
coi_extracted <- extract_gene(coi_long,
                               coi_start,
                               coi_end, 
                               coi_ref[[1]],)

#Convert to tidy dataframes as tibbles
df_coi <- data.frame(id = names(coi_clean), sequence = as.character(coi_extracted), label = "COI", stringAsFactors = FALSE)
df_cytb <- tibble(id = names(cytb_clean), sequence = as.character(cytb_clean), label = "cytb", stringAsFactors = FALSE)

#Concatenate the genes with the extracted genes
combined_cytb <- c(cytb_short, cytb_extracted)
combined_coi <- c(coi_short,coi_extracted)

#Check lengths match up with original
length(combined_coi) #Should return 136
length(combined_cytb) #Should return 2455

#Convert to tidy dataframes and add a lengths column to check sequence lengths
df_coi <- data.frame(
  id = names(combined_coi), 
  sequence = as.character(combined_coi),
  length = width(combined_coi),
  label = "COI", stringAsFactors = FALSE)
df_cytb <- data.frame(
  id = names(combined_cytb),
  sequence = as.character(combined_cytb),
  length = width(combined_cytb),
  label = "cytb", 
  stringAsFactors = FALSE)




# Part 6 BTS ====
# 
# extract_gene_by_pos <- DNAStringSet(lapply(cytb_clean, function(genome) {
#   subseq(genome, start = cytb_start, end = cytb_end)
# }))
# #Extract gene sequences from whole genomes (local alignment)
# dna_mat <- pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
# min_score <- 200  # adjust based on your gene length
# 
# extract_gene <- function(genome, ref) {
#   aln <- pwalign::pairwiseAlignment(pattern = ref,
#                                     subject = genome,
#                                     type = "local",
#                                     substitutionMatrix = dna_mat)
#   if(score(aln) < min_score) return(NA)
#   start = start(subject(aln))
#   end = end(subject(aln))
#   subseq(genome, start, end)
# }
#   
# cytb_extracted <- lapply(cytb_long, function(genome) {
#   extract_gene(genome, cytb_ref[[1]])
# })
# 
# cytb_extracted_set <- DNAStringSet(cytb_extracted)
#   ref_seq <- ref[[1]]
#   seq <- genome[[1]]
#   aln <- DECIPHER::AlignSeqs(DNAStringSet(c(ref_seq, seq)))
#   
#   aln_ref <- as.character(aln$AlignedSeqs[[1]])
#   aln_seq <- as.character(aln$AlignedSeqs[[2]])
#   
#   #Get location for where reference has bases, not gaps
#   ref_positions <- which(strsplit(ref_aln_char, "")[[1]] != "-")
#   
#   #Extract corresponding position from aligned subject
#   seq_block <- strsplit(seq_aln_char, "")[[1]][ref_positions]
#   seq_block <- paste(seq_block, collapse = "")
#    DNAString(seq_block)
# }
# cytb_extracted <- lapply(cytb_long, extract_gene, ref = cytb_ref)
# 
#   
  # aln <- pwalign::pairwiseAlignment(pattern = ref_seq,
  #                                   subject = genome,
  #                                   type = "local", 
  #                                   substitutionMatrix = NULL)
  # 
  # start <- start(subject(aln))
  # end <- end(subject(aln))
  # 
  # return(subseq(seq, start,))
  # dna_mat <- pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
  # align <- pwalign::pairwiseAlignment(ref, genome, type = "local", substitutionMatrix = dna_mat)
  # if (score(align) < min_score) {
  #   return(NA)
  # }
  # 
  # rng <- subject(align)
  # subseq(genome, start = start(rng), end = end(rng))
  # 


#Shorten sequences to only the gene of interest, only do this for sequences that are very long to avoid going through all sequences and speed the process
#Lengths based on the lengths of the genes
# coi_extracted <- lapply(coi_long, extract_gene, ref = coi_ref)
# cytb_extracted <- lapply(cytb_long, extract_gene, ref = cytb_ref)

#Convert to tidy dataframes as tibbles
df_coi <- data.frame(id = names(coi_extracted), sequence = as.character(coi_extracted), label = "COI", stringAsFactors = FALSE)
df_cytb <- tibble(id = names(cytb_extracted), sequence = as.character(cytb_extracted), label = "cytb", stringAsFactors = FALSE)

#Default parameters (can change)
k_mer = 3 #K-mer threshold, can change with speed and sample size
n_count = 0.05 #Max allowed N's per sequence (percent)
set.seed(916)

#Remove sequences with non-ACGT characters
clean_sequence <- function(df, n_count, gene_length) {
  # Clean the sequences so that it removes non-AGCT characters
  # Removes sequences with too many N's based on n_count
  # Remove sequences out of length range based on std_length
  # Remove sequence duplicates
  
  #Make column with information
  df <- df %>% 
    mutate(len = nchar(sequence), 
           n_counts = str_count(sequence, "[Nn]"),
           perc_n = n_counts / len,
           non_ACGT = str_detect(sequence, "[^ACGTacgtNn]"))
   
  #Apply thresholds for quality control
  df_clean <- df %>% 
    filter(!non_ACGT) %>% 
    filter(perc_n < n_count) %>% 
    distinct(sequence, .keep_all = TRUE) %>% 
    filter(nchar(sequence) < gene_length)
    
  return(df_clean)
}

df_coi_clean <- clean_sequence(df_coi, n_count, 800) #Gene length of COI ~ 700bp
df_cytb_clean <- clean_sequence(df_cytb, n_count, 1500 ) #Gene length of cytb ~ 1200bp

#Checkpoint for cleaning and QC
nrow(df_coi_clean) #Check that it is 52
dim(df_coi_clean) #Check that it returns 52 x 7
nrow(df_cytb_clean) #Check that it returns 943
dim(df_cytb_clean) #Check that it returns 943 x 7



# PART 4: PRE-PROCESSING ----

# Feature extraction
#Extract kmer frequencies (start with 4)
#Add GC% as features
#Convert raw counts -> Relative frequencies

#Feature QC / filtering
#remove kmer features with near-zero variance (caret::nearZeroVar) or very low mean across sequences

#Balance classes -> downsample majority class to match minority (simple and robust) OR upsample or use class weights

#Train / test / validation split
#stratified split using caret::createDataPartition (70% train, 15% test, 15% validation)

#Diagnostics prior to modeling 
#PCA on filtered features --> plot colored by label
#Correlation heatmap of top features

#save intermediate outputs


# Visualizations
# 1 = histogram of sequence length per label (QC / diagnostics)
#Plot length distributions




# 2 = bar plot of counts per label (before / after balancing)
# 3 = table = number sequencse per label after QC
