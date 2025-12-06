# Software Tools - Assignment 4 ----
## Script: BINF6210_A4.R
## Date: December 5, 2025
## Taxonomy of interest: Mus musculus
## Student: Stephanie Saab
## Research Questions: Gene classification fo mitochondrial DNA in mice. It builds a Random Forest classifier to predict mitochondrial genes (COI vs CYTB) based on their trinucleotide frequencies. Data are split into training (70%) and validation (30%) datasets to assess the model's generalization.
##***************************

# PART 1: LOAD PACKAGES & SETUP ----

# Loop checks if the packages are already installed, if not it installs them, if they are it loads them.
required_packages <- c("gridExtra", "BiocManager", "conflicted", "rentrez", "randomForest", "tidyverse", "dplyr", "viridis", "ggplot2")

for (package in required_packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package)
  } else {
    library(package, character.only = TRUE)
  }
}

#If these packages are not installed they must be installed with BiocManager (e.g. BiocManager::install("pwalign"))
library("pwalign")
library("Biostrings")
library("DECIPHER")
# Adjust for conflicts
conflict_prefer("filter", "dplyr")

# PART 2: DOWNLOAD & IMPORT DATA ----
# Ensure working directory is set to source file location (".../BINF6210_A4/R")

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
write(cytb_seq, file = "../Data/CYTB_ref.fasta")

#Inspecting the data
max(width(st_coi)) #Returns 16302
max(width(st_cytb)) #Returns 16379
min(width(st_coi)) #Returns 132
min(width(st_cytb)) #Returns 193
median(width(st_coi)) #Returns 653
median(width(st_cytb)) #Returns 1140

# PART 3: PRE-PROCESSING & QC ----

#Clean the sequences to remove non-base characters (non ACGTN)
clean_dna <- function(dna_str) {
  #Make all bases uppercase
  seq <- toupper(as.character(dna_str))

  #Replace all ambiguity codes with N, and remove gaps
  #These cannot be interpreted later by classifiers and alignment
  seq <- gsub("-", "", seq)
  seq <- gsub("[^ACGTN]", "N", seq)
  
  return(seq)
}

##Clean the fasta, remove headers and formatting characters
clean_fasta <- function(fasta) {
  
  #Read file
  raw <- readDNAStringSet(fasta)
  
  #Remove header and non-base characters (formatting characters)
  cleaned_seq <- lapply(seq_along(raw), function(x) {
      seq <- clean_dna(raw[[x]])
    
    #Remove empty of NA sequences
    if (nchar(seq) == 0 || is.na(seq)) 
      return(NULL)
    dna <- DNAString(seq)
    return(dna)
  })
  
  #Remove nulls
  cleaned_seq <- cleaned_seq[!sapply(cleaned_seq, is.null)]

  #Make DNAStringSet
  cleaned_seq <- do.call(DNAStringSet, cleaned_seq)
  
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
length(coi_short) #Should return 121

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
df_coi <- data.frame(id = names(coi_extracted), sequence = as.character(coi_extracted), label = "COI", stringAsFactors = FALSE)
df_cytb <- tibble(id = names(cytb_extracted), sequence = as.character(cytb_extracted), label = "cytb", stringAsFactors = FALSE)

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
  label = "CYTB", 
  stringAsFactors = FALSE)

#Balance the datasets with Random Undersampling
set.seed(123)

# Under-sample CYTB (2455) to match COI size (136)
df_cytb <- df_cytb %>% 
  slice_sample(n = nrow(df_coi))

# Combine them into one dataframe
df_balanced <- bind_rows(df_coi, df_cytb)

# Shuffle the rows to sample the training, testing, validation splits
df_balanced <- df_balanced %>% sample_frac(1)

# Check that the sampling is balanced, count should be the same for each gene (expect 136)
df_balanced %>%
  group_by(label) %>%
  dplyr::count()


# Part 4: Feature selection: Kmers, GC, length, AT/GC skew ====

#Sequence length was already added during pre-processing to check sequence lengths after gene extraction
#Convert sequences to DNAStringSet object 
df_balanced$sequence2 <- DNAStringSet(df_balanced$sequence)
#Check it's class DNAStringSet: 
class(df_balanced$sequence2)

#Calculating the nucleotide frequencies and appending onto our dataframe using cbind(). Using letterFrequency instead of a
df_balanced <- cbind(df_balanced, as.data.frame(letterFrequency(df_balanced$sequence2, letters = c("A", "C", "G", "T", "N"))))

#Check the new columns were added
names(df_balanced) #Should return added rows for nucleotides

#Adding A, T, and G proportions in relation to total nucleotides
df_balanced$Aprop <- (df_balanced$A) / (df_balanced$A + df_balanced$T + df_balanced$C + df_balanced$G)

df_balanced$Tprop <- (df_balanced$T) / (df_balanced$A + df_balanced$T + df_balanced$C + df_balanced$G)

df_balanced$Gprop <- (df_balanced$G) / (df_balanced$A + df_balanced$T + df_balanced$C + df_balanced$G)

#Adding multi-nucleotide frequencies (k-mers of length 3) (64 features)
df_balanced <- cbind(df_balanced, as.data.frame(trinucleotideFrequency(df_balanced$sequence2, as.prob = TRUE)))

#Adding GC content (%)
df_balanced$GC_cont <- (str_count(df_balanced$sequence, "G") + 
                           str_count(df_balanced$sequence, "C")) / df_balanced$length

#Adding AT/GC Skew to capture directional bias per gene
df_balanced$AT_skew <- (str_count(df_balanced$sequence, "A") - str_count(df_balanced$sequence, "T")) /
                          (str_count(df_balanced$sequence, "A") + str_count(df_balanced$sequence, "T"))
df_balanced$GC_skew <- (str_count(df_balanced$sequence, "G") - str_count(df_balanced$sequence, "C")) / 
  (str_count(df_balanced$sequence, "G") + str_count(df_balanced$sequence, "C"))

# Part 5: TRAINING THE CLASSIFICATION MODEL ====

#Split into training + validation sets
#Set seed because this step uses randomization
set.seed(217)

#Convert sequences back to string to use dataframe as a tibble
df_balanced$sequence2 <- as.character(df_balanced$sequence2)

#Make validation dataset, 30% of samples
df_validation <- df_balanced %>% 
  group_by(label) %>% 
  sample_n(0.3 * nrow(df_balanced)) %>% 
  ungroup()

#Training dataset is rest of data (70%)
df_training <- df_balanced %>% 
  filter(!id %in% df_validation$id)

#Check sample size for each marker is same
table(df_validation$label) #Should be 54 each
table(df_training$label) #Should be 82 each

#Combine labels and features into ML tables
feature_cols <- c("length", "Aprop", "Tprop", "Gprop", "GC_cont", "AT_skew", "GC_skew")
kmer_cols <- names(df_balanced)[grep("^[AGCT]{3}$", names(df_balanced))]
all_features <- c(feature_cols, kmer_cols)

#Make ML tables
df_training_ml <- df_training[, c("label", all_features)]
df_validation_ml <- df_validation[, c("label", all_features)]

#Make sure the label column is a factor
df_training_ml$label <- base::as.factor(df_training_ml$label)
df_validation_ml$label <- base::as.factor(df_validation_ml$label)

#Train Random Forest on all featuers
set.seed(123)
rf_all <- randomForest(
  x = df_training_ml[, all_features],
  y = df_training_ml$label,
  ntree = 200,
  importance = TRUE
)

#View classifier summary
rf_all
rf_all$importance
rf_all$confusion

#Predict on validation set
predict_validation <- predict(rf_all, df_validation_ml[, all_features])

#Confusion matrix
table(observed = df_validation_ml$label, predicted = predict_validation)
# Part 6: VISUALIZATIONS ====
# Visualization 1: Histograms of sequence lengths
df_before <- rbind(
  data.frame(label = "COI", length = width(st_coi)),
  data.frame(label = "CYTB", length = width(st_cytb))
)

df_after_coi <- df_balanced %>% 
  filter(label == "COI") %>% 
  select(label, length)
df_after_cytb <- df_balanced %>% 
  filter(label == "CYTB") %>% 
  select(label, length)
df_after <- rbind(df_after_coi, df_after_cytb)
hist_plot1 <- ggplot(df_before, aes(x = length, fill = label)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 40) +
  scale_fill_viridis_d() +
  scale_x_continuous(trans = "log10")+
  theme_classic() +
  labs(
    title = "Sequence Length Distribution \n Before Cleaning & Extraction",
    x = "Sequence Length (bp)",
    y = "Count",
  )

hist_plot2 <-ggplot(df_after, aes(x = length, fill = label)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 40) +
  scale_fill_viridis_d() +
  theme_classic() +
  labs(
    title = "Sequence Length Distribution \n After Cleaning & Extraction",
    x = "Sequence Length (bp)",
    y = "Count",
  )
hist_seq_lengths <- grid.arrange(hist_plot1, hist_plot2, ncol=2)

#Visualization 2: GC content as violin plots
violin_GC <- ggplot(df_balanced, aes(x = label, y = GC_cont, fill = label))+
  geom_violin(trim = FALSE, alpha = 0.7)+
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.4)+
  scale_fill_viridis_d()+
  theme_classic()+
  labs(
    title = "GC(%) Content Distribution for COI and CYTB mitochondrial \n genes",
    x = "Gene",
    y = "GC Content (%)"
  )

#Visualization 3: PCA of Kmers
kmers_only <- df_balanced[, kmer_cols]
pca <- prcomp(kmers_only, scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  label = df_balanced$label
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = label)) + 
  geom_point(size = 2, alpha = 0.7)+
  scale_color_viridis_d()+
  theme_classic()+
  labs(title = "PCA of Trinucleotide Frequencies (k=3) in mitochondrial \n genes")

#Visualization 4: Random Forest Feature Importance plot

#Get importance of each feature
rf_imp_df <- as.data.frame(importance(rf_all))
rf_imp_df$feature <- rownames(rf_imp_df)

#Sort by MeanDecreaseGini (impurity)
#Keep only top 20 by Gini
top20 <- rf_imp_df %>% 
  select(feature, MeanDecreaseAccuracy, MeanDecreaseGini) %>% 
  pivot_longer(cols = c(MeanDecreaseAccuracy, MeanDecreaseGini), 
               names_to = "Metric", 
               values_to = "Importance") %>% 
  filter(Metric == "MeanDecreaseGini") %>% 
  arrange(desc(Importance)) %>% 
  dplyr::slice(1:20) %>% 
  pull(feature)

#Build the dataframe to plot
plot_rf_df <- rf_imp_df %>% 
  pivot_longer(cols = c(MeanDecreaseAccuracy, MeanDecreaseGini), 
               names_to = "Metric", 
               values_to = "Importance") %>% 
  filter(feature %in% top20)
#Use ggplot over varImpPlot for control over aesthetics
rf_feature_plot <- ggplot(plot_rf_df, aes(x = reorder(feature, Importance), 
             y = Importance, fill = Metric))+
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_viridis_d() +
  theme_classic(base_size = 14) +
  labs(
    title = "Top 20 Feature Importance (Gini & Accuracy) in \n mitochondrial gene classification",
    x = "Feature", 
    y = "Importance Score", 
    fill = "Metric"
  )

#Table of Classification results
class_table <- table(observed = df_validation_ml$label, predicted = predict_validation)

#Save the plots to the figures folder
plot_list <- list(hist_seq_lengths = hist_seq_lengths, violin_GC = violin_GC, pca_plot = pca_plot, rf_feature_plot = rf_feature_plot)
for (nm in names(plot_list)) {
  ggsave(
    filename =  paste0("../Figs/", nm, ".png"),
    plot = plot_list[[nm]],
    width = 6,
    height = 4,
    dpi = 500
  )
  print(c("Saved",nm))
}
  
#Session info
sessionInfo()
