# ===============================================
# NX and Tau Calculation Pipeline with Batch Correction
# ===============================================
# Author: Michael Grapin
# Purpose: Compute normalized expression (NX) and tissue specificity (Tau)
#          for samples sequenced across multiple runs and tissues
# Following Protocol From: https://www.science.org/doi/10.1126/science.aax9198
# Protocol Github Repo: https://github.com/human-protein-atlas/BloodAtlas
# ===============================================
wkdir = "~/Downloads/TOAD_ANALYSIS/MultiTissue_RNAseq/Tissue_Expression"
setwd(wkdir)
# ------------------------------
# Load required libraries
# ------------------------------
library(edgeR)       # For TMM normalization
library(limma)       # For batch correction
library(tibble)      # For data frame manipulation
library(dplyr)       # For data wrangling
library(readr)       # For reading CSV files
library(tidyr)       # For unite() function
library(stringr)     # For str_extract() function 
library(matrixStats) # For matrix classifcation function


# Data prepocessing Steps

# ------------------------------
# Input files
# ------------------------------
TPM_file <- "reduc_gene_tpm_matrix.csv" 
design_file <- "../DEseq/design.csv"
outdir <- "Altas_Results"

# ------------------------------
# Step 0: Read in TPM and design files
# ------------------------------
TPM_data <- read_csv(TPM_file)
design <- read_csv(design_file)

# ------------------------------
# Step 0a: Prepare TPM matrix
# ------------------------------
TPM_data <- as.data.frame(TPM_data)
rownames(TPM_data) <- TPM_data[,1]  # Use first column as rownames
TPM_data <- TPM_data[,-1]           # Remove first column

# ------------------------------
# Step 0b: Prepare design file
# ------------------------------
design <- as.data.frame(design)
rownames(design) <- design$sample_name

# Keep only relevant columns
# Ensure design$individual and design$tissue are also defined
design$condition <- factor(unite(design, col = "condition", tissue, timepoint, sep = "_", remove = FALSE)$condition)
design$timepoint <- factor(design$timepoint, levels = c("con", "0h", "3h", "6h", "9h", "12h"))
design$tissue <- factor(design$tissue, levels = c("muscle", "liver", "parotoid"))
design <- design %>% mutate(individual = str_extract(sample_name, "[A-Z]{3}[0-9]{3}"))

# Temporary remove TJF235 
# design  <- design %>% filter(sample_name != "TJF235M_con")
# TPM_data <- TPM_data[, colnames(TPM_data) != "TJF235M_con"]


# Reorder TPM columns to match design rows
ordered <- match(rownames(design), colnames(TPM_data))
TPM_data_ordered <- TPM_data[, ordered]
###########################################################################################################
###########################################################################################################
# ------------------------------
# Step 1: TMM normalization within each sequencing run
# ------------------------------
# Reasoning:
# Runs can have composition bias, so we normalize separately within each run
TPM_tmm <- TPM_data_ordered  # initialize normalized matrix
run <- design$Run             # sequencing run for each sample

for (r in unique(run)) {
  idx <- run == r
  dge <- DGEList(counts = TPM_data_ordered[, idx])
  dge <- calcNormFactors(dge, method = "TMM")
  # Apply TMM normalization
  TPM_tmm[, idx] <- sweep(TPM_data_ordered[, idx], 2, dge$samples$norm.factors, FUN = "/")
}


# ------------------------------
# Step 2: Pareto scaling (gene-wise within each run)
# ------------------------------
# Stabilizes variance and compresses high-expression outlines
# Protocol doesn't subtract from mean(x) like in traditional Pareto Scaling formula
pareto <- function(x) {
  s <- sd(x, na.rm = TRUE)
  out <- x / sqrt(s)
  out[!is.finite(out)] <- 0
  out
}

TPM_tmm_pareto <- TPM_tmm

for (r in unique(run)) {
  idx <- run == r
  TPM_tmm_pareto[, idx] <- t(apply(TPM_tmm[, idx], 1, pareto))
}


# check which SD are making NA, NaN, Inf

# CANT REMOVE BATCH EFFECTS PROCEED TO TISSUE LEVEL AVERAGE

# 
# # ------------------------------
# # Step 3: Remove sequencing run effects (batch correction)
# # ------------------------------
# # Batch = sequencing run (technical)
# # Design = individual (biological replicate, preserved)
# 
# # Probably should remove TJF235 because after batch correction scales aren't the same
# # I feel like the batch correction is just distorting the singal and that run is minor compared to our
# # biological signal. 
# 
# # Check the overlap between Run (Batch) and Tissue (Biology)
# counts_table <- table(design$Run, design$tissue)
# print(counts_table) # confounded experiment, cant do batch correction
# 
# # I think the solution is just to check the similarity before and after and if its minor just keep it in the back 
# # pocket and use the TPM -> TMM -> Pareto. 
# 
# # 1. Define the sample to exclude
# exclude_sample <- "TJF235M_con"
# 
# # 2. Filter the metadata (design file)
# design_filtered <- design[design$sample_name != exclude_sample, ]
# 
# 
# # 4. Create the design matrix from the CLEANED metadata
# # We use as.factor to ensure levels are dropped if they no longer exist
# design_bio <- model.matrix(~ 0 + as.factor(design_filtered$condition))
# # Rename columns for clarity (optional)
# colnames(design_bio) <- levels(as.factor(design_filtered$condition))
# 
# 
# # Sequencing run is the batch effect
# # Fucntion
# batch_correct_by_individual <- function(expr_matrix, design,
#                                         sample_col = "sample_name",
#                                         individual_col = "individual",
#                                         batch_col = "Run") {
#   
#   # -------------------------------
#   # 1. Remove genes with any NA
#   # -------------------------------
#   expr_matrix_clean <- expr_matrix[
#     apply(expr_matrix, 1, function(x) !any(is.na(x))), 
#   ]
#   
#   # -------------------------------
#   # 2. Count how many batches each individual appears in
#   # -------------------------------
#   batch_table <- design %>%
#     distinct(.data[[individual_col]], .data[[batch_col]]) %>%
#     count(.data[[individual_col]])
#   
#   multi_batch_individuals <- batch_table %>%
#     filter(n >= 2) %>%
#     pull(.data[[individual_col]])
#   
#   single_batch_individuals <- batch_table %>%
#     filter(n == 1) %>%
#     pull(.data[[individual_col]])
#   
#   # -------------------------------
#   # 3. Map individuals to samples
#   # -------------------------------
#   samples_multi <- design %>%
#     filter(.data[[individual_col]] %in% multi_batch_individuals) %>%
#     pull(.data[[sample_col]])
#   
#   samples_single <- design %>%
#     filter(.data[[individual_col]] %in% single_batch_individuals) %>%
#     pull(.data[[sample_col]])
#   
#   # -------------------------------
#   # 4. Subset expression matrix
#   # -------------------------------
#   expr_multi  <- expr_matrix_clean[, samples_multi, drop = FALSE]
#   expr_single <- expr_matrix_clean[, samples_single, drop = FALSE]
#   
#   # -------------------------------
#   # 5. Batch vector (must match column order)
#   # -------------------------------
#   batch_vector <- design %>%
#     filter(.data[[sample_col]] %in% colnames(expr_multi)) %>%
#     arrange(match(.data[[sample_col]], colnames(expr_multi))) %>%
#     pull(.data[[batch_col]])
#   
#   # -------------------------------
#   # 6. Log-transform
#   # -------------------------------
#   expr_multi_log <- log10(expr_multi + 1)
#   
#   # -------------------------------
#   # 7. Remove batch effects
#   # -------------------------------
#   expr_multi_corrected <- removeBatchEffect(
#     expr_multi_log,
#     batch = batch_vector, 
#     design = design_bio
#   )
#   
#   # -------------------------------
#   # 8. Inverse log transform
#   # -------------------------------
#   expr_multi_corrected <- 10^expr_multi_corrected - 1
#   expr_multi_corrected[expr_multi_corrected < 0] <- 0
#   
#   # -------------------------------
#   # 9. Recombine matrix
#   # -------------------------------
#   expr_corrected <- cbind(expr_single, expr_multi_corrected)
# 
#   # Restore original column order
#   expr_corrected <- expr_corrected[, colnames(expr_matrix_clean)]
#   
#   return(expr_corrected)
# }
# 
# NX <- batch_correct_by_individual(
#   expr_matrix = TPM_tmm_pareto,
#   design = design,
#   sample_col = "sample_name",
#   individual_col = "individual",
#   batch_col = "Run"
# )



# ------------------------------
# Step 4: Collapse replicates to tissue-level NX
# ------------------------------

# Future Note: Could summarize by tissue x timepoint (col == condition)

# Average across individuals for each tissue 
tissue <- design$tissue

tissue_avg <- sapply(unique(tissue), function(t) {
  rowMeans(TPM_tmm_pareto[, tissue == t, drop = FALSE])
})
# row = genes, col = tissues
colnames(tissue_avg) <- unique(tissue)

# ------------------------------
# Step 5: Filter low-expression genes
# ------------------------------
# NX < 1 treated as "not expressed" (removes zeros and noise)
tissue_avg_filtered <- tissue_avg
tissue_avg_filtered[tissue_avg_filtered < 1] <- NA


# Logical vector: TRUE if row has no NAs
complete_rows <- complete.cases(tissue_avg_filtered)

# Number of rows with values in all columns
sum(complete_rows)

# ------------------------------
# Adding Classification Bins
# ------------------------------
classify_gene <- function(x) {
  # x = numeric vector of NX values for one gene
  max_nx <- max(x, na.rm = TRUE)
  mean_nx <- mean(x, na.rm = TRUE)
  n <- length(x)
  
  # counts
  detected <- sum(x >= 1, na.rm = TRUE)
  above_25 <- sum(x >= 0.25 * max_nx, na.rm = TRUE)
  
  # tissue enriched
  if (sum(x >= 4 * sort(x, decreasing = TRUE)[2], na.rm = TRUE) == 1) {
    return("Tissue enriched")
  }
  
  # group enriched
  if (above_25 == 2) {
    top2 <- sort(x, decreasing = TRUE)[1:2]
    rest <- sort(x, decreasing = TRUE)[-(1:2)]
    if (mean(top2) >= 4 * max(rest, na.rm = TRUE)) {
      return("Group enriched")
    }
  }
  
  # tissue enhanced
  if (any(x >= 4 * mean_nx, na.rm = TRUE)) {
    return("Tissue enhanced")
  }
  
  # low tissue specificity
  if (detected >= 1) {
    return("Low tissue specificity")
  }
  
  return("Not detected")
}



detect_class <- function(x) {
  n <- length(x)
  detected <- sum(x >= 1, na.rm = TRUE)
  frac <- detected / n
  
  if (detected == 0) return("Not detected")
  if (detected == n) return("Detected in all")
  if (frac >= 0.66) return("Detected in many")
  if (detected == 1) return("Detected in single")
}

# Do the classification
expr <- as_tibble(tissue_avg, rownames = "gene_id")

tissue_cols <- setdiff(colnames(expr), "gene_id")


expr_classified <- expr %>%
  rowwise() %>%
  mutate(
    specificity_class = classify_gene(c_across(all_of(tissue_cols))),
    detection_class   = detect_class(c_across(all_of(tissue_cols)))
  ) %>%
  ungroup()

# ------------------------------
# Save results
# ------------------------------
write.csv(expr, file = file.path(wkdir, "csv/NX_tissue_matrix_average.csv"))
write.csv(expr_classified, file = file.path(wkdir, "csv/NX_tissue_matrix_classified.csv"))


single <- expr_classified %>% filter(detection_class == "Detected in single" , specificity_class == "Tissue enriched")

single  <- single %>%
  rowwise() %>%
  mutate(tissue_specfic = tissue_cols[which.max(c_across(all_of(tissue_cols)))]) %>%
  ungroup()

counts <- single %>% dplyr::count(tissue_specfic)

counts

