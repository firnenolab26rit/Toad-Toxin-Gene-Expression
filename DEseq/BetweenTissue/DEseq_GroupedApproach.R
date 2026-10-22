################################################################################
#     Time Series Differential Expression Analysis: Toad Multi-Tissue
#     Author: Michael Grapin (Refactored/Expanded)
################################################################################
library(DESeq2)
library(apeglm)
library(dplyr)
library(tidyverse)
library(stringr)


# ==============================================================================
# Phase 1: Setup & Preprocessing
# ==============================================================================

  #-----------------------------
  # 1. Read in files
  #-----------------------------
  counts_file  <- "../fixed_reduc_gene_counts_matrix.csv"
  design_file  <- "../design.csv"
  
  count_data <- read_csv(counts_file) |> as.data.frame()
  design     <- read_csv(design_file) |> as.data.frame()

  #-----------------------------
  # 2. Format count matrix
  #-----------------------------
  
  # Set gene IDs as rownames and remove Geneid column
  rownames(count_data) <- count_data[[1]]
  count_data <- count_data[, -1]
  
  # Convert counts to integers (required for DESeq2)
  count_data[] <- lapply(count_data, as.integer)

  #-----------------------------
  # 3. Format design matrix
  #-----------------------------
  
  # Set sample IDs as rownames
  rownames(design) <- design$sample_name
  
  # Keep only relevant columns
  design <- design |> 
    select(tissue, timepoint)

  # Create combined condition column
  design <- design |>
    mutate(
      group = factor(paste(tissue, timepoint, sep = "_"), levels = c(
        "muscle_con",
        "liver_con",
        "liver_0h",
        "liver_3h",
        "liver_6h",
        "liver_9h",
        "liver_12h",
        "parotoid_con",
        "parotoid_0h",
        "parotoid_3h",
        "parotoid_6h",
        "parotoid_9h",
        "parotoid_12h"
      )),
      tissue    = factor(tissue, levels = c("muscle", "liver", "parotoid")),
      timepoint = factor(timepoint, levels = c("con", "0h", "3h", "6h", "9h", "12h")), 
      individual = str_extract(rownames(design), "[A-Z]{3}[0-9]{3}")
    )
  
  #-----------------------------
  # Set a cutoff for noise 
  #-----------------------------
  smallestGroupSize <- 3  # triplicates
  
  count_data <- count_data %>%
    filter(
      rowSums(across(everything(), ~ . >= 10)) >= smallestGroupSize
    )
  
  #-----------------------------
  # 4. Ensure matching sample order
  #-----------------------------
  
  count_data <- count_data[, rownames(design)]
  
  # Safety check
  stopifnot(all(colnames(count_data) == rownames(design)))


# Assuming 'counts_matrix' and 'colData' are already loaded in your environment.

# 1. Create a composite group variable for easy pairwise Wald testing
# This combines Tissue and Timepoint (e.g., "Liver_Control", "Liver_3hr")
#colData$group <- factor(paste0(colData$tissue, "_", colData$timepoint))

# Note on Muscle: Since Muscle only has control timepoints, including it in 
# full time-series LRT models can cause rank-deficiency errors. 
# Best practice: Subset the data for the global temporal LRTs.
colData_time <- design %>% filter(tissue %in% c("liver", "parotoid")) %>% droplevels()
counts_time <- count_data[, rownames(colData_time)]

# Safety check
stopifnot(all(colnames(counts_time) == rownames(colData_time)))

# ==============================================================================
# Phase 2: Global Dynamics via Likelihood Ratio Test (LRT)
# ==============================================================================

#' Run LRT to find global dynamic signatures
#' @param counts Matrix of raw counts
#' @param meta Metadata dataframe
#' @param full_model Formula for the full model
#' @param reduced_model Formula for the reduced model
#' @return DESeqDataSet object run with LRT
run_global_lrt <- function(counts, meta, full_model, reduced_model) {
  dds_lrt <- DESeqDataSetFromMatrix(countData = counts, 
                                    colData = meta, 
                                    design = full_model)
  dds_lrt <- DESeq(dds_lrt, test="LRT", reduced=reduced_model)
  return(dds_lrt)
}

# BIOLOGICAL REASONING: The Interaction (Tissue-Specific Responders)
# Question: Which genes respond differently to the stimulus over time depending on the tissue?
# (Q1: How does gene expression vary between tissue + timepoint?)

# Biology: Identifies genes linked to specialized tissue function (e.g., rapid toxin 
# synthesis in the parotoid gland vs. systemic metabolic buffering in the liver).
dds_interaction <- run_global_lrt(
  counts = counts_time, 
  meta = colData_time,
  full_model = ~  tissue + timepoint + tissue:timepoint,
  reduced_model = ~  tissue + timepoint
)

# Already did in th WithinTissue.R Script
# Probably want to do this with just parotoid and just liver because right now we can't separate liver vs parotoid DEG's.
# Want to look if its better to do two models or keep one and just use the tissue altas to separate out results. 

# Extracts a lists of list
res_interaction <- results(dds_interaction) # Genes significant here have different temporal slopes
# just makes a formal df
res_interaction <- as.data.frame(results(dds_interaction)) # Genes significant here have different temporal slopes


# BIOLOGICAL REASONING: Any Time Effect (Systemic Progression)
# Question: Which genes change over time consistently across both Liver and Parotoid?
# (Q2: How does gene expression vary across timepoint?)


## Don't think we can do this because indv and tissue are dependent
# Biology: Captures core systemic responses, such as global stress pathways, 
# circadian rhythms, or immune mobilization triggered by the experimental baseline.
dds_temporal <- run_global_lrt(
  counts = counts_time,
  meta = colData_time,
  full_model = ~ tissue + timepoint,
  reduced_model = ~ tissue
)

res_temporal <- as.data.frame(results(dds_temporal)) 

# Think through this 
# and this 
# full_model = ~ tissue + timepoint,
# reduced_model = ~ timepoint
###########################################################################################
## This is the LTR model that is used in the GSEA ranking parameter.
# How does tissue explain DEGs patterns overall
dds_tissue <- run_global_lrt(
  counts = counts_time, 
  meta = colData_time,
  full_model = ~ tissue ,
  reduced_model = ~ 1 
)

source("../WithinTissue/WithinTissue_Functions.R")

# 1. Run Shrinkage (only once)
res <- lfcShrink(
  dds_tissue,
  coef = "tissue_parotoid_vs_liver",
  type = "ashr"
)

# 2. Generate the Volcano Plot using your ggplot function
# We add the title here so it's ready for the PDF
p <- plot_volcano(res, lfc_threshold = 1, padj_threshold = 0.05) + 
  labs(title = "tissue_parotoid_vs_liver")

print(p)
ggsave(filename = "DE_results/BetweenTissue_LTR/BetweenTissue_ParotoidLiver_LTR_VolcanoPlot.pdf", plot = p, dpi = 300, width = 20, height = 8.27, units = "in")

res <- as.data.frame(res) %>% rownames_to_column("Gene")
write_csv(res, file = "DE_results/BetweenTissue_LTR/BetweenTissue_ParotoidLiver_LTR.csv", )

# 3. Filter for significant results
res_sig <- subset(res, abs(log2FoldChange) > 1 & padj < 0.05)
res_sig <- as.data.frame(res_sig) %>% rownames_to_column("Gene")
write_csv(res_sig, file = "DE_results/BetweenTissue_LTR/BetweenTissue_ParotoidLiver_LTR_SigDEGS.csv")

# ==============================================================================
# Phase 3: Specific Pairwise Dynamics (Wald Tests using Group Variable)
# ==============================================================================
# Here, we use the FULL dataset (including Muscle) because the 'group' variable 
# sidesteps the unbalanced design issue.

dds_wald <- DESeqDataSetFromMatrix(countData = count_data, 
                                   colData = design, 
                                   design = ~ group)
dds_wald <- DESeq(dds_wald, test="Wald")

# ==============================================================================
# Phase 3: Downstream Filtering and Analysis
# ==============================================================================
source("../WithinTissue/WithinTissue_Functions.R")
# vectors
comp <- c("liver_con",
          "liver_0h",
          "liver_3h",
          "liver_6h",
          "liver_9h",
          "liver_12h",
          "parotoid_con",
          "parotoid_0h",
          "parotoid_3h",
          "parotoid_6h",
          "parotoid_9h",
          "parotoid_12h")

ref  <- "muscle_con"

# selected_tissues <- names(dds_list)[c(1,2)]

# Containers
res_list <- list()

# for (e in selected_tissues) {
# 
#   message("\n--- Processing tissue: ", e, " ---")
#   dds_obj <- dds_list[[e]]$wald
dds_obj <- dds_wald

  
  # 1. Map through comparisons
  # returns a list of lists (each containing $df and $plot)
  tissue_output <- map(comp, ~extract_results(dds_obj, "group", .x, ref))
  names(tissue_output) <- paste0(comp, "_vs_", ref)
  
  # 2. Extract components for this specific tissue
  tissue_dfs <- map(tissue_output, "df")
  tissue_plots <- map(tissue_output, "plot")
  
  # 3. Create directory
  out_dir <- file.path("DE_results", "BetweenTissue")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 4. Save CSVs for this tissue
  iwalk(tissue_dfs, function(df, comp_name) {
    write.csv(
      df, 
      file = file.path(out_dir, paste0("DEG_", "BetweenTissue", "_", comp_name, ".csv")),
      row.names = TRUE
    )
  })
  
  # 5. START PDF DEVICE (Specific to this tissue)
  pdf_path <- file.path(out_dir, paste0("Volcano_Plots_", "BetweenTissue", "WaldComparison", ".pdf"))
  message("Saving plots to: ", pdf_path)
  
  pdf(pdf_path, width = 11.69, height = 8.27)
  walk(tissue_plots, print)
  dev.off() # Close device for this tissue
  
  # Store the dataframes in the master list if needed for later R analysis
  res_list[["BetweenTissue_Wald"]] <- tissue_dfs
  
#   message(e, " processing complete.")
# }

# Check to make sure the names are correct (2 tissues x 5 timepoints == 10)
sum(map_int(res_list, length))



BetweenTissueWald_degs <- get_deg_counts(res_list$BetweenTissue_Wald)
colnames(BetweenTissueWald_degs)[1] <- "group"

BetweenTissueWald_degs$group <- factor(
  BetweenTissueWald_degs$group,
  levels = c(
    "liver_con",
    "liver_0h",
    "liver_3h",
    "liver_6h",
    "liver_9h",
    "liver_12h",
    "parotoid_con",
    "parotoid_0h",
    "parotoid_3h",
    "parotoid_6h",
    "parotoid_9h",
    "parotoid_12h"
))

# Mutate to add a fill athestic for tissue 
BetweenTissueWald_degs <- BetweenTissueWald_degs %>%
  mutate(
    tissue = stringr::word(group, 1, sep = "_"),
    tissue = factor(tissue, levels = c("parotoid", "liver")),
    timepoint = stringr::word(group, 2, sep = "_"),
    timepoint = factor(timepoint, levels = c("con","0h","3h","6h","9h","12h"))
  )

BetweenTissueWald_degs <- BetweenTissueWald_degs %>% filter(timepoint != "con")
write_csv(BetweenTissueWald_degs, file = 'BetweenTissuePairwiseWaldDEGsUpDown.csv')

## The plotting 
# The plot (Barplot code)
plot2 <- ggplot(
  BetweenTissueWald_degs,
  aes(
    x = timepoint,
    y = deg_counts_plot,
    fill = tissue,
    pattern = tissue
  )
) +

  geom_col_pattern(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.4,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.2,
    pattern_spacing = 0.05
  ) +

  geom_hline(
    yintercept = 0,
    linetype = "solid",
    color = "black",
    linewidth = 1.0
  ) +

  labs(
    x = "Between Tissue Comparisons",
    y = "Number of DEGs",
    title = "",
    fill = "Tissue",
    pattern = "Tissue"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1.0),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 14, face = "bold"),
    legend.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 16)
  ) +

  scale_fill_manual(
    values = c(
      "parotoid" = "#E65700",  # choose your color
      "liver"    = "#03548a"
    )
  ) +

  scale_pattern_manual(
    values = c(
      "parotoid" = "none",
      "liver"    = "none"
    )
  )

print(plot2)

saveRDS(plot2, file = "BetweenTissue_BarPlot.rds")

#' # I can write a function to iterate over all the combinations and save them as a list
#' 
#' #' Extract pairwise contrast safely
#' #' @param dds DESeq object run with Wald test
#' #' @param group1 The numerator of the contrast (e.g., "Liver_3hr")
#' #' @param group2 The denominator/baseline of the contrast (e.g., "Liver_Control")
#' get_pairwise_contrast <- function(dds, group1, group2) {
#'   # BIOLOGICAL REASONING: 
#'   # This tests "What is happening precisely at Time X vs Baseline within Tissue Y?"
#'   # It is highly sensitive but might capture noisy, timepoint-specific spikes 
#'   # that lack a broader trend.
#'   res <- results(dds, contrast=c("group", group1, group2))
#'   return(res)
#' }
#' 
#' # Examples of within-tissue temporal comparisons
#' res_liver_3hr_vs_ctrl <- get_pairwise_contrast(dds_wald, "Liver_3hr", "Liver_Control")
#' res_parotid_12hr_vs_ctrl <- get_pairwise_contrast(dds_wald, "Parotid_12hr", "Parotid_Control")
#' 
#' # Examples of baseline tissue differences (Using your null anchor)
#' # BIOLOGICAL REASONING: 
#' # Muscle control serves as the absolute somatic baseline (non-secretory, non-metabolic hub).
#' # Comparing Liver/Parotoid controls to Muscle control isolates genes structurally 
#' # required for liver/parotoid identity prior to any experimental stimulus.
#' res_liver_vs_muscle_ctrl <- get_pairwise_contrast(dds_wald, "Liver_Control", "Muscle_Control")
#' res_parotid_vs_muscle_ctrl <- get_pairwise_contrast(dds_wald, "Parotid_Control", "Muscle_Control")
#' 
#' # ==============================================================================
#' # Phase 4: LFC Shrinkage Automation
#' # ==============================================================================
#' # Shrinkage (apeglm) is essential for ranking genes for GSEA or interpreting 
#' # biological magnitude, as it penalizes high fold-changes driven by low read counts.
#' 
#' #' Apply apeglm shrinkage to all non-intercept coefficients
#' #' @param dds A run DESeq object
#' #' @return A named list of shrunk results objects
#' shrink_all_coefs <- function(dds) {
#'   coef_names <- resultsNames(dds)[-1] # Remove Intercept
#'   
#'   shrunk_list <- lapply(coef_names, function(coef_name) {
#'     message("Applying apeglm shrinkage to: ", coef_name)
#'     lfcShrink(dds, coef = coef_name, type = "apeglm")
#'   })
#'   
#'   names(shrunk_list) <- coef_names
#'   return(shrunk_list)
#' }
#' 
#' # Run the automation on our Wald dds
#' shrunken_wald_results <- shrink_all_coefs(dds_wald)
#' 
#' # ==============================================================================
#' # Phase 5: Filtering and Intersecting with Tissue Atlas
#' # ==============================================================================
#' 
#' #' Extract significant DE genes based on padj and LFC thresholds
#' #' @param res_obj DESeq results object
#' #' @param p_val Adjusted p-value threshold (default 0.05)
#' #' @param lfc Log2 Fold Change threshold (default 1)
#' #' @return Character vector of significant Ensembl/Gene IDs
#' extract_sig_genes <- function(res_obj, p_val = 0.05, lfc = 1) {
#'   res_df <- as.data.frame(res_obj)
#'   sig_genes <- rownames(res_df %>% 
#'                           filter(padj < p_val & abs(log2FoldChange) >= lfc))
#'   return(sig_genes)
#' }
#' 
#' # Example usage:
#' liver_early_genes <- extract_sig_genes(res_liver_3hr_vs_ctrl)
#' liver_baseline_genes <- extract_sig_genes(res_liver_vs_muscle_ctrl)

# BIOLOGICAL REASONING: 
# If a gene is identified in 'liver_early_genes' but NOT in your pre-existing Tissue Atlas,
# it means the gene is strictly "stress/stimulus-induced" rather than a constitutive 
# marker of liver identity.

# ------------------------------------------------------------------------------
