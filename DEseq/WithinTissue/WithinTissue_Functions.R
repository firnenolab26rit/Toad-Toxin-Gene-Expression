# within tissue functions

# Run the deseq pipeline
run_deseq_by_tissue <- function(tissue_name, design_df, count_mat) {
  
  # Subset metadata
  colData_sub <- design_df %>%
    filter(tissue == tissue_name) %>%
    droplevels()
  
  # Subset counts
  counts_sub <- count_mat[, rownames(colData_sub)]
  
  # Safety check
  stopifnot(all(colnames(counts_sub) == rownames(colData_sub)))
  
  # Set reference level
  colData_sub$timepoint <- relevel(factor(colData_sub$timepoint), ref = "con")
  
  # Build DESeq object
  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData   = colData_sub,
    design    = ~ timepoint
  )
  
  # Run Wald
  dds_wald <- DESeq(dds, test = "Wald")
  
  # Run LRT (need to re-run from raw dds, not from dds_wald)
  dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
  
  return(list(
    wald = dds_wald,
    lrt  = dds_lrt
  ))
}

##============================================================================##

extract_results <- function(deseq_object, design_var, comparison, reference) {
  
  coef_name <- paste0(design_var, "_", comparison, "_vs_", reference)
  message("Processing: ", coef_name)
  
  # 1. Run Shrinkage (only once)
  res <- lfcShrink(
    deseq_object,
    coef = coef_name,
    type = "ashr"
  )
  
  # 2. Generate the Volcano Plot using your ggplot function
  # We add the title here so it's ready for the PDF
  p <- plot_volcano(res, lfc_threshold = 1, padj_threshold = 0.05) + 
    labs(title = coef_name)
  
  # 3. Filter for significant results
  res_sig <- subset(res, abs(log2FoldChange) > 1 & padj < 0.05)
  
  message(comparison, " vs ", reference, " had ", nrow(res_sig), " DE genes")
  
  # 4. Return as a named list
  return(list(
    df = as.data.frame(res_sig),
    plot = p
  ))
}
##============================================================================##

# # Extract pairwise wald test significant results
# extract_results <- function(deseq_object, design_var, comparison, reference) {
#   
#   coef_name <- paste0(design_var, "_", comparison, "_vs_", reference)
#   
#   # lfc shrinkage applied. (method = apeglm)
#     res <- lfcShrink(
#     deseq_object,
#     # if using model coefficents
#     coef = coef_name,
#     
#     # if want to do custom comparisons
#     #contrast = c(design_var, comparison, reference),
#     
#     # Use with coef
#     type = "apeglm"
#     # Use with contrast
#     #type = "ashr"
#   )
#   
#   # Extract results to df without any shrinkage
#   # res <- results(
#   #   deseq_object,
#     
#     # chose coef or contrast method based on comparison
#   #   #contrast = c(design_var, comparison, reference),
#   #   coef = coef_name,
#   #   alpha = 0.05
#   # )
#   
#   message("Extracting ", comparison, " vs ", reference, "...")
#   
#   res_sig <- subset(res, abs(log2FoldChange) > 1 & padj < 0.05)
#   
#   message(comparison, " vs ", reference, 
#           " had ", nrow(res_sig), " DE genes")
#   
#   return(as.data.frame(res_sig))
# }

# extract inputs to accept if wanted. (will leave it out for now)
# c = c("coef" == "apeglm", "contrast" == "ashr)

##============================================================================##
# Get Plotting DF for DEG's on comparisons
get_deg_counts <- function(res_list, lfc_threshold = 1) {
  
  map_dfr(names(res_list), function(nm) {
    
    df <- res_list[[nm]]
    tp <- sub("_vs_.*", "", nm)
    
    tibble(
      timepoint = tp,
      reg_dir = c("up", "down"),
      deg_counts = c(
        sum(df$log2FoldChange >  lfc_threshold, na.rm = TRUE),
        sum(df$log2FoldChange < -lfc_threshold, na.rm = TRUE)
      )
    )
    
  }) %>%
    
    # Add negative plotting column
    mutate(
      deg_counts_plot = deg_counts * ifelse(reg_dir == "down", -1, 1)
    )
}

##============================================================================##

# Extract LTR model DEG's 
extract_LRT_DEGs <- function(LRT_df,
                             tissue,
                             model = "LRT",
                             lfc_threshold = 1,
                             alpha = 0.05) {
  
  # Filter significant genes
  sig_df <- LRT_df %>%
    filter(abs(log2FoldChange) > 1 & padj < 0.05)
  
  # Count up/down
  up_count <- sum(sig_df$log2FoldChange >  lfc_threshold, na.rm = TRUE)
  down_count <- sum(sig_df$log2FoldChange < -lfc_threshold, na.rm = TRUE)
  
  # Build summary dataframe
  out_df <- tibble(
    tissue = tissue,
    model = model,
    reg_dir = c("up", "down"),
    deg_counts = c(up_count, down_count)
  ) %>%
    mutate(
      deg_counts_plot = deg_counts * ifelse(reg_dir == "down", -1, 1)
    )
  
  return(out_df)
}


##============================================================================##

# 1. VOLCANO PLOT FUNCTION (With specific linewidths and classic theme)
plot_volcano <- function(resultsObject, lfc_threshold = 1, padj_threshold = 0.05) {
  df <- as.data.frame(resultsObject)
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  
  df$Significance <- "Not Significant"
  sig_idx <- df$padj < padj_threshold & abs(df$log2FoldChange) > lfc_threshold
  df$Significance[sig_idx] <- "Significant"
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("Not Significant" = "grey70", "Significant" = "firebrick")) +
    
    # Vertical line at 0 (linewidth 0.5, dotted)
    geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
    
    # Threshold lines (linewidth 1, dotdash)
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), 
               linetype = "dotdash", color = "black", linewidth = 1) +
    geom_hline(yintercept = -log10(padj_threshold), 
               linetype = "dotdash", color = "black", linewidth = 1) +
    
    # Updated Y-axis label
    labs(
      title = " ",
      x = expression(Log[2]~fold~change),
      y = expression(-log[10]~P~Adjusted)
    ) +
    
    # The requested classic theme and sizing
    theme_classic() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 14),
      legend.position = "none"
    )
}

# 2. PCA FROM VOLCANO HITS FUNCTION
plot_pca_from_sig <- function(vsd, resultsObject, lfc_threshold = 1, padj_threshold = 0.05, intgroup = "timepoint") {
  res_df <- as.data.frame(resultsObject)
  sig_genes <- rownames(subset(res_df, padj < padj_threshold & abs(log2FoldChange) > lfc_threshold))
  
  if(length(sig_genes) < 3) {
    stop("Too few significant genes found to perform PCA.")
  }
  
  mat <- t(assay(vsd)[sig_genes, ])
  pca_res <- prcomp(mat, scale. = TRUE)
  percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)))
  
  pca_df <- as.data.frame(pca_res$x)
  pca_df$Sample <- rownames(pca_df)
  metadata <- as.data.frame(colData(vsd))
  pca_df <- merge(pca_df, metadata, by.x = "Sample", by.y = "row.names")
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[intgroup]])) +
    geom_point(size = 4) +
    stat_ellipse(aes(fill = .data[[intgroup]]), alpha = 0.1, geom = "polygon") +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(title = "PCA of Significant DE Genes",
         subtitle = paste(length(sig_genes), "genes used for clustering")) +
    theme_classic() + # Consistency with Volcano plot theme
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 14)
    )
}

# --- EXAMPLE EXECUTION ---
# volcano_p <- plot_volcano(res)
# pca_p <- plot_pca_from_sig(vsd, res, intgroup = "Treatment")
# library(patchwork); volcano_p + pca_p


##============================================================================##

