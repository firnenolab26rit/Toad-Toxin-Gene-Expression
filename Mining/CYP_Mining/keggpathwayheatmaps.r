library(tidyverse)
library(ggplot2)
library(ggpattern)

# ================================
# 1. Get genes from a KEGG term
# ================================
get_genes_from_term <- function(gsea, eggnog, term) {
  
  df <- gsea@result %>%
    filter(Description == term)
  
  pattern <- df %>%
    pull(core_enrichment) %>%
    str_split("/") %>%
    unlist()
  
  eggnog %>%
    filter(Reduce(`|`, lapply(pattern, function(p) grepl(p, trinity_id))))
}


# ================================
# 2. Prepare TPM summary
# ================================
prepare_tpm_summary <- function(TPM_path, design_path, genes_df) {
  
  TPM <- read_csv(TPM_path)
  colnames(TPM)[1] <- "trinity_id"
  
  filtered_TPM <- TPM %>%
    filter(trinity_id %in% genes_df$trinity_id) %>%
    filter(rowSums(select(., -trinity_id) >= 1) > 0)
  
  TPM_long <- filtered_TPM %>%
    pivot_longer(-trinity_id, names_to = "sample_name", values_to = "TPM")
  
  design <- read_csv(design_path)
  
  TPM_meta <- TPM_long %>%
    left_join(design, by = "sample_name") %>%
    select(trinity_id, TPM, tissue, timepoint)
  
  TPM_meta %>%
    group_by(trinity_id, tissue, timepoint) %>%
    summarise(mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
    mutate(log2_mean_TPM = log2(mean_TPM + 1))
}


# ================================
# 3. Add DEG flags per tissue
# ================================
add_deg_flags <- function(TPM_summary, eggnog, tissues, deg_base_path,
                          padj_cutoff = 0.05, lfc_cutoff = 1) {
  
  map_dfr(tissues, function(i) {
    
    deg_files <- list.files(
      path = file.path(deg_base_path, i),
      pattern = "DEG_.*\\.csv",
      full.names = TRUE
    )
    
    deg_df <- deg_files %>%
      set_names() %>%
      map_dfr(~ read_csv(.x), .id = "file") %>%
      mutate(
        timepoint = str_extract(file, "[0-9]+h"),
        tissue = i
      )
    
    colnames(deg_df)[2] <- "trinity_id"
    
    deg_flag <- deg_df %>%
      mutate(is_DEG = padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff) %>%
      filter(is_DEG) %>%
      distinct(trinity_id, timepoint, tissue) %>%
      mutate(DEG_flag = TRUE)
    
    TPM_summary %>%
      filter(tissue == i) %>%
      left_join(deg_flag, by = c("trinity_id", "timepoint", "tissue")) %>%
      mutate(DEG_flag = coalesce(DEG_flag, FALSE)) %>%
      left_join(
        eggnog %>% distinct(trinity_id, .keep_all = TRUE) %>%
          select(trinity_id, Preferred_name),
        by = "trinity_id"
      )
  })
}


# ================================
# 4. Plot function
# ================================
plot_heatmap <- function(df, time_levels, tissue_levels) {
  
  df <- df %>%
    mutate(
      tissue = factor(tissue, levels = tissue_levels),
      timepoint = factor(timepoint, levels = time_levels),
      label = if_else(
        !is.na(Preferred_name) & trimws(Preferred_name) != "",
        paste0(Preferred_name, " (", trinity_id, ")"),
        trinity_id
      )
    )
  
  ggplot(df, aes(x = timepoint, y = label, fill = log2_mean_TPM)) +
    geom_tile() +
    
    geom_tile_pattern(
      data = subset(df, DEG_flag),
      pattern = "stripe",
      pattern_fill = "black",
      pattern_angle = 45,
      pattern_density = 0.1,
      pattern_spacing = 0.02,
      fill = NA,
      color = "black",
      linewidth = 0.6
    ) +
    
    scale_fill_gradient(low = "lightgrey", high = "red") +
    facet_wrap(~tissue) +
    
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_text(size = 16, face = "bold", color = "black"),
      axis.text  = element_text(size = 14, face = "bold", color = "black"),
      legend.text = element_text(face = "bold", size = 14, color = "black"),
      legend.title = element_text(face = "bold", size = 16, color = "black"),
      axis.ticks.y = element_line(color = "black"), 
      strip.text.x = element_text(size = 16, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      legend.position = "top"
    ) +
    labs(
      x = "Timepoint",
      fill = expression(log[2](bar(TPM) + 1))
    )
}

# The pipeline 
# Need to load in gsea, eggnog, and terms
eggnog <- read_tsv("../../Annotation/Gene_Level_Annotations_Wide.tsv")
gsea <- read_rds("../../GSEA/results/Tissue/Deseq_Modules_Input_GSEA.rds")
terms <- c(
  "Glycine, serine and threonine metabolism",
  "Estrogen signaling pathway",
  "Metabolism of xenobiotics by cytochrome P450",
  "Retinol metabolism",
  "Cholesterol metabolism",
  "Steroid hormone biosynthesis",
  "Drug metabolism",
  "Microbial metabolism in diverse environments",
  "Ascorbate and aldarate metabolism",
  "Drug metabolism",
  "Biosynthesis of secondary metabolites",
  "Protein export",
  "Histidine metabolism",
  "Fat digestion and absorption",
  "Arginine biosynthesis",
  "Carbon metabolism",
  "Glycosphingolipid biosynthesis",
  "Biosynthesis of amino acids",
  "beta-Alanine metabolism",
  "One carbon pool by folate",
  "Glyoxylate and dicarboxylate metabolism"
)

plots <- map(terms, function(term) {
  
  genes <- get_genes_from_term(gsea, eggnog, term)
  
  TPM_summary <- prepare_tpm_summary(
    "../../Tissue_Expression/reduc_gene_tpm_matrix.csv",
    "../../WGCNA/Design_Formated.csv",
    genes
  )
  
  combined <- add_deg_flags(
    TPM_summary,
    eggnog,
    tissues = c("parotoid", "liver"),
    deg_base_path = "~/Downloads/TOAD_ANALYSIS/MultiTissue_RNAseq/DEseq/WithinTissue/DE_results/"
  )
  
  plot_heatmap(
    combined,
    time_levels = c("con","0h","3h","6h","9h","12h"),
    tissue_levels = c("muscle","parotoid","liver")
  ) +
    ggtitle(term)
})



# Save them to a pdf
pdf("KeggPathwayGenesByTissue_heatmaps.pdf", width = 12, height = 24)

walk(plots, print)

dev.off()

