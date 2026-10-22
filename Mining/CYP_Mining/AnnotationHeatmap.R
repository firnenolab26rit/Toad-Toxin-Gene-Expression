library(tidyverse)
library(ggplot2)
library(ggpattern)

# Could also look at broader pathway gens in these two tissues
# ================================
# 1. Load annotation + filter genes
# ================================
eggnog <- read_tsv("../../Annotation/Gene_Level_Annotations_Wide.tsv")

# Have to look up other names for these genes
# Give a pattern of gene  names
pattern <- c(
  "CYP11A1",
  "CYP46A1",
  "CYP39A1",
  "CYP7B1",
  "CYP8B1",
  "CH25H",
  "CYP27A1",
  "AMACR",
  "HSD17B4",
  "SCP2",
  "BAAT",
  "ACOX2"
)

pattern <- c(
  "CYP11A1",
  "CYP46A1",
  "CYP39A1",
  "CYP7B1",
  "CYP8B1",
  "CH25H",
  "CYP27A1",
  "AMACR",
  "HSD17B4",
  "SCP2",
  "BAAT",
  "ACOX2", 
  "ELOVL7", 
  "PTPLAD2",
  "PTPLA", 
  "ABCD3", 
  "AKR1C4", 
  "ACOT"
)


pattern <- c("CYP7A1", "HSD3B7", "CYP8B1", "AKR1D1", "AKR1C4", "CYP27A1")

pattern <- "[Pp]450"
mine <- eggnog %>%
  filter(Reduce(`|`, lapply(pattern, function(p) grepl(p, PFAMs))))
#mine <- eggnog %>%
#filter(if_any(everything(), ~ str_detect(.x, pattern)))


mine <- eggnog %>%
  filter(Reduce(`|`, lapply(pattern, function(p) grepl(p, Preferred_name))))

mine <- eggnog %>%
  filter(Reduce(`|`, lapply(pattern, function(p) grepl(p, trinity_id))))

# ================================
# 2. Load TPM matrix
# ================================
TPM <- read_csv("../../Tissue_Expression/reduc_gene_tpm_matrix.csv")
colnames(TPM)[1] <- "trinity_id"

# Subset genes of interest
filtered_TPM <- TPM %>%
  filter(trinity_id %in% mine$trinity_id)

# Keep genes with at least one TPM >= 1
filtered_TPM <- filtered_TPM %>%
  filter(rowSums(select(., -trinity_id) >= 1) > 0)

# ================================
# 3. Pivot + join metadata
# ================================
TPM_long <- filtered_TPM %>%
  pivot_longer(
    cols = -trinity_id,
    names_to = "sample_name",
    values_to = "TPM"
  )

design <- read_csv("../../WGCNA/Design_Formated.csv")

TPM_meta <- TPM_long %>%
  left_join(design, by = "sample_name") %>%
  select(trinity_id, sample_name, TPM, tissue, timepoint, group) %>%
  dplyr::rename(tissue_timepoint = group)

# ================================
# 4. Compute mean TPM + log2
# ================================
TPM_summary <- TPM_meta %>%
  group_by(trinity_id, tissue, timepoint) %>%
  summarise(mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(log2_mean_TPM = log2(mean_TPM + 1)
  )


gene_order <- TPM_summary %>%
  group_by(trinity_id) %>%
  summarise(var_all = var(log2_mean_TPM)) %>%
  arrange(desc(var_all)) %>%
  pull(trinity_id)

# Factor levels
timepoint_levels <- c("con", "0h", "3h", "6h", "9h", "12h")
tissue_levels <- c("muscle", "parotoid", "liver")

TPM_summary <- TPM_summary %>%
  mutate(
    tissue = factor(tissue, levels = tissue_levels),
    timepoint = factor(timepoint, levels = timepoint_levels)
  )

##################################################################

# ================================
# 1. USER INPUT
# ================================

tissues <- c("parotoid", "liver")

combined <- map_dfr(tissues, function(i) {
  
  # -------------------------------
  # 1. Read DEG files
  # -------------------------------
  deg_files <- list.files(
    path = paste0("~/Downloads/TOAD_ANALYSIS/MultiTissue_RNAseq/DEseq/WithinTissue/DE_results/", i),
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
  
  # -------------------------------
  # 2. DEG flag
  # -------------------------------
  deg_flag <- deg_df %>%
    mutate(is_DEG = padj < 0.05 & abs(log2FoldChange) > 1) %>%
    filter(is_DEG) %>%
    distinct(trinity_id, timepoint, tissue) %>%
    mutate(DEG_flag = TRUE)
  
  # -------------------------------
  # 3. TPM + join
  # -------------------------------
  plot_df_specific <- TPM_summary %>%
    filter(tissue == i) %>%
    left_join(deg_flag, by = c("trinity_id", "timepoint", "tissue")) %>%
    mutate(DEG_flag = coalesce(DEG_flag, FALSE)) %>%
    left_join(
      eggnog %>% distinct(trinity_id, .keep_all = TRUE) %>% select(trinity_id, Preferred_name),
      by = "trinity_id"
    )
  
  # -------------------------------
  # 4. Order genes (within tissue)
  # -------------------------------
  plot_df_specific <- plot_df_specific %>%
    group_by(trinity_id) %>%
    mutate(var_gene = var(log2_mean_TPM, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(var_gene)) %>%
    mutate(trinity_id = factor(trinity_id, levels = unique(trinity_id)))
  
  return(plot_df_specific)
})

# Factor levels
timepoint_levels <- c("con", "0h", "3h", "6h", "9h", "12h")
tissue_levels <- c("muscle", "parotoid", "liver")

combined <- combined %>%
  mutate(
    tissue = factor(tissue, levels = tissue_levels),
    timepoint = factor(timepoint, levels = timepoint_levels)
  )


# Alternative to make a unique preferred anme 
combined <- combined  %>%
mutate(label = paste0(Preferred_name, " (", trinity_id, ")"))
## ggplot is just taking the first gene name or last and plotting it
combined <- combined %>%
  mutate(label = if_else(
    !is.na(Preferred_name) & trimws(Preferred_name) != "",
    paste0(Preferred_name, " (", trinity_id, ")"),
    trinity_id
  ))


# ================================
# 7. Heatmap with DEG overlay
# ================================
# Lets make this plot pretty 

p <- ggplot(combined, aes(x = timepoint, y = label, fill = log2_mean_TPM)) +
  
  # Base heatmap
  geom_tile() +
  
  # Bold rectangle around differentially expressed cells
  geom_tile(
    data = subset(combined, !is.na(DEG_flag) & DEG_flag),
    aes(x = timepoint, y = label),
    inherit.aes = FALSE,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Fill color scale
  scale_fill_gradient(
    low = "lightgrey",
    high = "red"
  ) +
  
  # Separate heatmaps by tissue
  facet_wrap(~tissue) +
  
  # Color legend
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  
  # Theme
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    
    panel.grid = element_blank(),
    
    # Border around each heatmap panel
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 2
    ),
    
    # Horizontal legend at the top
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    # Axis titles and text
    axis.title = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    
    axis.text = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    
    # Legend text
    legend.text = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 16,
      color = "black"
    ),
    
    # Axis ticks
    axis.ticks.y = element_line(color = "black"),
    
    # Tissue facet labels
    strip.text.x = element_text(
      size = 16,
      face = "bold"
    )
  ) +
  
  # Axis labels and legend title
  labs(
    x = "Timepoint",
    fill = expression(log[2](bar(TPM) + 1))
  )

# Display plot
print(p)

# Number of DEGS
combined %>% filter(DEG_flag == TRUE) %>% nrow()
combined %>% filter(DEG_flag == TRUE, tissue == "liver") %>% nrow()
combined %>% filter(DEG_flag == TRUE, tissue == "parotoid") %>% nrow()


# p <- ggplot(combined, aes(x = timepoint, y = label, fill = log2_mean_TPM)) +
#   geom_tile() +
#   
#   # Cross-hatch overlay for DEGs
#   geom_tile_pattern(
#     data = subset(combined, DEG_flag),
#     pattern = "stripe",
#     pattern_fill = "black",
#     pattern_angle = 45,
#     pattern_density = 0.1,
#     pattern_spacing = 0.02,
#     fill = NA,
#     color = "black",
#     linewidth = 0.6
#   ) +
#   
#   scale_fill_gradient(low = "lightgrey", high = "red") +
#   facet_wrap(~tissue) +
#   
#   guides(
#     fill = guide_colorbar(title.position = "top", title.hjust = 0.5)
#   ) +
#   
#   theme_minimal(base_size = 14) +
#   theme(
#     axis.title.y = element_blank(),
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#     panel.grid = element_blank(),
#     
#     # Border around each heatmap
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
#     
#     # Horizontal legend on top
#     legend.position = "top",
#     legend.direction = "horizontal",
#     legend.box = "horizontal",
#     
#     axis.title = element_text(size = 16, face = "bold", color = "black"),
#     axis.text  = element_text(size = 14, face = "bold", color = "black"),
#     legend.text = element_text(face = "bold", size = 14, color = "black"),
#     legend.title = element_text(face = "bold", size = 16, color = "black"),
#     axis.ticks.y = element_line(color = "black"), 
#     strip.text.x = element_text(size = 16, face = "bold")
#   ) +
#   labs(
#     x = "Timepoint",
#     fill = expression(log[2](bar(TPM) + 1))
#   )
# 
# 
# p


ggsave(filename = "FigS3LongCypHeatap.pdf", plot = p, width = 10, height = 72, units = "in", dpi = 600, limitsize = FALSE)
ggsave(filename = "FigS3LongCypHeatap.jpeg", plot = p, width = 10, height = 72, units = "in", dpi = 600, limitsize = FALSE)


# taking away the deg overlay (and now just plotting the deg == TRUE)
# combined <- combined %>% filter(DEG_flag == TRUE)
# v <- combined %>% filter(DEG_flag == TRUE) %>% pull(trinity_id)
# combined <- combined %>% filter(trinity_id %in% v)
# 
# 
# p <- ggplot(combined, aes(x = timepoint, y = label, fill = log2_mean_TPM)) +
#   geom_tile() +
#   scale_fill_gradient(low = "lightgrey", high = "red") +
#   facet_wrap(~tissue) +
#   
#   guides(
#     fill = guide_colorbar(title.position = "top", title.hjust = 0.5)
#   ) +
#   
#   theme_minimal(base_size = 14) +
#   theme(
#     axis.title.y = element_blank(),
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#     panel.grid = element_blank(),
#     
#     # Border around each heatmap
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
#     
#     # Horizontal legend on top
#     legend.position = "top",
#     legend.direction = "horizontal",
#     legend.box = "horizontal",
#     
#     axis.title = element_text(size = 16, face = "bold", , color = "black"),
#     axis.text  = element_text(size = 14, face = "bold", color = "black"),
#     legend.text = element_text(face = "bold", size = 14, color = "black"),
#     legend.title = element_text(face = "bold", size = 16, color = "black"),
#     axis.ticks.y = element_line(color = "black"), 
#     strip.text.x = element_text(size = 16, face = "bold")
#   ) +
#   labs(
#     x = "Timepoint",
#     fill = expression(log[2](bar(TPM) + 1))
#   )
# 
# 
# p
# 
# 
# 
# write_rds(p,file = "HeatMapGenesOFInterest.rds")
# 9 x 11.69
ggsave(filename = "FigureS3.HeatmapCGLong.pdf", plot = p, width = 10, height = 16, units = "in", dpi = 600)
ggsave(filename = "FigureS3.HeatmapCGLong.pdf", plot = p, width = 10, height = 16, units = "in", dpi = 600)

# I just merged the rows for the main text, and the supplementary long file can be found there.
combined_unique <- combined %>% distinct(Preferred_name, timepoint, tissue, .keep_all = TRUE) %>%
  filter(Preferred_name != "CYP46A1,ITGA2")

p_reduce <- ggplot(combined_unique, aes(x = timepoint, y = Preferred_name, fill = log2_mean_TPM)) +
  
  # Base heatmap
  geom_tile() +
  
  # Bold rectangle around differentially expressed cells
  geom_tile(
    data = subset(combined_unique, !is.na(DEG_flag) & DEG_flag),
    aes(x = timepoint, y = Preferred_name),
    inherit.aes = FALSE,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  
  # Fill color scale
  scale_fill_gradient(
    low = "lightgrey",
    high = "red"
  ) +
  
  # Separate heatmaps by tissue
  facet_wrap(~tissue) +
  
  # Color legend
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  
  # Theme
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    
    panel.grid = element_blank(),
    
    # Border around each heatmap panel
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 2
    ),
    
    # Horizontal legend at the top
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    # Axis titles and text
    axis.title = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    
    axis.text = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    
    # Legend text
    legend.text = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 16,
      color = "black"
    ),
    
    # Axis ticks
    axis.ticks.y = element_line(color = "black"),
    
    # Tissue facet labels
    strip.text.x = element_text(
      size = 16,
      face = "bold"
    )
  ) +
  
  # Axis labels and legend title
  labs(
    x = "Timepoint",
    fill = expression(log[2](bar(TPM) + 1))
  )

# Display plot
print(p_reduce)

ggsave(filename = "Fig3B.Heatmap.pdf", p_reduce, units = "in", width = 10, height = 12, dpi = 600)
write_rds(p_reduce,file = "HeatMapGenesOFInterest.rds")


# Make a multipanel figure 

library(cowplot)
library(readr)

A <- read_rds(file = "../../GSEA/results/Tissue/GSEA_Tissue_KeggPathwaysOfInterest_plot.rds")

B <- read_rds(file = "HeatMapGenesOFInterest.rds")


P <- plot_grid(A, B, labels = c('A', 'B'), label_size = 30)

ggsave2(filename = "Figure3Annotation.pdf", P, units = "in", width = 20, height = 10, dpi = 600)
ggsave2(filename = "Figure3Annotation.jpeg", P, units = "in", width = 20, height = 10, dpi = 600)
ggsave2(filename = "Figure3Annotation.svg", P, units = "in", width = 20, height = 10, dpi = 600)
