pdf(file = "TissueAltas_PCAs.pdf", width = 11.69, height = 8.27)

pca_theme <- theme_classic(base_size = 12, base_family = "sans") +
  theme(
    # Text
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 11,
      color = "grey30",
      margin = margin(b = 10)
    ),
    axis.title = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    
    # Axes
    axis.line = element_line(
      color = "black",
      linewidth = 0.5
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.length = unit(2, "mm"),
    
    # Legend
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.position = "right",
    legend.key = element_blank(),
    legend.background = element_blank(),
    
    # Panel and figure
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )

#----------------------------------
# Helper function to run + plot PCA
#----------------------------------
run_pca_plot <- function(mat, title, scale_flag) {
  
  # remove zero-variance genes
  keep <- apply(mat, 1, sd) > 0
  mat_filt <- mat[keep, , drop = FALSE]
  
  # PCA (DO NOT CHANGE PARAMETERS)
  pca <- prcomp(t(mat_filt), center = TRUE, scale. = scale_flag)
  
  # explained variance
  percent_var <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100, 1)
  
  # dataframe for plotting
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    run = as.character(run),
    individual = design$individual,
    condition = design$condition,
    tissue = design$tissue
  )
  
  # create plot object
  p <- ggplot(pca_df, aes(PC1, PC2, color = tissue, shape = run)) +
    geom_point(size = 3) +
    labs(
      title = title,
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance")
    ) +
    pca_theme
  
  # display plot
  print(p)
  
  # print variance table to console
  cat("\n", title, "\n")
  print(percent_var[1:10])
  
  # return ggplot object
  return(p)
}

####################################
# 1) TPM RAW DATA
####################################
run_pca_plot(
  mat = TPM_data_ordered,
  title = "TPM RAW DATA",
  scale_flag = TRUE
)

FigS1 <- run_pca_plot(
  mat = TPM_data_ordered,
  title = "",
  scale_flag = TRUE
)
ggsave(filename = "FigS1.BioRepPCA.pdf", plot = FigS1, dpi = 600, width = 8, height = 6, units = "in", device = cairo_pdf)
ggsave(filename = "FigS1.BioRepPCA.jpeg", plot = FigS1, dpi = 600, width = 8, height = 6, units = "in")
ggsave(filename = "FigS1.BioRepPCA.svg", plot = FigS1, dpi = 600, width = 8, height = 6, units = "in")


####################################
# 2) TPM -> TMM -> PARETO (before batch correction)
####################################
run_pca_plot(
  mat = TPM_tmm_pareto,
  title = "TPM -> TMM -> Pareto Data",
  scale_flag = FALSE
)

####################################
# 3) AFTER BATCH CORRECTION
####################################
run_pca_plot(
  mat = NX,
  title = "TPM -> TMM -> Pareto -> Remove Batch Effect Data",
  scale_flag = FALSE
)

dev.off()