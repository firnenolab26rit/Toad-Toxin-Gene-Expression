################################################################################
#     Time Series Differential Expression: Within Tissue
#     Author: Michael Grapin 
################################################################################
library(DESeq2)
library(dplyr)
library(tidyverse)
library(stringr)
library(purrr)
library(ggpattern)
library(cowplot)
library(ashr)
library(apeglm)


source("WithinTissue_Functions.R")
# ==============================================================================
# Phase 1: Setup & Reference Levels
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

# Create combined condition column
design <- design |>
  mutate(
    group = factor(paste(tissue, timepoint, sep = "_")),
    tissue    = factor(tissue, levels = c("muscle", "liver", "parotoid")),
    timepoint = factor(timepoint, levels = c("con", "0h", "3h", "6h", "9h", "12h")), 
    individual = str_extract(rownames(design), "[A-Z]{3}[0-9]{3}")
  ) |>
  select(tissue, timepoint)

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

# ==============================================================================
# Phase 2: Model Execution
# ==============================================================================
# Design: Expression ~ tissue-specific time effects
tissues <- c("parotoid", "liver")

dds_list <- map(
  tissues,
  ~ run_deseq_by_tissue(.x, design, count_data)
)

names(dds_list) <- tissues

# Generates
# dds_list$parotoid$wald
# dds_list$parotoid$lrt
# dds_list$liver$wald
# dds_list$liver$lrt

# If I want the list flattened
# dds_flat <- purrr::flatten(dds_list)

# ==============================================================================
# Phase 3: Downstream Filtering and Analysis
# ==============================================================================
# vectors
comp <- c("0h", "3h", "6h", "9h", "12h")
ref  <- "con"

selected_tissues <- names(dds_list)[c(1,2)]

# Containers
res_list <- list()

for (e in selected_tissues) {
  
  message("\n--- Processing tissue: ", e, " ---")
  dds_obj <- dds_list[[e]]$wald
  
  # 1. Map through comparisons
  # returns a list of lists (each containing $df and $plot)
  tissue_output <- map(comp, ~extract_results(dds_obj, "timepoint", .x, ref))
  names(tissue_output) <- paste0(comp, "_vs_", ref)
  
  # 2. Extract components for this specific tissue
  tissue_dfs <- map(tissue_output, "df")
  tissue_plots <- map(tissue_output, "plot")
  
  # 3. Create directory
  out_dir <- file.path("DE_results", e)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 4. Save CSVs for this tissue
  iwalk(tissue_dfs, function(df, comp_name) {
    write.csv(
      df, 
      file = file.path(out_dir, paste0("DEG_", e, "_", comp_name, ".csv")),
      row.names = TRUE
    )
  })
  
  # 5. START PDF DEVICE (Specific to this tissue)
  pdf_path <- file.path(out_dir, paste0("Volcano_Plots_", e, "WaldComparison", ".pdf"))
  message("Saving plots to: ", pdf_path)
  
  pdf(pdf_path, width = 11.69, height = 8.27)
  walk(tissue_plots, print)
  dev.off() # Close device for this tissue
  
  # Store the dataframes in the master list if needed for later R analysis
  res_list[[e]] <- tissue_dfs
  
  message(e, " processing complete.")
}

# Check to make sure the names are correct (2 tissues x 5 timepoints == 10)
sum(map_int(res_list, length))

## Need to install R  4.5 to get the most up to date version of "apeglm" 1.32
## Using 'ashr' shrinkage method because 'apeglm' didnt converge so some timepoints

# Plotting Number DEG's up and down
# A.) WithinTissue (Parotoid and Liver) Line plot

# 1.) Generate the Dataframe
liver_degs <- get_deg_counts(res_list$liver)
liver_degs$timepoint <- factor(
  liver_degs$timepoint,
  levels = c("0h", "3h", "6h", "9h", "12h")
)


parotoid_degs <- get_deg_counts(res_list$parotoid)
parotoid_degs$timepoint <- factor(
  parotoid_degs$timepoint,
  levels = c("0h", "3h", "6h", "9h", "12h")
)

# Joining the tissues for the plotting

liver_degs <- liver_degs %>%
  mutate(tissue = "liver")

parotoid_degs <- parotoid_degs %>%
  mutate(tissue = "parotoid")

combined_degs <- bind_rows(parotoid_degs, liver_degs)


combined_degs$tissue <- factor(
  combined_degs$tissue,
  levels = c("parotoid", "liver")
)
write_csv(combined_degs, file = 'WithinTissuePairwiseWaldDEGsUpDown.csv')


plot1 <- ggplot(
  combined_degs,
  aes(
    x = timepoint,
    y = deg_counts_plot,
    color = tissue,
    linetype = reg_dir,
    group = interaction(tissue, reg_dir)
  )
) +
  
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  
  geom_hline(
    yintercept = 0,
    linetype = "solid",
    color = "black",
    linewidth = 1.0
  ) +
  
  labs(
    x = "Time Points",
    y = "Number of DEGs",
    color = "Tissue"
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
  
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 10)
  ) +
  
  scale_color_manual(
    values = c(
      "parotoid" = "#E65700",
      "liver"    = "#03548a"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "up" = "solid",
      "down" = "dotted"),
    guide = "none"      # removes the regulation legend
  ) 


print(plot1)

saveRDS(plot1, file = "WithinTissue_LinePlot.rds")


# Filtering for the LTR models
liver_lrt <- as.data.frame(results(dds_list$liver$lrt, alpha = 0.05))

liver_LRT_DEGs <- extract_LRT_DEGs(
  LRT_df = liver_lrt,
  tissue = "liver"
)


parotoid_lrt <- as.data.frame(results(dds_list$parotoid$lrt, alpha = 0.05))

parotoid_LRT_DEGs <- extract_LRT_DEGs(
  LRT_df = parotoid_lrt,
  tissue = "parotoid"
)


combined_LTR_degs <- bind_rows(parotoid_LRT_DEGs, liver_LRT_DEGs)


combined_LTR_degs$tissue <- factor(
  combined_LTR_degs$tissue,
  levels = c("parotoid", "liver")
)

# The plot
plot2 <- ggplot(
  combined_LTR_degs,
  aes(
    x = model,
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
    x = "Likelihood Ratio Test",
    y = "Number of DEGs",
    fill = "Tissue",
    pattern = "Tissue"
  ) +
  
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 10)   # increase number of ticks
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1.0),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),     # remove x tick labels
    axis.ticks.x = element_blank(),    # remove x tick marks
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  
  scale_fill_manual(
    values = c(
      "parotoid" = "#E65700",
      "liver"    = "#BA0600"
    )
  ) +
  
  scale_pattern_manual(
    values = c(
      "parotoid" = "pch",
      "liver"    = "pch"
    )
  )

print(plot2)

plot1 + theme(plot.margin = margin(40, 10, 10, 10))
plot2 + theme(plot.margin = margin(40, 10, 10, 10))

withinTissue <- plot_grid(
  plot1,
  plot2,
  labels = c("A", "B"),
  label_fontface = "bold",
  label_size = 30,
  hjust = -0.2,
  vjust = 1.1
)


ggsave(filename = "WithinTissueFigure.png", plot = withinTissue, dpi = 300, width = 20, height = 8.27, units = "in")


# To do's 
# Decide if we want a gene intersection for this figure (We don't)
# Decide if we like this plot (Line plot now)
# Make Tables from foramted for latex that can ve added to the supplementary
# Do we have to apply shrinkage for these results?? (apply shrinkage)