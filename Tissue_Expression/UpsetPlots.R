library(dplyr)
library(tidyr)
library(ComplexUpset)
library(ggplot2)
library(readr)

#############################################################################################################
## This script makes upset plots for the intersections based on tissue specificity and tissue detection
#############################################################################################################
# Read in the data from generated from the tissue altas (atlas.R)
expr_classified <- read_csv("csv/NX_tissue_matrix_classified.csv") 

# -----------------------------
# Define parameters
# -----------------------------
# tissues of interest for intersection
tissues <- c("muscle", "liver", "parotoid")
# Colors c("#B1B3B3FF", "#888B8DFF","#E6D400", "#E69F00", "#E65700")

# Detection colors
my_colors <- rev(c("#B1B3B3FF", "#888B8DFF","#E6D400", "#E65700"))
# Specicifity colors
my_colors2 <- rev(c(
  "Not detected"            = "#B1B3B3FF",
  "Low tissue specificity"  = "#888B8DFF",
  "Tissue enhanced"         = "#E69F00",
  "Group enriched"          = "#E6D400",
  "Tissue enriched"         = "#E65700"
))


# -----------------------------
# Long format + detection flag
# -----------------------------
# make long format and te
expr_long <- expr_classified %>%
  pivot_longer(
    cols = all_of(tissues),
    names_to = "tissue",
    values_to = "NX"
  ) %>%
  mutate(detected = NX >= 1)

# -----------------------------
# Gene–tissue membership table
# -----------------------------
gene_sets <- expr_long %>%
  select(gene_id, tissue, detected) %>%
  pivot_wider(
    names_from = tissue,
    values_from = detected,
    values_fill = FALSE
  )

# -----------------------------
# Add classification columns
# -----------------------------
upset_df <- gene_sets %>%
  left_join(
    expr_classified %>% 
      select(gene_id, specificity_class, detection_class),
    by = "gene_id"
  )

# Remove the 0 groups 
#upset_df <- upset_df %>% filter(specificity_class != "Tissue enhanced")

# Specify Detection and Specificity Levels 
s <- rev(c("Not detected", "Low tissue specificity", "Tissue enhanced", "Group enriched", "Tissue enriched"))
d <- rev(c("Not detected", "Detected in all", "Detected in many", "Detected in single"))

upset_df <- upset_df %>% mutate(
  specificity_class = factor(specificity_class, levels = s), 
  detection_class = factor(detection_class, levels = d)
)

# ----------------------------
# UpSet plot 1: Specificity class
# -----------------------------
# 1. Plot for Specificity
plot1 <- upset(
  upset_df,
  tissues,
  name = "Tissue Intersections",
  width_ratio = 0.1,
  
  # Matrix styling
  matrix = intersection_matrix(
    segment = geom_segment(linewidth = 1.2)
  ),
  
  # Annotation styling
  base_annotations = list(
    'Intersection size' = intersection_size( 
      counts = FALSE,
      mapping = aes(fill = specificity_class),
      color = "black",
      linewidth = 0.6
    ) + 
      scale_fill_manual(values = my_colors2, name = "Specificity") + 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "bold", size = 14), # Bold Y axis of the bar plot
        axis.title.y = element_text(face = "bold", size = 14)
        
      )
  ),
  set_sizes = FALSE
) + 
  # Global theme for the matrix and set labels
  theme(
    axis.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14), 
    axis.title.x = element_text(face = "bold", size = 14)
  )

print(plot1)
# -----------------------------
# UpSet plot 2: Detection class
# -----------------------------
# 2. Plot for Detection
plot2 <- upset(
  upset_df,
  tissues,
  name = "Tissue Intersections",
  width_ratio = 0.1,
  
  # Matrix styling
  matrix = intersection_matrix(
    segment = geom_segment(linewidth = 1.2)
  ),
  
  # Annotation styling
  base_annotations = list(
    'Intersection size' = intersection_size( 
      counts = FALSE,
      mapping = aes(fill = detection_class),
      color = "black",
      linewidth = 0.6
    ) + 
      scale_fill_manual(values = my_colors, name = "Detection") + 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "bold", size = 14), # Bold Y axis of the bar plot
        axis.title.y = element_text(face = "bold", size = 14)
        
      )
  ),
  set_sizes = FALSE
) + 
  # Global theme for the matrix and set labels
  theme(
    axis.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 14), 
    axis.title.x = element_text(face = "bold", size = 14)
  )

print(plot2)

ggsave(filename = "Upset_Specificity.png", plot = plot1, dpi = 300, width = 11.69, height = 8.27, units = "in")
ggsave(filename = "Upset_Detection.png", plot = plot2, dpi = 300, width = 11.69, height = 8.27, units = "in")

# Venn Diagram 
# Make a class that is interest (Ex. detect single and group enriched) and then plot that as a venn diagram across 
# the three tissues

# can do that with multiple categories or just a single like (group enriched)
