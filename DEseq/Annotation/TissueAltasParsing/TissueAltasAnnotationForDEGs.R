library(tidyverse)
library(ggplot2)

TA <- read_csv("~/Downloads/TOAD_ANALYSIS/MultiTissue_RNAseq/Tissue_Expression/csv/NX_tissue_matrix_classified.csv") %>%
  select(-1)


# Need to get a list of all unique DEG'S we found in the comparisons

dir_paths <- c("~/Downloads/TOAD_ANALYSIS/MultiTissue_RNAseq/DEseq/WithinTissue/DE_results/parotoid",
              "~/Downloads/TOAD_ANALYSIS/MultiTissue_RNAseq/DEseq/WithinTissue/DE_results/liver",
              "~/Downloads/TOAD_ANALYSIS/MultiTissue_RNAseq/DEseq/BetweenTissue/DE_results/BetweenTissue")
# Function (Accepts directory path)
# Reads in all *.csv extracts specific column1 to vector and then does unique()

get_unique_column_values <- function(dir_paths, column_name) {
  dir_paths |>
    map(~list.files(.x, pattern = "\\.csv$", full.names = TRUE)) |>
    unlist() |>
    map_dfr(read_csv, show_col_types = FALSE) |>
    pull({{column_name}}) |>
    unique()
}

genes <- get_unique_column_values(dir_paths, 1)



# Join with the altlas to get there status (Actually a filter)
TA_DEG <- TA %>% filter(gene_id %in% genes)

# Some general stats
test <- as.data.frame(table(TA_DEG$specificity_class))
table(TA_DEG$detection_class)

# Can use mutate and count to make percentages with group by


# Filter Out DEG's that are in single tissue and enriched
specific <- TA_DEG %>% filter(specificity_class == "Tissue enriched", detection_class == "Detected in single")

# ~ 16% of DEGs are specific to a single tissue 
specific[specific < 1] <- NA


# 1 = expressed, 0 = not expressed
presence <- !is.na(specific)

colnames(presence)


muscle <- which(presence[,2])
liver <- which(presence[,3])
parotoid <- which(presence[,4])

length(muscle)  # tissue 1
length(liver)  # tissue 2
length(parotoid)  # tissue 3

length(intersect(muscle, liver))
length(intersect(muscle, parotoid))
length(intersect(liver, parotoid))
length(Reduce(intersect, list(muscle, liver, parotoid)))


overlap_table <- data.frame(
  Category = c(
    "Muscle only",
    "Liver only",
    "Parotoid only",
    "Muscle ∩ Liver",
    "Muscle ∩ Parotoid",
    "Liver ∩ Parotoid",
    "Muscle ∩ Liver ∩ Parotoid",
    "N Genes Total"
  ),
  Count = c(
    length(setdiff(muscle, union(liver, parotoid))),       # Muscle only
    length(setdiff(liver, union(muscle, parotoid))),       # Liver only
    length(setdiff(parotoid, union(muscle, liver))),       # Parotoid only
    length(intersect(muscle, liver)) - length(Reduce(intersect, list(muscle, liver, parotoid))),
    length(intersect(muscle, parotoid)) - length(Reduce(intersect, list(muscle, liver, parotoid))),
    length(intersect(liver, parotoid)) - length(Reduce(intersect, list(muscle, liver, parotoid))),
    length(Reduce(intersect, list(muscle, liver, parotoid))), 
    nrow(specific)
  )
)


total  = nrow(specific)

overlap_table <- overlap_table %>% mutate(Percent = round(Count / total * 100))

# Keep only the "only" categories
exclusive <- overlap_table %>% 
  filter(Category %in% c("Muscle only", "Liver only", "Parotoid only")) %>%
  mutate(Category = factor(Category, levels = c("Muscle only", "Liver only", "Parotoid only")))

plot_title <- bquote("Tissue Enriched" ~ union() ~ "Single Tissue")


# Barplot
plot3 <- ggplot(exclusive, aes(x = Category, y = Count, fill = Category)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.6) +
  geom_text(aes(label = Count), vjust = 3.0, size = 10, fontface = "bold") +
  scale_fill_manual(values = c("#888B8DFF", "#03548a", "#E65700")) +
  labs(
    title = plot_title,
    y = "Intersection size",
    x = ""
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(face = "bold", size = 14),   # <-- set axis text size here
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

plot3



ggsave(filename = "FigS2Tissue_SpecificDEGs.pdf", plot = plot3, dpi = 600, width = 11.69, height = 8.27, units = "in", device = cairo_pdf)
ggsave(filename = "FigS2Tissue_SpecificDEGs.jpeg", plot = plot3, dpi = 600, width = 11.69, height = 8.27, units = "in")
ggsave(filename = "FigS2Tissue_SpecificDEGs.svg", plot = plot3, dpi = 600, width = 11.69, height = 8.27, units = "in")


specific_long <- specific %>%
  pivot_longer(
    cols = c(muscle, liver, parotoid),
    names_to = "tissue",
    values_to = "NX"
  ) %>%
  filter(!is.na(NX)) %>%
  select(gene_id, tissue, NX)

write_tsv(specific_long, "specific_tissue_NX.tsv")

# To get the specific tissues which they are in 
specific_long_tissue <- specific_long %>% group_by(gene_id) %>% filter(NX == max(NX))

# Filter out DEG's tissue enriched


# Filter out DEG's tissue enhanced

## Think of comparisons that are useful


## Get annotations for DEG'S filtered and can do ORA/GESA
specific_long <- read_tsv("specific_tissue_NX.tsv")
colnames(specific_long)[1] <- "trinity_id"


anno <- read_tsv("../A.Americanus._Transcriptome_Annotations_new.tsv")
# Readname column 1 to eviene_id 
colnames(anno)[1] <- "evigene_id"



gm <- read_table("../blast_search/gene_map.tsv", col_names = c("trinity_id", "evigene_id"))

# Lets see how many of our seqeuences had annotations 
nrow(anno) / nrow(gm) * 100 

joined <- left_join(anno, gm, by = "evigene_id")

tissue_specific_DEGs <- joined %>% filter(trinity_id %in% specific_long$trinity_id) %>% full_join(specific_long, by = "trinity_id")

tissue_specific_DEGs %>% group_by(tissue) %>% summarise(Count = n())


# Lets see the breakdown of groups
tissue_specific_DEGs_summary <- tissue_specific_DEGs %>% group_by(tissue, Description) %>% summarise(n())




### Example looking at specific timepoint 
# p_timepoint <- read_csv("~/Downloads/TOAD_ANALYSIS/MultiTissue_RNAseq/DEseq/WithinTissue/DE_results/parotoid/DEG_parotoid_6h_vs_con.csv")
# colnames(p_timepoint)[1] <- "trinity_id"
# 
# p_timepoint_annotations <- joined %>% filter(trinity_id %in% p_timepoint$trinity_id)
# 
# test <- inner_join(p_timepoint_annotations, p_timepoint, by = "trinity_id")