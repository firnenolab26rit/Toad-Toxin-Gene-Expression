# Parses the NX classified expression matrix into a table format
library(tidyverse)

expr_classified <- read_csv("csv/NX_tissue_matrix_classified.csv") 

# Tissue columns used in the UpSet
tissues <- c("muscle", "liver", "parotoid")


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

# Overall Intersections
# Group by your tissue columns and count the sizes of each intersection
intersection_sums_table <- upset_df %>%
  group_by(across(all_of(tissues))) %>%
  tally(name = "Count") %>%
  arrange(desc(Count)) # Sorts from largest to smallest intersection

# View the result
print(intersection_sums_table)

# Intersections
# 1. Identify which tissues are "TRUE" for each gene
# 2. Group by those tissues AND the specificity/detection classification
# 3. Count the members in each subgroup

# Specificity Intersections
S_intersection_stacked_counts <- upset_df %>%
  group_by(across(all_of(tissues)), specificity_class) %>%
  tally(name = "Count") %>%
  ungroup() %>%
  # Optional: Filter out rows where the gene wasn't detected in ANY tissue
  #filter(if_any(all_of(tissues), ~ . == TRUE)) %>%
  arrange(desc(Count)) %>% mutate(Percent = Count / sum(Count) * 100)

# View the table
print(S_intersection_stacked_counts)


# Detection Intersections
D_intersection_stacked_counts <- upset_df %>%
  group_by(across(all_of(tissues)), detection_class) %>%
  tally(name = "Count") %>%
  ungroup() %>%
  # Optional: Filter out rows where the gene wasn't detected in ANY tissue
  #filter(if_any(all_of(tissues), ~ . == TRUE)) %>%
  arrange(desc(Count)) %>% mutate(Percent = Count / sum(Count) * 100)

# View the table
print(D_intersection_stacked_counts)



# Table for "specific" intersection (detection_class == "Detected in single" , specificity_class == "Tissue enriched")
Specific_intersection_stacked_counts <- upset_df %>%
  filter(detection_class == "Detected in single" , specificity_class == "Tissue enriched") %>%
  group_by(across(all_of(tissues))) %>%
  tally(name = "Count") %>%
  ungroup() %>%
  # Optional: Filter out rows where the gene wasn't detected in ANY tissue
  #filter(if_any(all_of(tissues), ~ . == TRUE)) %>%
  arrange(desc(Count)) %>% mutate(Percent = Count / sum(Count) * 100)

# View the table
print(Specific_intersection_stacked_counts)

write_tsv(Specific_intersection_stacked_counts, file = "TissueAltasTables/SpecificIntersections.tsv")
write_tsv(S_intersection_stacked_counts, file = "TissueAltasTables/SpecificityIntersections.tsv")
write_tsv(D_intersection_stacked_counts, file = "TissueAltasTables/DetectionIntersections.tsv")
