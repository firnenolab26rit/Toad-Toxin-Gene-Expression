# See what annotatations of intereest are tissue specific. (single tissue & Tissue enriched)

library(tidyverse)

# Load in the tissue altas classifications
expr_classified <- read_csv("../../Tissue_Expression/csv/NX_tissue_matrix_classified.csv") 

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

colnames(expr_long)[2] <- "trinity_id"


# Load in the eggnog
# ================================
# 1. Load annotation + filter genes
# ================================
eggnog <- read_tsv("../../Annotation/Gene_Level_Annotations_Wide.tsv")

# All search pattern
pattern <- c("[Pp]450", 
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

# Just to find the number of cyps
#pattern <- c("[Pp]450")
search_cols <- c("PFAMs", "Description", "Preferred_name")

# Combine patterns with regex OR
pattern_regex <- paste(pattern, collapse = "|")

mine <- eggnog %>%
  filter(
    if_any(
      all_of(search_cols),
      ~ str_detect(as.character(.x), regex(pattern_regex, ignore_case = TRUE))
    )
  )


# Join the atlas data and the annotation data
joined_anno <- inner_join(mine, expr_long, by = "trinity_id")



# Get the numbers by filtering and summarizing
# Muscle specific 
joined_anno %>% filter(tissue == "muscle", 
                       detected == TRUE,
                       specificity_class == "Tissue enriched", 
                       detection_class == "Detected in single") %>% nrow()

# Liver specific 
joined_anno %>% filter(tissue == "liver", 
                       detected == TRUE,
                       specificity_class == "Tissue enriched", 
                       detection_class == "Detected in single") %>% nrow()

# Parotoid specific 
joined_anno %>% filter(tissue == "parotoid", 
                       detected == TRUE,
                       specificity_class == "Tissue enriched", 
                       detection_class == "Detected in single") %>% nrow()




