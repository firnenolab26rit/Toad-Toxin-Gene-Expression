library(tidyverse)
library(stringr)

# 1. Read in LTR DEseq2 expr ~ tissue 
Deseq_df <- read_csv(file = "../DEseq/BetweenTissue/DE_results/BetweenTissue_LTR/BetweenTissue_ParotoidLiver_LTR.csv")
# Ensure the first column is named "Gene" to match our downstream joins
if(!"Gene" %in% colnames(Deseq_df)) colnames(Deseq_df)[1] <- "Gene"

# 2. Extract Modules of interest from the lm model 
#Tissue_Models <- read_csv(file = "../WGCNA/All_Modules_Tissue_LM_Results.csv")

# Pull significant modules and strip the "ME" prefix (e.g., "MEblue" -> "blue")
# Modules_Significant <- Tissue_Models %>% 
#   filter(p.value < 0.05) %>% 
#   pull(module) %>%
#   str_remove("^ME")

# 3. Identify genes belonging to the modules
# Using geneInfo0 gives us the hard moduleColor assignments
#geneInfo <- read.csv("../WGCNA/MG_AA_wgcna_April4_2025_geneInfo0.csv", row.names = 1)

# genes_of_interest <- geneInfo %>%
#   pull(Gene)

# 4. Subset the Deseq_df for those genes (Have to subset some we have the same genes bewtween the analysis)
# I think the prefiltering was slightly different. 
# Deseq_df_subset <- Deseq_df %>%
#   filter(Gene %in% genes_of_interest)

# 5. Save the subset to a temporary or intermediate file
# callr::rscript requires a file path string, it cannot read a live R environment object
subset_input_file <- "./results_NOSubset/Tissue/Deseq_KeggPathway.csv"
write_csv(Deseq_df, file = subset_input_file)

# 6. Run GSEA 
# Corrected callr syntax
callr::rscript(
  "../GSEA/GSEA.R",
  c(
    "-i", subset_input_file,
    "-g", "../Annotation/Enrichment_Inputs_subset/KEGG/term2gene_kegg_pathway.tsv",
    "-n", "../Annotation/Enrichment_Inputs_subset/KEGG/term2name_kegg_pathway.tsv",
    "-r", "log2FoldChange",
    "-o", "results_NOSubset/Tissue"
  )
)

