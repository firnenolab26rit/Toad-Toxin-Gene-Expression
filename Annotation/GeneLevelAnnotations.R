# Puprose of this is script is join annotation data per gene with a gene. 
# ==============================================================================
# Functional Annotation Aggregation: Transcript to Gene Level
# ==============================================================================
library(tidyverse)
# orginal
eggnog="../DEseq/Annotation/A_Americanus_Transcriptome_Annotations_new.tsv"

# Subseted
eggnog="../DEseq/Annotation/Annotations_subset.tsv"

# 1. Import Data
cat("Importing annotations and gene map...\n")
anno <- read_delim(eggnog, 
                   show_col_types = FALSE) %>% 
  dplyr::rename(evigene_id = 1)

gm <- read_delim("../DEseq/Annotation/blast_search/gene_map.tsv", 
                 col_names = c("trinity_id", "evigene_id"), 
                 show_col_types = FALSE, col_types = "cc")

# Join and establish baseline dataset
anno_data <- dplyr::left_join(anno, gm, by = "evigene_id")

# 2. Define Columns of Interest
focus_cols <- c(
  "trinity_id", "COG_category", "Description", "Preferred_name",
  "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
  "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"
)

# Subset data and remove any unmapped transcripts
anno_subset <- anno_data %>%
  select(any_of(focus_cols)) %>%
  filter(!is.na(trinity_id))

# ==============================================================================
# 3. Modular Aggregation Function
# ==============================================================================
# This function takes a vector of strings (from multiple transcripts), 
# splits any comma-separated values, removes NAs/blanks, and returns a 
# single string of unique, collapsed values for the parent gene.
collapse_unique <- function(x) {
  # Remove NAs, empty strings, and eggNOG's default "-" empty character
  x <- x[!is.na(x) & x != "" & x != "-"]
  
  if (length(x) == 0) return(NA_character_)
  
  # Split strings by comma, flatten to a single list, trim whitespace
  split_vals <- unlist(strsplit(as.character(x), split = ","))
  clean_vals <- trimws(split_vals)
  
  # Extract unique terms and collapse back to a single string
  unique_vals <- unique(clean_vals[clean_vals != ""])
  
  if (length(unique_vals) == 0) return(NA_character_)
  
  return(paste(unique_vals, collapse = ","))
}

# ==============================================================================
# 4. Execute Gene-Level Summarization
# ==============================================================================
cat("Aggregating transcript annotations to gene level...\n")

gene_level_annotations <- anno_subset %>%
  group_by(trinity_id) %>%
  # Apply the collapse_unique function to all functional columns
  summarise(across(everything(), collapse_unique), .groups = 'drop')

# 5. Export Wide Format (For DEseq2 / WGCNA merging)
write_tsv(gene_level_annotations, "Gene_Level_Annotations_Wide.tsv")
cat("Wide format saved. Rows:", nrow(gene_level_annotations), "\n")

# ==============================================================================
# Functional Annotation: Prep for clusterProfiler (ORA / GSEA)
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(AnnotationDbi)
  library(GO.db)
})

# Assuming 'gene_level_annotations' is already in your environment from the previous script
# If not: gene_level_annotations <- read_tsv("Gene_Level_Annotations_Wide.tsv")

# ==============================================================================
# 1. Universal Long-Format & Dictionary Join Function
# ==============================================================================
# This function creates the exact 2-column dataframes clusterProfiler needs:
# term2gene: [TermID, GeneID]
# term2name: [TermID, TermDescription]
build_enrichment_tables <- function(df, target_col, ref_dict = NULL, prefix_to_remove = NULL, deduplicate_kegg = FALSE) {
  
  # 1. Extract and Pivot Long
  long_df <- df %>%
    dplyr::select(trinity_id, all_of(target_col)) %>%
    tidyr::drop_na() %>%
    tidyr::separate_rows(all_of(target_col), sep = ",") %>%
    dplyr::mutate(!!target_col := trimws(get(target_col)))
  
  # 2. Clean Prefixes (e.g., stripping 'ko:' from 'ko:K00001')
  if (!is.null(prefix_to_remove)) {
    long_df <- long_df %>%
      dplyr::mutate(!!target_col := stringr::str_remove(get(target_col), prefix_to_remove))
  }
  
  # 3. Optional KEGG Deduplication (keeps koXXXXX, drops mapXXXXX)
  if (deduplicate_kegg) {
    long_df <- long_df %>%
      dplyr::mutate(temp_num_id = stringr::str_extract(get(target_col), "\\d+")) %>%
      dplyr::group_by(trinity_id, temp_num_id) %>%
      dplyr::arrange(get(target_col), .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(-temp_num_id)
  }
  
  # 4. Create term2gene [Term, Gene]
  term2gene <- long_df %>%
    dplyr::select(all_of(target_col), trinity_id) %>%
    dplyr::distinct()
  
  # 5. Create term2name [Term, Name]
  term2name <- NULL
  if (!is.null(ref_dict)) {
    ref_colnames <- colnames(ref_dict)
    
    term2name <- long_df %>%
      dplyr::select(all_of(target_col)) %>%
      dplyr::distinct() %>%
      dplyr::left_join(ref_dict, by = setNames(ref_colnames[1], target_col)) %>%
      tidyr::drop_na()
  }
  
  return(list(term2gene = term2gene, term2name = term2name))
}

# ==============================================================================
# 2. Gene Ontology (GO) Processing using GO.db
# ==============================================================================
message("Processing GO terms...")

# Build reference dictionary from GO.db
go_meta <- AnnotationDbi::select(GO.db, 
                                 keys = keys(GO.db), 
                                 columns = c("TERM", "ONTOLOGY"), 
                                 keytype = "GOID") %>%
  dplyr::rename(ID = GOID, Name = TERM, Ontology = ONTOLOGY)

# Extract GO term2gene (No prefix removal needed for GOs)
go_tables <- build_enrichment_tables(gene_level_annotations, "GOs", ref_dict = go_meta)

# Split term2gene by Ontology (BP, MF, CC)
go_term2gene_annotated <- go_tables$term2gene %>%
  dplyr::left_join(go_meta, by = c("GOs" = "ID"))

dir.create("Enrichment_Inputs_subset/GO", recursive = TRUE, showWarnings = FALSE)

for (ont in c("BP", "MF", "CC")) {
  # Subset for ontology and format strictly as [Term, Gene]
  t2g_sub <- go_term2gene_annotated %>%
    dplyr::filter(Ontology == ont) %>%
    dplyr::select(GOs, trinity_id)
  
  # Format strictly as [Term, Name]
  t2n_sub <- go_meta %>%
    dplyr::filter(Ontology == ont) %>%
    dplyr::select(ID, Name)
  
  write_tsv(t2g_sub, paste0("Enrichment_Inputs_subset/GO/term2gene_GO_", ont, ".tsv"))
  write_tsv(t2n_sub, paste0("Enrichment_Inputs_subset/GO/term2name_GO_", ont, ".tsv"))
}

# ==============================================================================
# 3. KEGG Reference Download & Management
# ==============================================================================
get_kegg_references <- function(cache_dir = "kegg_cache") {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define file paths
  pathway_file <- file.path(cache_dir, "kegg_pathway.tsv")
  ko_file <- file.path(cache_dir, "kegg_ko.tsv")
  module_file <- file.path(cache_dir, "kegg_module.tsv")
  
  # Download if missing
  if (!file.exists(pathway_file)) {
    message("Downloading KEGG references...")
    
    # Base pathways (map)
    map_path <- read.delim("https://rest.kegg.jp/list/pathway", header = FALSE)
    # KO specific pathways (ko) - As requested!
    ko_path <- read.delim("https://rest.kegg.jp/list/pathway/ko", header = FALSE)
    
    kegg_pathway <- bind_rows(map_path, ko_path) %>%
      setNames(c("id", "name")) %>%
      mutate(id = gsub("path:", "", id),
             name = str_remove(name, " - .*")) %>%
      distinct()
    
    kegg_ko <- read.delim("https://rest.kegg.jp/list/ko", header = FALSE) %>%
      setNames(c("id", "name")) %>%
      mutate(id = gsub("ko:", "", id),
             name = str_remove(name, ";.*"))
    
    kegg_module <- read.delim("https://rest.kegg.jp/list/module", header = FALSE) %>%
      setNames(c("id", "name")) %>%
      mutate(id = gsub("md:", "", id))
    
    write_tsv(kegg_pathway, pathway_file)
    write_tsv(kegg_ko, ko_file)
    write_tsv(kegg_module, module_file)
  }
  
  # Load and return
  list(
    pathway = read_tsv(pathway_file, show_col_types = FALSE),
    ko = read_tsv(ko_file, show_col_types = FALSE),
    module = read_tsv(module_file, show_col_types = FALSE)
  )
}

kegg_refs <- get_kegg_references()

# ==============================================================================
# 4. KEGG Processing (Pathways, Modules, KOs)
# ==============================================================================
message("Processing KEGG annotations...")
dir.create("Enrichment_Inputs_subset/KEGG", recursive = TRUE, showWarnings = FALSE)

# --- KEGG Pathways ---
# Eggnog usually outputs pathways like "ko00010,map00010". No prefix removal needed.
pathway_tables <- build_enrichment_tables(gene_level_annotations, "KEGG_Pathway", ref_dict = kegg_refs$pathway, deduplicate_kegg = TRUE)
write_tsv(pathway_tables$term2gene, "Enrichment_Inputs_subset/KEGG/term2gene_kegg_pathway.tsv")
write_tsv(pathway_tables$term2name, "Enrichment_Inputs_subset/KEGG/term2name_kegg_pathway.tsv")

# --- KEGG Modules ---
# Eggnog usually outputs modules with prefixes or cleanly. We'll strip "md:" just in case.
module_tables <- build_enrichment_tables(gene_level_annotations, "KEGG_Module", ref_dict = kegg_refs$module, prefix_to_remove = "^md:")
write_tsv(module_tables$term2gene, "Enrichment_Inputs_subset/KEGG/term2gene_kegg_module.tsv")
write_tsv(module_tables$term2name, "Enrichment_Inputs_subset/KEGG/term2name_kegg_module.tsv")

# --- KEGG KO ---
# Strip the 'ko:' prefix that eggNOG applies (e.g., 'ko:K00001' -> 'K00001')
ko_tables <- build_enrichment_tables(gene_level_annotations, "KEGG_ko", ref_dict = kegg_refs$ko, prefix_to_remove = "^ko:")
write_tsv(ko_tables$term2gene, "Enrichment_Inputs_subset/KEGG/term2gene_kegg_ko.tsv")
write_tsv(ko_tables$term2name, "Enrichment_Inputs_subset/KEGG/term2name_kegg_ko.tsv")

message("Complete! All term2gene and term2name tables are ready for clusterProfiler.")



# End of Script 
# Extra old logic










































# # ==============================================================================
# # OPTIONAL: 6. Long Format for Enrichment Analysis (e.g., clusterProfiler)
# # ==============================================================================
# # Enrichment tools usually require a 2-column data frame: [Term, Gene]
# # This cleanly pivots a specific target column (like GOs) into that format.
# 
# extract_long_format <- function(df, target_col) {
#   df %>%
#     select(trinity_id, all_of(target_col)) %>%
#     drop_na() %>%
#     # Separate the collapsed strings into individual rows
#     separate_rows(all_of(target_col), sep = ",") %>%
#     # Reorder columns for standard term2gene format: [Pathway/GO, Gene]
#     select(all_of(target_col), trinity_id)
# }
# 
# # Example: Extract GO terms for enrichment
# term2gene_GO <- extract_long_format(gene_level_annotations, "GOs")
# write_tsv(term2gene_GO, "Term2Gene_GO_Long.tsv")
# 
# term2gene_KEGG_pathway <- extract_long_format(gene_level_annotations, "KEGG_Pathway")
# 
# 
# term2gene_COG_catergory <- extract_long_format(gene_level_annotations, "COG_category")
# 
# 
# # Make all this generalizable so I can do the same with other columsn like BRITE, or Kegg_rclass if I want to.
# ## Need to get the Descriptions (GO, KEGG)
# 
# # I want the columns
#  "trinity_id", "Gos", "Name"
#  
# # So I can make a term2gene and term2name by just selecting different columns (ie. "GeneID", "TERM", "NAME")
# # I want to make this generalizabke to use with any of my annotation columns (I.e KEGG_pathway, KEGG_Module)
# 
# # Go 
# # Using the GO.db much faster and should be up to date enough 
# # Extract a dataframe of ALL GO terms, their names, and their ontology (BP, MF, CC)
# go_meta <- AnnotationDbi::select(GO.db, 
#                                  keys = keys(GO.db), 
#                                  columns = c("TERM", "ONTOLOGY"), 
#                                  keytype = "GOID")
# term2gene_GO <- go_meta %>% full_join()
# 
# # Split and Write to a results GO/
# ont_BP <- go_meta[go_meta$ONTOLOGY == "BP", ]
# ont_MF <- go_meta[go_meta$ONTOLOGY == "MF", ]
# ont_CC <- go_meta[go_meta$ONTOLOGY == "CC", ]
# 
# # Kegg
# # ============================================================
# # KEGG Annotation Utilities
# # - Pathways (mapXXXXX) and (koXXXXXX) need to add (https://rest.kegg.jp/list/pathway/ko) for ko
# # - KOs (KXXXXX)
# # - Modules (MXXXXX)
# # Designed for eggNOG-mapper outputs
# # ============================================================
# 
# suppressPackageStartupMessages({
#   library(tidyverse)
# })
# 
# # ============================================================
# # 1. Download KEGG reference tables
# # ============================================================
# 
# download_kegg_data <- function(cache_dir = NULL, overwrite = FALSE) {
#   
#   if (!is.null(cache_dir)) {
#     dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
#   }
#   
#   # Add https://rest.kegg.jp/list/pathway/ko as well to get full annotations
#   message("Downloading KEGG pathway data...")
#   kegg_pathway <- read.delim("https://rest.kegg.jp/list/pathway", header = FALSE) %>%
#     setNames(c("kegg_id", "description")) %>%
#     mutate(
#       kegg_id = gsub("path:", "", kegg_id),
#       description = str_remove(description, " - .*")
#     )
#   
#   message("Downloading KEGG KO data...")
#   kegg_ko <- read.delim("https://rest.kegg.jp/list/ko", header = FALSE) %>%
#     setNames(c("ko_id", "description")) %>%
#     mutate(
#       ko_id = gsub("ko:", "", ko_id),
#       description = str_remove(description, ";.*")
#     )
#   
#   message("Downloading KEGG module data...")
#   kegg_module <- read.delim("https://rest.kegg.jp/list/module", header = FALSE) %>%
#     setNames(c("module_id", "description")) %>%
#     mutate(
#       module_id = gsub("md:", "", module_id)
#     )
#   
#   # Save if caching enabled
#   if (!is.null(cache_dir)) {
#     write_tsv(kegg_pathway, file.path(cache_dir, "kegg_pathway.tsv"))
#     write_tsv(kegg_ko, file.path(cache_dir, "kegg_ko.tsv"))
#     write_tsv(kegg_module, file.path(cache_dir, "kegg_module.tsv"))
#   }
#   
#   return(list(
#     pathway = kegg_pathway,
#     ko = kegg_ko,
#     module = kegg_module
#   ))
# }
# 
# # ============================================================
# # 2. Load cached KEGG data
# # ============================================================
# 
# load_kegg_data <- function(cache_dir) {
#   list(
#     pathway = read_tsv(file.path(cache_dir, "kegg_pathway.tsv"), show_col_types = FALSE),
#     ko = read_tsv(file.path(cache_dir, "kegg_ko.tsv"), show_col_types = FALSE),
#     module = read_tsv(file.path(cache_dir, "kegg_module.tsv"), show_col_types = FALSE)
#   )
# }
# 
# # ============================================================
# # 3. Annotation functions (vector input)
# # ============================================================
# 
# annotate_kegg_pathway <- function(ids, ref) {
#   tibble(kegg_id = ids) %>%
#     left_join(ref, by = "kegg_id")
# }
# 
# annotate_kegg_ko <- function(ids, ref) {
#   tibble(ko_id = ids) %>%
#     left_join(ref, by = "ko_id")
# }
# 
# annotate_kegg_module <- function(ids, ref) {
#   tibble(module_id = ids) %>%
#     left_join(ref, by = "module_id")
# }
# 
# 
# # ============================================================
# # 3. Annotation functions (column Input)
# # ============================================================
# 
# annotate_kegg_column <- function(df, column, ref, ref_key) {
#   
#   col_sym <- rlang::sym(column)
#   
#   df %>%
#     tidyr::separate_longer_delim(!!col_sym, delim = ",") %>%
#     dplyr::select(!!col_sym) %>%
#     dplyr::rename(!!ref_key := !!col_sym) %>%
#     dplyr::left_join(ref, by = ref_key) %>%
#     tidyr::drop_na()
# }
# 
# 
# # ============================================================
# # 6. Example usage (uncomment to run)
# # ============================================================
# 
# # ---- Step 1: Download or load KEGG ----
# # refs <- download_kegg_data(cache_dir = "kegg_cache")
# refs <- load_kegg_data("kegg_cache")
# 
# # ---- Step 2: Example vectors ----
# # annotate_kegg_ko(c("K00001", "K00844"), refs$ko)
# 
# # ---- Step 3: Example eggNOG table ----
# eggnog <- read_tsv("../../A.Americanus._Transcriptome_Annotations_new.tsv")
# 
# 
# annotated_module <- annotate_kegg_column(
#   df = eggnog,
#   column = "KEGG_Module",
#   ref = refs$module,
#   ref_key = "module_id"
# )
# 
# 
# annotated_pathway <- annotate_kegg_column(
#   df = eggnog,
#   column = "KEGG_Pathway",
#   ref = refs$pathway,
#   ref_key = "kegg_id"
# )
# 
# 
# df = eggnog
# column = "KEGG_ko"
# ref = refs$ko
# ref_key = "ko_id"
# 
# 
# 
# # Had to mutate the KO: out of the id
# col_sym <- rlang::sym(column)
# 
# annotated_ko <- df %>%
#   tidyr::separate_longer_delim(!!col_sym, delim = ",") %>%
#   dplyr::select(!!col_sym) %>%
#   dplyr::mutate(!!col_sym := stringr::str_remove(!!col_sym, "^ko:")) %>%
#   dplyr::rename(!!ref_key := !!col_sym)%>%
#   dplyr::left_join(ref, by = ref_key) %>%
#   tidyr::drop_na()
# 
# 
# 
# # ---- Step 4: Write outputs ----
# write_tsv(annotated_ko, "eggnog_ko_annotated.tsv")
# write_tsv(annotated_pathway, "eggnog_pathway_annotated.tsv")
# write_tsv(annotated_module, "eggnog_module_annotated.tsv")




