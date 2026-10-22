#!/usr/bin/env Rscript

# ============================================================
# KEGG Annotation Utilities
# - Pathways (mapXXXXX)
# - KOs (KXXXXX)
# - Modules (MXXXXX)
# Designed for eggNOG-mapper outputs
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================
# 1. Download KEGG reference tables
# ============================================================

download_kegg_data <- function(cache_dir = NULL, overwrite = FALSE) {
  
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  message("Downloading KEGG pathway data...")
  kegg_pathway <- read.delim("https://rest.kegg.jp/list/pathway", header = FALSE) %>%
    setNames(c("kegg_id", "description")) %>%
    mutate(
      kegg_id = gsub("path:", "", kegg_id),
      description = str_remove(description, " - .*")
    )
  
  message("Downloading KEGG KO data...")
  kegg_ko <- read.delim("https://rest.kegg.jp/list/ko", header = FALSE) %>%
    setNames(c("ko_id", "description")) %>%
    mutate(
      ko_id = gsub("ko:", "", ko_id),
      description = str_remove(description, ";.*")
    )
  
  message("Downloading KEGG module data...")
  kegg_module <- read.delim("https://rest.kegg.jp/list/module", header = FALSE) %>%
    setNames(c("module_id", "description")) %>%
    mutate(
      module_id = gsub("md:", "", module_id)
    )
  
  # Save if caching enabled
  if (!is.null(cache_dir)) {
    write_tsv(kegg_pathway, file.path(cache_dir, "kegg_pathway.tsv"))
    write_tsv(kegg_ko, file.path(cache_dir, "kegg_ko.tsv"))
    write_tsv(kegg_module, file.path(cache_dir, "kegg_module.tsv"))
  }
  
  return(list(
    pathway = kegg_pathway,
    ko = kegg_ko,
    module = kegg_module
  ))
}

# ============================================================
# 2. Load cached KEGG data
# ============================================================

load_kegg_data <- function(cache_dir) {
  list(
    pathway = read_tsv(file.path(cache_dir, "kegg_pathway.tsv"), show_col_types = FALSE),
    ko = read_tsv(file.path(cache_dir, "kegg_ko.tsv"), show_col_types = FALSE),
    module = read_tsv(file.path(cache_dir, "kegg_module.tsv"), show_col_types = FALSE)
  )
}

# ============================================================
# 3. Annotation functions (vector input)
# ============================================================

annotate_kegg_pathway <- function(ids, ref) {
  tibble(kegg_id = ids) %>%
    left_join(ref, by = "kegg_id")
}

annotate_kegg_ko <- function(ids, ref) {
  tibble(ko_id = ids) %>%
    left_join(ref, by = "ko_id")
}

annotate_kegg_module <- function(ids, ref) {
  tibble(module_id = ids) %>%
    left_join(ref, by = "module_id")
}


# ============================================================
# 3. Annotation functions (column Input)
# ============================================================

annotate_kegg_column <- function(df, column, ref, ref_key) {
  
  col_sym <- rlang::sym(column)
  
  df %>%
    tidyr::separate_longer_delim(!!col_sym, delim = ",") %>%
    dplyr::select(!!col_sym) %>%
    dplyr::rename(!!ref_key := !!col_sym) %>%
    dplyr::left_join(ref, by = ref_key) %>%
    tidyr::drop_na()
}


# ============================================================
# 6. Example usage (uncomment to run)
# ============================================================

# ---- Step 1: Download or load KEGG ----
# refs <- download_kegg_data(cache_dir = "kegg_cache")
refs <- load_kegg_data("kegg_cache")

# ---- Step 2: Example vectors ----
# annotate_kegg_ko(c("K00001", "K00844"), refs$ko)

# ---- Step 3: Example eggNOG table ----
eggnog <- read_tsv("../../A.Americanus._Transcriptome_Annotations_new.tsv")


annotated_module <- annotate_kegg_column(
  df = eggnog,
  column = "KEGG_Module",
  ref = refs$module,
  ref_key = "module_id"
)


annotated_pathway <- annotate_kegg_column(
  df = eggnog,
  column = "KEGG_Pathway",
  ref = refs$pathway,
  ref_key = "kegg_id"
)


df = eggnog
column = "KEGG_ko"
ref = refs$ko
ref_key = "ko_id"



# Had to mutate the KO: out of the id
col_sym <- rlang::sym(column)

annotated_ko <- df %>%
  tidyr::separate_longer_delim(!!col_sym, delim = ",") %>%
  dplyr::select(!!col_sym) %>%
  dplyr::mutate(!!col_sym := stringr::str_remove(!!col_sym, "^ko:")) %>%
  dplyr::rename(!!ref_key := !!col_sym)%>%
  dplyr::left_join(ref, by = ref_key) %>%
  tidyr::drop_na()



# ---- Step 4: Write outputs ----
write_tsv(annotated_ko, "eggnog_ko_annotated.tsv")
write_tsv(annotated_pathway, "eggnog_pathway_annotated.tsv")
write_tsv(annotated_module, "eggnog_module_annotated.tsv")

# ============================================================
# END
# ============================================================

# Join the Annotations Back with out Genes:

