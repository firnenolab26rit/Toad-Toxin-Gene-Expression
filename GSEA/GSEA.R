# =========================================================================================== #
# GSEA ANALYSIS SCRIPT (clusterProfiler::GSEA)
#
# DESCRIPTION:
# This script performs Gene Set Enrichment Analysis (GSEA) using a ranked gene list
# and custom TERM2GENE / TERM2NAME mappings.
#
# ------------------------------------------------------------------------------------------- #
# INPUT REQUIREMENTS:
#
# 1. Input table (--input):
#    - Must contain:
#        • A gene identifier column (default = first column)
#        • A numeric ranking column (specified with --rank_col)
#    - Example format:
#        gene_id    log2FoldChange
#        geneA      2.3
#        geneB     -1.5
#
# 2. TERM2GENE (--term2gene):
#    - Two-column table:
#        term_id    gene_id
#    - Example:
#        GO:0008150    geneA
#        GO:0008150    geneB
#
# 3. TERM2NAME (--term2name):
#    - Two-column table:
#        term_id    term_description
#
# ------------------------------------------------------------------------------------------- #
# IMPORTANT ASSUMPTIONS:
#
# • Gene column:
#   - By default, the FIRST column of the input file is used as gene IDs.
#   - If your gene column has a specific name (e.g., "trinity_id"), modify:
#       gene_col <- colnames(df)[1]
#
# • Ranking direction:
#   - Gene list is sorted in decreasing order:
#       sort(decreasing = TRUE)
#   - Assumes higher values = stronger upregulation.
#
# • Duplicate genes:
#   - Duplicate gene IDs are removed using distinct().
#   - This keeps the first occurrence only.
#   - Alternative strategies (recommended for RNA-seq):
#       group_by(gene) %>% summarise(rank = max(rank))
#
# ------------------------------------------------------------------------------------------- #
# RECOMMENDED RANKING METRICS (RNA-seq):
#
# Choose a biologically meaningful ranking statistic:
#
# • log2FoldChange
#     - Simple and interpretable
#
# • stat (from DESeq2)
#     - Often preferred (accounts for variance)
#
# • Combined metric (recommended for stronger signal):
#     sign(log2FoldChange) * -log10(pvalue)
#
# ------------------------------------------------------------------------------------------- #
# NOTES ON UNIVERSE:
#
# • The "universe" argument is typically NOT used in GSEA.
# • GSEA operates on ranked lists and does not require a background gene set.
#
# ------------------------------------------------------------------------------------------- #
# OUTPUT:
#
# • CSV file: enrichment results table
# • RDS file: full GSEA object (for downstream plotting/analysis)
#
# ------------------------------------------------------------------------------------------- #
# EXTENSIONS (future improvements):
#
# • Batch GSEA across multiple contrasts
# • Automatic DESeq2 input parsing
# • Visualization (dotplot, ridgeplot, gseaplot2)
# • Integration with WGCNA modules
#
# =========================================================================================== #


library(clusterProfiler)
library(optparse)
library(readr)
library(dplyr)

# ----------------------------
# CLI options
# ----------------------------
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              help = "Input table with gene IDs and ranking column"),
  
  make_option(c("-g", "--term2gene"),
              type = "character",
              help = "Term2Gene TSV file"),
  
  make_option(c("-n", "--term2name"),
              type = "character",
              help = "Term2Name TSV file"),
  
  make_option(c("-r", "--rank_col"),
              type = "character",
              help = "Column name for ranking (e.g., log2FoldChange)"),
  
  make_option(c("-o", "--output"),
              type = "character",
              default = "./",
              help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ----------------------------
# Validate inputs
# ----------------------------
if (is.null(opt$input) || is.null(opt$term2gene) || 
    is.null(opt$term2name) || is.null(opt$rank_col)) {
  stop("Error: --input, --term2gene, --term2name, and --rank_col are required.", call. = FALSE)
}

# ----------------------------
# Load data
# ----------------------------

# Input table (must contain gene IDs + ranking column)
df <- read_csv(opt$input, show_col_types = FALSE)

# Check ranking column exists
if (!(opt$rank_col %in% colnames(df))) {
  stop(paste("Ranking column not found:", opt$rank_col))
}

# Assume first column is gene ID (adjust if needed)
gene_col <- colnames(df)[1]

# Create ranked gene list
gene_list <- df %>%
  dplyr::select(all_of(gene_col), all_of(opt$rank_col)) %>%
  dplyr::filter(!is.na(.data[[opt$rank_col]])) %>%
  distinct() %>%
  tibble::deframe()

# Sort decreasing (required for GSEA)
gene_list <- sort(gene_list, decreasing = TRUE)

# Load TERM2GENE / TERM2NAME
t2g <- read_tsv(opt$term2gene, col_names = FALSE)
t2n <- read_tsv(opt$term2name, col_names = FALSE)

# ----------------------------
# Run GSEA
# ----------------------------
gsea_results <- GSEA(
  geneList = gene_list,
  TERM2GENE = t2g,
  TERM2NAME = t2n,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# ----------------------------
# Output setup
# ----------------------------
outdir <- opt$output
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Might need to give it a better prefix name
prefix <- tools::file_path_sans_ext(basename(opt$input))

out_csv <- file.path(outdir, paste0(prefix, "_GSEA_results.csv"))
out_rds <- file.path(outdir, paste0(prefix, "_GSEA.rds"))

# ----------------------------
# Save results
# ----------------------------
saveRDS(gsea_results, file = out_rds)

if (!is.null(gsea_results) && nrow(gsea_results@result) > 0) {
  write.csv(gsea_results@result, file = out_csv, row.names = FALSE)
}

# ----------------------------
# Logging
# ----------------------------
cat("#========================================================#\n")
cat("GSEA Results saved to:", out_csv, "\n")
cat("Total genes ranked:", length(gene_list), "\n")

if (!is.null(gsea_results) && nrow(gsea_results@result) > 0) {
  sig_hits <- gsea_results@result %>%
    dplyr::filter(p.adjust <= 0.05)
  
  if (nrow(sig_hits) == 0) {
    cat("No significant enrichment found.\n")
  } else {
    cat("Number of significant pathways:", nrow(sig_hits), "\n")
  }
} else {
  cat("GSEA returned NULL or empty results.\n")
}

cat("#========================================================#\n")
