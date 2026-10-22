# Functions 


# ==============================================================================
# Annotation Joining Utility
# ==============================================================================

annotate_results <- function(result_df, t2g, t2n) {
  # IDENTIFY COLUMN NAMES DYNAMICALLY
  # 1. Find the Gene ID column in the results (e.g., 'trinity_id', 'gene_id', or 'Gene')
  # 2. Find the Term ID column in t2g (the first column)
  # 3. Find the Term Name column in t2n (the second column)
  
  res_gene_col <- colnames(result_df)[1]  # Assumes 1st column is Gene ID
  t2g_term_col <- colnames(t2g)[1]       # e.g., 'GOs' or 'kegg_id'
  t2g_gene_col <- colnames(t2g)[2]       # e.g., 'trinity_id'
  t2n_name_col <- colnames(t2n)[2]       # e.g., 'Name' or 'description'
  
  message(paste0("Joining ", res_gene_col, " with ", t2g_gene_col, "..."))
  
  # EXECUTE JOINS
  annotated_df <- result_df %>%
    # Join results to term2gene map
    # Uses setNames to handle: result_df$ID == t2g$trinity_id
    dplyr::left_join(t2g, by = setNames(t2g_gene_col, res_gene_col)) %>%
    
    # Join result to term2name dictionary
    # Uses setNames to handle: t2g$TermID == t2n$TermID
    dplyr::left_join(t2n, by = t2g_term_col) %>%
    
    # REORGANIZE COLUMNS
    # Move the annotation info to the front, followed by original result data
    dplyr::relocate(all_of(c(res_gene_col, t2g_term_col, t2n_name_col)))
  
  return(annotated_df)
}

# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

# Example with your GO BP tables:
# final_go_results <- annotate_results(
#   result_df = res_degs, # Your DEGs or WGCNA membership table
#   t2g = t2g_go_bp,      # From the previous script
#   t2n = t2n_go_bp       # From the previous script
# )