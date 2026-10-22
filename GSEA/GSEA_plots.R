library(readr)
library(clusterProfiler)

# Load the data
GSEA_obj <- read_rds("results/Tissue/Deseq_Modules_Input_GSEA.rds")
# Make a df for reading 
GSEA_df <- GSEA_obj@result

# Write a minimal table to a csv for Sup Table S9
GSEA_df %>% select(ID, Description, setSize, NES, p.adjust) %>% write_csv(file = "TableS9.GSEA_Results.csv" )

# Global options 
options(enrichplot.colours = c("orange","grey"))


# Keep 
vec <- c(
  "Glycine, serine and threonine metabolism",
  "Estrogen signaling pathway",
  "Metabolism of xenobiotics by cytochrome P450",
  "Retinol metabolism",
  "Cholesterol metabolism",
  "Steroid hormone biosynthesis",
  "Drug metabolism",
  "Microbial metabolism in diverse environments",
  "Ascorbate and aldarate metabolism",
  "Drug metabolism",
  "Biosynthesis of secondary metabolites",
  "Protein export",
  "Histidine metabolism",
  "Fat digestion and absorption",
  "Arginine biosynthesis",
  "Carbon metabolism",
  "Glycosphingolipid biosynthesis",
  "Biosynthesis of amino acids",
  "beta-Alanine metabolism",
  "One carbon pool by folate",
  "Glyoxylate and dicarboxylate metabolism"
)


p <- dotplot(GSEA_obj, showCategory=vec, label_format=NULL)

ggsave(filename = "./results/Tissue/GSEA_Tissue_KeggPathwaysOfInterest.pdf", plot = p, dpi = 300, width = 8.27, height = 11.69, units = "in")
write_rds(p, file = "./results/Tissue/GSEA_Tissue_KeggPathwaysOfInterest_plot.rds")


cnetplot(GSEA_obj)


p1 <- cnetplot(GSEA_obj, node_label="category") 
p2 <- cnetplot(GSEA_obj, node_label="gene") 
p3 <- cnetplot(GSEA_obj, node_label="all") 
p4 <- cnetplot(GSEA_obj, node_label="none", 
               color_category='firebrick', 
               color_item='steelblue') 
plot_list(p1, p2, p3, p4, ncol=2, tag_levels = 'A')

heatplot(GSEA_obj)


library(ggtree)

GSEA_obj <- enrichplot::pairwise_termsim(GSEA_obj)
p1 <- enrichplot::treeplot(GSEA_obj, cladelab_offset=8, tiplab_offset=.3, fontsize_cladelab =5) + 
  hexpand(.2)
p2 <- enrichplot::treeplot(GSEA_obj, cluster_method = "average", 
               cladelab_offset=14, tiplab_offset=.3, fontsize_cladelab =5) + 
  hexpand(.3)
aplot::plot_list(p1, p2, tag_levels='A', ncol=2)


ridgeplot(GSEA_obj, stat = bin)
