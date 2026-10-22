######################################################################################
#        Author: Gianna Leone

#        Goal: Generate counts file from Salmon quant files with tximport
######################################################################################

library(tximport)
library(tidyverse)


######################################################################################
#        ASSIGN VARIABLES
######################################################################################

QUANTS_DIR <- "/home/firneno/LEONE_ANALYSIS/SALMON_QUANT_reduc"
TX2GENE <- "/home/firneno/LEONE_ANALYSIS/GENE_MAP/reduc_gene_map.tsv"
OUTDIR <- "/home/firneno/LEONE_ANALYSIS/TXIMPORT_COUNTS_reduc"


######################################################################################
#        COMMANDS
######################################################################################

#Create the output directory
dir.create(OUTDIR)

#Read in the gene-to-transcript map
tx2gene <- read.delim(TX2GENE, header = FALSE)

#I need the transcript column first and the gene column second
tx2gene <- tx2gene[, c(2, 1)]

#Assign column names to tx2gene
colnames(tx2gene) <- c("transcript", "gene")


# FIND QUANT.SF FILES
#############################

samples <- list.dirs(QUANTS_DIR, full.names = FALSE, recursive = FALSE)
files <- file.path(QUANTS_DIR, samples, "quant.sf")
files <- files[file.exists(files)]

#Make sure that the columns in matrix after running tximport correlate to my sample names
names(files) <- samples


# RUN TXIMPORT
#############################

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

write.csv(txi$counts,    file.path(OUTDIR, "fixed_reduc_gene_counts_matrix.csv"))
write.csv(txi$abundance, file.path(OUTDIR, "fixed_reduc_gene_tpm_matrix.csv"))
write.csv(txi$length,    file.path(OUTDIR, "fixed_reduc_gene_length_matrix.csv"))

message("tximport complete! Output saved to: ", OUTDIR)

