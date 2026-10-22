#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --job-name=SalmonQuant
#SBATCH --error=./ERR/%x.%J.err
#SBATCH --output=./OUT/%x.%J.out
#SBATCH --nodes=1
#SBATCH --ntasks=1              
#SBATCH --cpus-per-task=4
#SBATCH --partition=batch-high,batch-medium,batch-low

# Conda ENV
set +eou pipefail
source ~/.bashrc
conda activate "/home/firneno/LEONE_ANALYSIS/conda_env/RNAseq"
set -eou pipefail

#==================================================================#
# Quantification from a de novo assembled transcriptome with Salmon

# AUTHOR: Gianna Leone
# NOVEMBER 17, 2025
#==================================================================#

# This script uses the 'quasi-mapping mode' of Salmon and therefore has two steps:
# 1). Index the reference transcriptome 
# 2). Quantify reads

####################################################################
#     Assign Variables
####################################################################

#Working directory
WKDIR="/home/firneno/LEONE_ANALYSIS"

#Output directory
OUTDIR="${WKDIR}/SALMON_QUANT_reduc"

#Path to annotated transcriptome
TRANSCRIPTOME="${WKDIR}/TRINITY_TRANSCRIPTOME/MultiTissue_RNAseq/DeNovo_Transcriptome/Transcript_Reduction/A_americanus_trinity.okay.mrna.gz"

#Path to trimmed reads
READS="${WKDIR}/TRIM_QC2/02_TRIMMED_DATA"

#Path to Index
INDEX="${WKDIR}/A.Americanus_reduc_index"

####################################################################
#     Commands
####################################################################

#Set working directory
cd "$WKDIR"

#Make Output Directory
mkdir -p "$OUTDIR" 

######################################
#     Create Indicies:
######################################

#IF INDEX DOES NOT ALREADY EXIST, index the annotated transcriptome
if [ ! -d "$INDEX" ]; then
    salmon index -t "$TRANSCRIPTOME" -i "$INDEX" -k 31
else
    echo "Index already exists — skipping indexing."
fi

######################################
#     Quantify:
######################################

cd "$READS"

#Find matching paired-end reads
for R1 in *R1*.fastq.gz; do
    
    #Infer matching R2 file name
    R2="${R1/_R1_/_R2_}"

    #Confirm inferred R2 file exists
    if [ ! -f "$R2" ]; then
        echo "WARNING: R2 file not found for $R1 — skipping this pair."
        continue
    fi

    # Define sample name from filename prefix
    SAMPLE=$(echo "$R1" | sed 's/_R1_.*//')

    echo "==============================="
    echo " Quantifying sample: $SAMPLE"
    echo " R1: $R1"
    echo " R2: $R2"
    echo "==============================="

    #Quanitfy Reads
    salmon quant -i "$INDEX" -l A -1 "$R1" -2 "$R2" --validateMappings -o "${OUTDIR}/${SAMPLE}"
done

echo "All samples processed successfully!"
