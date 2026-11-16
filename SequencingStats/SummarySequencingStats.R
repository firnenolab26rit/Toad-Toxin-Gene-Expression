
library(dplyr)
library(stringr)
library(readr)

#### LOGIC FOR WRANGLING HISAT2 LOG #####

# Input file containing concatenated HISAT2 logs
input_file <- "Alignment_Counts_Hisat2Job_86658.out"

# Read entire file
log_lines <- readLines(input_file)

# Identify the beginning of each sample
sample_starts <- grep("^===== .* =====$", log_lines)

# Extract each sample's information
results <- lapply(seq_along(sample_starts), function(i) {
  
  start <- sample_starts[i]
  
  # Determine end of current sample block
  end <- if (i < length(sample_starts)) {
    sample_starts[i + 1] - 1
  } else {
    length(log_lines)
  }
  
  block <- log_lines[start:end]
  
  # Extract sample name
  sample_name <- str_match(
    block[1],
    "^===== (.*) =====$"
  )[, 2]
  
  # Extract total read pairs
  total_pairs <- as.numeric(
    str_match(
      paste(block, collapse = "\n"),
      "(\\d+) reads; of these:"
    )[, 2]
  )
  
  # Find the mates section
  mates_line <- grep("mates make up the pairs", block)
  
  if (length(mates_line) == 0) {
    stop(paste("Could not find mates section for", sample_name))
  }
  
  # Extract number of individual reads that aligned 0 times
  unaligned_reads <- as.numeric(
    str_match(
      paste(block[(mates_line + 1):length(block)],
            collapse = "\n"),
      "(\\d+) \\([0-9.]+%\\) aligned 0 times"
    )[, 2]
  )
  
  # Extract reported HISAT2 overall alignment rate
  alignment_rate <- as.numeric(
    str_match(
      paste(block, collapse = "\n"),
      "([0-9.]+)% overall alignment rate"
    )[, 2]
  )
  
  # Convert read pairs to individual reads
  total_reads <- total_pairs * 2
  
  # Calculate number of individual reads mapped at least once
  mapped_reads <- total_reads - unaligned_reads
  
  # Calculate alignment rate from read counts
  calculated_rate <- (mapped_reads / total_reads) * 100
  
  data.frame(
    Sample = sample_name,
    Read_pairs = total_pairs,
    Total_reads = total_reads,
    Mapped_reads = mapped_reads,
    Unmapped_reads = unaligned_reads,
    HISAT2_alignment_rate = alignment_rate,
    Calculated_alignment_rate = calculated_rate
  )
})

# Combine sample results
results <- bind_rows(results)

# Calculate overall statistics across all samples
overall <- results %>%
  summarise(
    Total_samples = n(),
    Total_read_pairs = sum(Read_pairs),
    Total_reads = sum(Total_reads),
    Total_mapped_reads = sum(Mapped_reads),
    Total_unmapped_reads = sum(Unmapped_reads),
    
    # Weighted overall alignment rate
    Overall_alignment_rate =
      (Total_mapped_reads / Total_reads) * 100,
    
    # Simple mean of individual sample alignment rates
    Mean_sample_alignment_rate =
      mean(HISAT2_alignment_rate)
  )

# Display results
print(results)
print(overall)

print(paste(
  "A total of",
  round(overall$Overall_alignment_rate, 2),
  "% (",
  format(round(overall$Total_mapped_reads), big.mark = ","),
  "reads) mapped to the assembled transcriptome."
))

# Save results
#write_csv(results, "HISAT2_sample_alignment_results.csv")
#write_csv(overall, "HISAT2_overall_alignment_summary.csv")


### QC1 Summary Stats ###
QC <- read_tsv("MG_Toad_QC1_MULTQC/general_stats_table_QC1.tsv")

# Read in `Seqs` in millions of sequnces so * 1,000,000 puts in the full number (sequences)
QC <- QC %>%  mutate(sequences = Seqs * 1e6)
  
t <- sum(QC$sequences)
a <- mean(QC$sequences)
s <- sd(QC$sequences)

### QC2 Summary Stats ###

QC2 <- read_tsv("MG_Toad_MULTIQC/general_stats_table.tsv")

# Read in `Seqs` in millions of sequnces so * 1,000,000 puts in the full number (sequences)
QC2 <- QC2 %>%  mutate(sequences = Seqs * 1e6)

t <- sum(QC2$sequences)
a <- mean(QC2$sequences)
s <- sd(QC2$sequences)

