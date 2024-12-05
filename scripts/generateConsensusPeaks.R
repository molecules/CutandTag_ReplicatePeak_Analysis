# filterConsensusPeaks.R
# Author: Kevin Boyd (modified from an original by Chris Sansam)
# Date Modified: 12/4/2024
# Purpose: This script filters consensus peaks based on the minimum overlap count.
#          The input is an RDS file containing merged peaks with overlap metadata,
#          and the output is a BED file containing the filtered consensus peaks.

library(magrittr)
library(GenomicRanges)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input <- args[-c(length(args) - 1, length(args))]    # All but last two arguments are input files
min <- args[length(args) - 1]                        # Second-to-last argument is the BAM directory
output <- args[length(args)]                         # Last argument is the output directory

# Read in the RDS file with metadata
gr <- readRDS(input)

# Filter peaks based on overlap counts
gr2 <- gr %>%
  .[(gr %>%
      mcols %>%
      as.matrix %>%
      rowSums) >= min]

# Convert filtered GRanges object to BED format
df <- data.frame(seqnames = seqnames(gr2),
                 starts = start(gr2) - 1,
                 ends = end(gr2),
                 names = c(rep(".", length(gr2))),
                 scores = c(rep(".", length(gr2))),
                 strands = strand(gr2))

# Write output to BED file
write.table(df, file = output, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
