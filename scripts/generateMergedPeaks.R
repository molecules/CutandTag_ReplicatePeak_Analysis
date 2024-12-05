# generateMergedPeaks.R
# Author: Kevin Boyd (modified from an original by Chris Sansam)
# Date Modified: 12/5/2024
# Purpose: Generate a merged peaks file by combining individual peak files.
#          Overlap information is calculated for each replicate, and the merged
#          peaks are saved as an RDS file with metadata indicating overlaps.

library(magrittr)
library(GenomicRanges)
library(purrr)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[-length(args)]  # All but last argument are input files
output <- args[length(args)]        # Last argument is the output file

# Create a vector of sample names from input file paths
sample_names <- input_files %>%
  basename() %>%
  gsub("_peaks.narrowPeak", "", .)

# Load individual peak files into GRanges objects and combine them
peak_list <- input_files %>%
  map(~ read.table(.x, header = FALSE) %>%
        GenomicRanges::makeGRangesFromDataFrame(
          ignore.strand = TRUE,
          seqnames.field = "V1",
          start.field = "V2",
          end.field = "V3"
        ))

# Merge peaks across all files
merged_peaks <- Reduce(function(x, y) GenomicRanges::union(x, y), peak_list)

# Calculate overlaps for each file and create an overlap matrix
overlap_matrix <- map_dfc(peak_list, ~ GenomicRanges::countOverlaps(merged_peaks, .x))

# Assign meaningful column names based on sample names
colnames(overlap_matrix) <- sample_names

# Add overlap matrix as metadata to the merged peaks object
GenomicRanges::mcols(merged_peaks) <- overlap_matrix

# Save the merged peaks with metadata to an RDS file
saveRDS(merged_peaks, file = output)
