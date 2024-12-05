# generateMergedPeaks.R
# Author: Kevin Boyd (modified from an original by Chris Sansam)
# Date Modified: 12/5/2024
# Purpose: This script merges MACS2 peak files into a single consensus peak set
#          and calculates the number of overlaps across all input peak files.

library(magrittr)
library(GenomicRanges)
library(stringr)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
peak_files <- args[-length(args)]  # All but last argument are input files
output <- args[length(args)]      # Last argument is the output file

# Load MACS2 peak files as GRanges objects
gr_list <- lapply(peak_files, function(file) {
  peaks <- read.table(file, header = FALSE)
  GenomicRanges::makeGRangesFromDataFrame(
    peaks,
    ignore.strand = TRUE,
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3"
  )
})

# Merge all peaks into a single GRanges object
merged_peaks <- reduce(do.call(c, gr_list))  # Merge all ranges into a single GRanges

# Count overlaps for each merged peak across input files
overlap_counts <- sapply(gr_list, function(gr) countOverlaps(merged_peaks, gr))

# Add overlap counts as metadata
GenomicRanges::mcols(merged_peaks) <- data.frame(overlap_counts)

# Save the merged peaks as an RDS file
saveRDS(merged_peaks, file = output)
