# generateConsensusPeaks.R
# Author: Kevin Boyd (modified from an original by Chris Sansam)
# Date Modified: 12/4/2024
# Purpose: This script calculates overlaps between individual peak files
#          to generate a consensus peaks file.

# Load required libraries
library(magrittr)
library(GenomicRanges)
library(stringr)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
peak_files_input <- args[-length(args)]  # All arguments except the last are input peak files
output <- args[length(args)]            # Last argument is the output file

# Process input peak files
peaks_list <- lapply(peak_files_input, function(file) {
    read.table(file, header = FALSE) %>%
    GenomicRanges::makeGRangesFromDataFrame(
        ignore.strand = TRUE,
        seqnames.field = "V1",
        start.field = "V2",
        end.field = "V3"
    )
})

# Combine all peaks into a single GRanges object
merged_peaks <- reduce(do.call(c, peaks_list))

# Count overlaps across input peak files
overlap_counts <- countOverlaps(merged_peaks, peaks_list)

# Filter peaks based on overlap threshold
consensus_peaks <- merged_peaks[overlap_counts >= 2]  # Adjust threshold as needed (e.g., 2 out of 3)

# Save consensus peaks to RDS file
saveRDS(consensus_peaks, file = output)
