# midpointAndPeakOverlaps.R
# Author: Kevin Boyd
# Purpose: Generate midpoint and full overlap BED files from a list of reproducible peaks,
#          and additionally compute unique sample midpoints (unique to only one sample).
# Date Modified: 12/5/2024

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# All but the last argument are input BED files; the last argument is the output directory for overlaps.
input_beds <- args[-length(args)]
output_dir <- args[length(args)]

# Ensure the output directory ends with a slash
if (!grepl("/$", output_dir)) {
    output_dir <- paste0(output_dir, "/")
}

# -----------------------------
# Original Overlap and Midpoint Calculation
# -----------------------------

# Load BED files
Beds.df <- lapply(input_beds, read.table, header = FALSE)

# Create GRanges objects from each BED file
Beds.gr <- lapply(Beds.df, makeGRangesFromDataFrame,
                  seqnames.field = "V1", start.field = "V2", end.field = "V3")

# Merge all GRanges and reduce them to form a unified set of ranges
AllBeds.gr <- reduce(do.call("c", Beds.gr))

# Generate overlap metadata:
# For each input GRanges, check which ranges in AllBeds.gr overlap it.
# Combine the resulting logical vectors into a data frame and attach as metadata.
peakOverlaps <- lapply(Beds.gr, function(gr) { AllBeds.gr %over% gr }) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  # Use the file basenames (with "_consensus_peaks.bed" removed) as column names
  set_names(gsub("_consensus_peaks.bed", "", basename(input_beds))) %>%
  { mcols(AllBeds.gr) <- .; AllBeds.gr }

# Calculate midpoints for the consensus overlaps.
midpointOverlaps <- peakOverlaps %>%
  { start(.) <- rowMeans(cbind(start(.), end(.))); . } %>%
  { end(.) <- start(.) + 1; . }

# Export the overall midpoint and full overlaps files to the specified output directory.
export.bed(midpointOverlaps, paste0(output_dir, "MidpointOverlaps.bed"))
export.bed(peakOverlaps, paste0(output_dir, "PeakOverlaps.bed"))

# -----------------------------
# Additional: Calculate Unique Sample Midpoints
# -----------------------------
# For each sample (i.e. each column in the overlap metadata), extract the consensus peaks
# that are contributed by only that sample (where the row sum is 1 and that sample is TRUE).
# Then recalculate the midpoint for those peaks and export them as a BED file in the consensus peaks folder

# Assume the consensus peaks folder is the same as the folder of the input BED files.
consensus_dir <- dirname(input_beds[[1]])

# Convert the metadata (overlap matrix) to a data frame for easy subsetting.
overlap_df <- as.data.frame(mcols(peakOverlaps))

# Loop over each sample column.
for(sample in colnames(overlap_df)) {
    # Find indices where this sample is TRUE and the total number of TRUE values in that row is 1.
    unique_idx <- which(overlap_df[[sample]] & (rowSums(overlap_df) == 1))
    
    # Subset the consensus GRanges to obtain peaks unique to this sample.
    unique_peaks <- peakOverlaps[unique_idx]
    
    # If there are any unique peaks, recalc their midpoint.
    if (length(unique_peaks) > 0) {
        new_midpoints <- round((start(unique_peaks) + end(unique_peaks)) / 2)
        start(unique_peaks) <- new_midpoints
        end(unique_peaks) <- new_midpoints + 1
    }
    
    # Define the output filename. The file will be placed in the consensus peaks folder.
    out_file <- file.path(consensus_dir, paste0(sample, "_unique_MP.bed"))
    
    # Export the unique midpoints BED file.
    export.bed(unique_peaks, out_file)
}
