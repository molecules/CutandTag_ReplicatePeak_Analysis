# midpointAndPeakOverlaps.R
# Author: Kevin Boyd
# Purpose: Generate midpoint and full overlap BED files from a list of reproducible peaks.
# Date Modified: 12/5/2024

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_beds <- args[-length(args)]  # All but the last argument are input BED files
output_prefix <- args[length(args)]  # The last argument is the output prefix

# Load BED files
Beds.df <- lapply(input_beds, read.table)

# Create genomic ranges objects and reduce them
Beds.gr <- lapply(Beds.df, GenomicRanges::makeGRangesFromDataFrame, 
                  seqnames.field = "V1", start.field = "V2", end.field = "V3")
AllBeds.gr <- GenomicRanges::reduce(do.call("c", Beds.gr))

# Generate overlap metadata
peakOverlaps <- lapply(Beds.gr, function(gr) { AllBeds.gr %over% gr }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame() %>% 
  set_names(gsub("_0.05.*", "", basename(input_beds))) %>% 
  {mcols(AllBeds.gr) <- .; AllBeds.gr}

# Calculate midpoints for overlaps
midpointOverlaps <- peakOverlaps %>%
  {start(.) <- rowMeans(cbind(start(.), end(.))); . } %>%
  {end(.) <- start(.) + 1; .}

# Export results
export.bed(midpointOverlaps, paste0(output_prefix, "_MidpointOverlaps.bed"))
export.bed(peakOverlaps, paste0(output_prefix, "_PeakOverlaps.bed"))
