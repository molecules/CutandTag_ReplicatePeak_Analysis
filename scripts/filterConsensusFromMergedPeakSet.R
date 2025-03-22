#!/usr/bin/env Rscript

# Author: Kevin Boyd
# Purpose: Filter merged peaks to retain only those supported by >= n replicates
# Usage:
# Rscript filterConsensusFromMergedPeakSet.R merged_peaks.bed rep1.bed rep2.bed rep3.bed 2 output_consensus.bed

library(GenomicRanges)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
input_files <- args[1:(length(args) - 2)]
min_overlap <- as.integer(args[length(args) - 1])
output_file <- args[length(args)]

# Load merged peak set (first file)
merged_peaks_file <- input_files[1]
merged_peaks <- import(merged_peaks_file, format = "BED")

# Load replicate peak sets
replicate_files <- input_files[-1]
replicate_grs <- lapply(replicate_files, function(f) import(f, format = "BED"))

# Count how many replicates each merged peak overlaps
overlap_counts <- sapply(replicate_grs, function(gr) countOverlaps(merged_peaks, gr) > 0)
overlap_matrix <- as.data.frame(overlap_counts)
support_counts <- rowSums(overlap_matrix)

# Filter peaks supported by >= min_overlap replicates
consensus_peaks <- merged_peaks[support_counts >= min_overlap]

# Export to BED
export.bed(consensus_peaks, output_file)
