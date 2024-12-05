# generateConsensusPeaks.R
# Author: Kevin Boyd (modified from an original by Chris Sansam)
# Date Modified: 12/4/2024
# Purpose: This script generates a consensus peaks file by calculating the overlap
#          between a merged peaks file and a list of individual peak files. The
#          resulting metadata is saved in an RDS file.

library(magrittr)
library(GenomicRanges)
library(stringr)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
merged_peaks_file <- args[-c(length(args) - 1, length(args))]  # All but last two arguments are input files
peak_files_input <- args[length(args) - 1]                        # Second-to-last argument is the BAM directory
output <- args[length(args)]                         # Last argument is the output directory

# Read in merged peaks file and convert to GRanges object
merged_peaks <- merged_peaks_file %>%
  read.table %>%
  GenomicRanges::makeGRangesFromDataFrame(.,
                                          ignore.strand = TRUE,
                                          seqnames.field = "V1",
                                          start.field = "V2",
                                          end.field = "V3")

# Process each peak file in input parameter
overlaps <- peak_files_input %>%
  stringr::str_split(., pattern = ",") %>% # Split input parameter by comma
  .[[1]] %>% # Extract filenames
  lapply(., read.table) %>% # Read in each peak file
  magrittr::set_names(stringr::str_split(peak_files_input, pattern = ",") %>%
                        .[[1]] %>%
                        gsub("^.*/", "", .) %>%
                        gsub("_peaks.narrowPeak", "", .) %>%
                        gsub("_(?=[^_]*$).*", "", ., perl = TRUE)) %>% # Extract sample names from filenames
  lapply(., GenomicRanges::makeGRangesFromDataFrame,
         ignore.strand = TRUE,
         seqnames.field = "V1",
         start.field = "V2",
         end.field = "V3") %>% # Convert each peak file to GRanges object
  lapply(., function(gr) { merged_peaks %over% gr }) %>% # Calculate overlap between each peak file and merged peaks
  do.call(cbind, .) # Convert list of overlaps to matrix

# Add overlaps matrix as metadata to merged peaks object
GenomicRanges::mcols(merged_peaks) <- overlaps

# Save merged peaks object in RDS format
saveRDS(merged_peaks, file = output)
