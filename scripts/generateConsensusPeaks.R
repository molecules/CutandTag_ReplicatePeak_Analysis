# Author: Inistial Chris Sansam script modified by Kevin Boyd
# Date Modified: 12/4/2024
# Purpose: This script calculates the overlap between a merged peaks file 
#          and a list of additional peak files, storing the overlap information 
#          as metadata in the merged peaks object. The merged peaks object 
#          is saved in RDS format as the output file.

# Load necessary packages
library(magrittr)
library(GenomicRanges)
library(stringr)

# Get input and output files and parameter
merged_peaks_file <- snakemake@input[[1]] # Input merged peaks file
output <- snakemake@output[[1]] # Output file
peak_files_input <- snakemake@params[[1]] # Input peak files parameter

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
