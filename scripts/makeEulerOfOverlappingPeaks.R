###############################################
# Script name: makeEulerOfOverlappingPeaks.R   
# Author: Kevin Boyd (modified from Chris Sansam)
# Date Modified: 12/5/2024                        
# Purpose: Generate an Euler plot of overlapping peaks that have been merged.
#          Input is an RDS file containing the merged peaks with scores 
#          indicating overlap with each of the replicate peak sets.
#          Output is an RDS and PDF file with a plot. Parameters are provided 
#          via command-line arguments.                          
###############################################

# Load required packages.
library(magrittr)
library(eulerr)
library(GenomicRanges)
library(ggplotify)
library(stringr)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]        # Input RDS file
output_rds <- args[2]        # Output RDS file for the plot
output_pdf <- args[3]        # Output PDF file for the plot
font_size <- as.numeric(args[4])  # Font size for the plot
colors <- str_split(args[5], pattern = ",")[[1]]  # Colors for the plot
pdf_width <- as.numeric(args[6])  # Width of the PDF
pdf_height <- as.numeric(args[7])  # Height of the PDF

# Load the merged peaks input file
merged_peaks <- readRDS(input_file)

# Set options for the Euler plot
eulerr_options(labels = list(fontsize = font_size),
               quantities = list(fontsize = font_size - 2, padding = grid::unit(100, "mm")),
               legend = list(fontsize = font_size, vgap = 0.01))

# Create the Euler plot
EulerPlot <- GenomicRanges::mcols(merged_peaks) %>%
    as.matrix() %>%
    euler(., shape = "ellipse") %>%
    plot(., quantities = TRUE, legend = TRUE, adjust_labels = TRUE, fills = colors) %>%
    as.ggplot()

# Save the Euler plot as an RDS file
saveRDS(EulerPlot, output_rds)

# Save the Euler plot as a PDF file
pdf(output_pdf, width = pdf_width, height = pdf_height)
print(EulerPlot)
dev.off()
