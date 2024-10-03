# Load necessary libraries
library(SeuratDisk)
library(Seurat)
library(argparse)

# Function to parse command line arguments
parse_args <- function() {
  parser <- ArgumentParser(description = 'Process and convert Seurat object')
  
  # Adding arguments
  parser$add_argument('--input', required = TRUE, help = 'Path to the input RDS file')
  parser$add_argument('--output', required = TRUE, help = 'Output file name without extension')
  
  # Parsing arguments
  args <- parser$parse_args()
  return(args)
}

# Main function
main <- function() {
  # Parse arguments
  args <- parse_args()
  
  # Read the Seurat object from the RDS file
  sc <- readRDS(args$input)
  
  # Apply SCTransform
  seurat <- SCTransform(sc, verbose = FALSE)
  
  # Define output file names
  h5seurat_filename <- paste0(args$output, ".h5Seurat")
  h5ad_filename <- paste0(args$output, ".h5ad")
  
  # Save as H5Seurat
  SaveH5Seurat(seurat, filename = h5seurat_filename)
  
  # Convert to h5ad
  Convert(h5seurat_filename, dest = "h5ad")
  
  cat("Processing and conversion completed.\n")
}

# Run the main function
main()
