# Define arguments to be passed from command line
args <- commandArgs(trailingOnly = TRUE)

# Check for required arguments
if (length(args) < 3) {
  stop("Error: Please provide data directory, barcode file, and output directory as arguments.")
}

# Extract arguments
data_dir <- args[1]
barcode_file <- args[2]
output_dir <- args[3]

# Load libraries
library(Seurat)
library(Matrix)
library(data.table)

# Read the 10x data assuming Cell Ranger output format
my_data <- Read10X(data.dir = data_dir)

# Read barcode filter file (one barcode per line, no header)
barcode_list <- read.table(barcode_file, header = FALSE, stringsAsFactors = FALSE)$V1

# Create Seurat object
seurat_object <- CreateSeuratObject(counts = my_data)

# Subset data based on barcode filter
filtered_data <- subset(seurat_object, cells = barcode_list)

# Extract expression matrix (features) using the `layer` argument
filtered_matrix <- GetAssayData(filtered_data, layer = "counts")

# Extract cell barcodes
filtered_barcodes <- colnames(filtered_matrix)

# Extract gene names
filtered_features <- rownames(filtered_matrix)

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save the filtered matrix in the Matrix Market format
matrix_file <- file.path(output_dir, "matrix.mtx")
writeMM(filtered_matrix, file = matrix_file)

# Add metadata to the matrix file
metadata <- '%%MatrixMarket matrix coordinate integer general\n%metadata_json: {"software_version": "cellranger-8.0.1", "format_version": 2}'
matrix_data <- readLines(matrix_file)
matrix_data <- c(metadata, matrix_data)
writeLines(matrix_data, matrix_file)
system(glue::glue("gzip -f {matrix_file}"))

# Save barcodes and features in the 10x format as gzipped TSV files
write.table(filtered_barcodes, file = gzfile(file.path(output_dir, "barcodes.tsv.gz")), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Create features.tsv with appropriate format
feature_data <- data.frame(
  V1 = filtered_features,
  V2 = filtered_features, # Assuming the second column is the same as the first, adjust if needed
  V3 = "Gene Expression"
)
write.table(feature_data, file = gzfile(file.path(output_dir, "features.tsv.gz")), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

print("Subsetting and file generation completed!")

