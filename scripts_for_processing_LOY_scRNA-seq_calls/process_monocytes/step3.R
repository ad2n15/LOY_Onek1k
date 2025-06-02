library(phateR)
library(reticulate)
library(Seurat)
library(dplyr)
library(ggplot2)


male_cells <- readRDS("male_mono_cells_191124_2.rds")
variable_features <- VariableFeatures(male_cells)


# Subset the normalized expression matrix to include only variable features
data_matrix <- as.matrix(GetAssayData(male_cells, slot = "scale.data")[variable_features, ])

# Transpose the matrix
transposed_data_matrix <- t(data_matrix)

# Save the transposed matrix as a CSV file
write.csv(transposed_data_matrix, "data_matrix_mono_cells_for_phate.csv", row.names = TRUE)

