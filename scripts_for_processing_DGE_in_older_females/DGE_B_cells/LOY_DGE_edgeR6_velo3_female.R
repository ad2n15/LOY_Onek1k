library(Seurat)
library(dplyr)
library(ggplot2)
library(limma)
library(ggrepel)
library(tidyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(grid)





# Load the Seurat object
data <- readRDS("/mainfs/ddnb/Ahmed/Data/OneK1K/b3e8792e-a31b-404b-a866-250c43bc06d5.rds")

# Subset male cells only
male_cells <- subset(data, subset = sex == "female")

# Define cell types to include
cell_types_to_include <- c("B naive", "B memory", "B intermediate")

# Filter for these cell types
male_cells <- subset(male_cells, subset = predicted.celltype.l2 %in% cell_types_to_include)

# Step 4: Define T cell markers and CD34 marker
t_cell_markers <- c("ENSG00000010610", "ENSG00000160654", "ENSG00000167286", 
                    "ENSG00000198851", "ENSG00000153563", "ENSG00000172116")
cd34_marker <- "ENSG00000174059"

# Remove cells expressing T cell or CD34 markers
# Check expression of the markers and subset to exclude cells with expression > 0
marker_genes <- c(t_cell_markers, cd34_marker)
marker_expression <- GetAssayData(male_cells, assay = "RNA", slot = "data")[marker_genes, ]

# Identify cells with non-zero expression of any marker
cells_to_exclude <- colnames(marker_expression)[colSums(marker_expression > 0) > 0]

# Subset to keep only cells not expressing these markers
male_cells <- subset(male_cells, cells = setdiff(Cells(male_cells), cells_to_exclude))

# Final filtered dataset
male_cells

umap_plot <- DimPlot(male_cells, reduction = "umap", group.by = "predicted.celltype.l2") +
  ggtitle("UMAP Clustering of Male Cells")
ggsave(filename = "UMAP_Clustering_FeMale_B_Cells_original.png", plot = umap_plot,
       width = 8, height = 6, dpi = 300)



#############################
male_cells@assays$RNA@counts <- male_cells@assays$RNA@data

male_cells <- NormalizeData(male_cells, normalization.method = "LogNormalize", scale.factor = 10000)
male_cells <- FindVariableFeatures(male_cells, selection.method = "vst", nfeatures = 500)
male_cells <- ScaleData(male_cells , vars.to.regress = c('percent.mt', 'pool_number'))

# Run PCA on the SCT scaled data (Pearson residuals)
male_cells <- RunPCA(male_cells, verbose = TRUE)

# Use the top 30 principal components
num_pcs <- 30

# Perform nearest-neighbor graph-based clustering based on the PCA result
male_cells <- FindNeighbors(male_cells, dims = 1:num_pcs)
male_cells <- FindClusters(male_cells, resolution = 0.5)

# Check the clustering results
print(male_cells)

# Optionally, plot PCA variance ratio to inspect the contribution of each PC
ElbowPlot(male_cells, ndims = 50)

# Run UMAP based on the PCA results and using the SCT data
male_cells <- RunUMAP(male_cells, dims = 1:num_pcs)

# save
# Save the Seurat object
saveRDS(male_cells, file = "female_B_cells_191124.rds")



# Visualize UMAP with clusters
umap_plot <- DimPlot(male_cells, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("UMAP Clustering of Male Cells")
ggsave(filename = "UMAP_Clustering_FeMale_B_Cells_instructions.png", plot = umap_plot,
       width = 8, height = 6, dpi = 300)

# Check the cluster assignments
table(male_cells$seurat_clusters)


# Visualize UMAP with predicted_celltypes
umap_plot <- DimPlot(male_cells, reduction = "umap", group.by = "predicted.celltype.l2") +
  ggtitle("UMAP Clustering of Male Cells")
ggsave(filename = "UMAP_Clustering_FeMale_B_Cells_azimuth_instructions.png", plot = umap_plot,
       width = 8, height = 6, dpi = 300)

##########################################################



# Find all markers
cluster_markers <- FindAllMarkers(
  object = male_cells,
  min.pct = 0.1,
  logfc.threshold = 0.1
)

#######################################################

# Gene dictionary with Ensembl IDs and symbols
gene_symbols <- list(
  'ENSG00000077238' = 'IL4R',
  'ENSG00000100721' = 'TCL1A',
  'ENSG00000104921' = 'FCER2',
  'ENSG00000139193' = 'CD27',
  'ENSG00000156738' = 'MS4A1',
  'ENSG00000240505' = 'TNFRSF13B'
)

# Get Ensembl IDs of interest
genes_of_interest <- names(gene_symbols)

# Filter cluster_markers for genes of interest
filtered_markers <- cluster_markers[cluster_markers$gene %in% genes_of_interest, ]

# Replace Ensembl IDs with gene symbols in the 'gene' column
filtered_markers$gene <- sapply(filtered_markers$gene, function(x) {
  if(x %in% names(gene_symbols)) {
    return(gene_symbols[[x]])  # Return gene symbol if found
  } else {
    return(x)  # Return Ensembl ID if not found in the dictionary
  }
})

# Rename the column from 'gene' to 'gene_name'
colnames(filtered_markers)[colnames(filtered_markers) == "gene"] <- "gene_name"

# Print the filtered markers with gene names
print(filtered_markers)

# Optionally, format the output to show only certain columns, such as p-value, logFC, and adjusted p-value
filtered_markers %>%
  select(gene_name, p_val, avg_log2FC, p_val_adj) %>%
  print()

# Save filtered markers to a CSV file
write.csv(filtered_markers, "filtered_Femarkers_B_cells_clusters_191124_female.csv", row.names = FALSE)

######################################################

