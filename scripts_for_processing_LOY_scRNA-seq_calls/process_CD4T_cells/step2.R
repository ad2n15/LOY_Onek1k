library(Seurat)
library(dplyr)
library(ggplot2)
#library(limma)
#library(ggrepel)
#library(tidyr)
#library(pheatmap)
#library(ComplexHeatmap)
#library(circlize)
#library(grid)
#library(MAST)
library(ggrepel)




# Load the Seurat object
male_cells <- readRDS("male_CD4T_cells_191124.rds")


# Check existing clusters
table(male_cells$seurat_clusters)

# Ensure clusters are character
male_cells$seurat_clusters <- as.character(male_cells$seurat_clusters)
male_cells <- subset(x = male_cells, idents = c("1", "5", "6", "12", "14", "16", "17", "18", "19", "20", "21"), invert = TRUE)
#male_cells <- subset(x = male_cells, idents = c("1", "5", "6", "12", "14", "16", "17", "18", "19", "20", "21", "8"), invert = TRUE)

# Classify clusters into "B-in" and "B-mem"
male_cells$new_classification <- male_cells$seurat_clusters
male_cells$new_classification[male_cells$seurat_clusters %in% c("2", "4", "7", "10", "11", "13", "15")] <- "CD4nc"
male_cells$new_classification[male_cells$seurat_clusters %in% c("0", "3")] <- "CD4et"
male_cells$new_classification[male_cells$seurat_clusters %in% c("8")] <- "Treg"
male_cells$new_classification[male_cells$seurat_clusters %in% c("9")] <- "TCL"



# Convert classification to factor
male_cells$new_classification <- factor(male_cells$new_classification)
table(male_cells$new_classification)

# Remove unwanted clusters

#Idents(object = male_cells) <- "seurat_clusters"
#male_cells <- subset(x = male_cells, idents = c("7", "9", "10", "11", "12"), invert = TRUE)
table(male_cells$seurat_clusters)
table(male_cells$new_classification)



saveRDS(male_cells, file = "male_CD4T_cells_191124_2.rds")


# Visualize UMAP with predicted_celltypes
umap_plot <- DimPlot(male_cells, reduction = "umap", group.by = "new_classification") +
  ggtitle("UMAP Clustering of Male Cells")
ggsave(filename = "UMAP_Clustering_Male_CD4T_Cells_newclassification.png", plot = umap_plot,
       width = 8, height = 6, dpi = 300)

##################################################

# Ensure 'B-in' cells are the focus
Idents(object= male_cells) <- "new_classification"  # Set cluster classification as identity
male_cells_B_in <- subset(male_cells, idents = "CD4nc")  # Subset only 'B-in' cells
male_cells_B_in

# Set identity to 'loy_status' for differential analysis
male_cells$LOY_status <- as.character(male_cells$LOY_status)
Idents(male_cells_B_in) <- "LOY_status"


# Step 1: Load the gene list
gene_list <- read.table("longest_transcripts_unique.txt", header = TRUE, sep = "\t")

# Step 2: Filter out chrY genes
non_chrY_genes <- gene_list[gene_list$chromosome != "chrY", "gene"]

# Step 3: Subset the Seurat object
# Assuming 'male_cells_B_in' is your Seurat object
filtered_seurat_object <- subset(
  male_cells_B_in, 
  features = non_chrY_genes
)


# Load the gene mapping file
gene_mapping <- read.csv("Onek1k_mart_export.txt", header = FALSE, sep = ",", col.names = c("gene_id", "gene_name"))

# Perform DGE between LOY (1) and non-LOY (0) within 'B-in' cells
dge_results <- FindMarkers(
  object = filtered_seurat_object,
  ident.1 = "1",  # LOY cells
  ident.2 = "0",  # Non-LOY cells
  test.use = "MAST",  # Default test
  verbose = FALSE
)

# Ensure rownames (gene IDs) are a column for mapping
dge_results <- cbind(gene_id = rownames(dge_results), dge_results)

# Merge DGE results with gene names
dge_results <- merge(dge_results, gene_mapping, by = "gene_id", all.x = TRUE)

# View top DGE results
head(dge_results)

# Save the results to a CSV file
write.csv(dge_results, file = "DGE_results_CD4nc_in_LOY_with_gene_names.csv", row.names = FALSE)

# Required libraries
library(ggplot2)
library(ggrepel)

# Ensure gene_name is available in the DGE results
dge_results$Significance <- ifelse(
  dge_results$p_val_adj < 0.05 & abs(dge_results$avg_log2FC) > 0.1,
  "Significant",
  "Not Significant"
)

# Create the volcano plot with gene names
volcano_plot <- ggplot(dge_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significance), alpha = 0.6) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  labs(
    title = "Volcano Plot: B-in Cells (LOY vs Non-LOY)",
    x = "Average Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_classic() +  # Use a white background theme
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

# Highlight top genes (optional, based on significance and p-value ranking)
top_genes <- head(dge_results[order(dge_results$p_val_adj), ], 10)  # Top 10 by p-value
volcano_plot <- volcano_plot +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name),  # Use the gene_name column for labels
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = 10
  )

# Display the plot
print(volcano_plot)

# Save the volcano plot
ggsave("volcano_plot_CD4nc_LOY.png", volcano_plot, width = 8, height = 6)
###################################################
# Ensure 'B-in' cells are the focus
Idents(object= male_cells) <- "new_classification"  # Set cluster classification as identity
male_cells_B_mem <- subset(male_cells, idents = "CD4et")  # Subset only 'B-in' cells
male_cells_B_mem

# Set identity to 'loy_status' for differential analysis
male_cells$LOY_status <- as.character(male_cells$LOY_status)
Idents(male_cells_B_mem) <- "LOY_status"


# Step 1: Load the gene list
gene_list <- read.table("longest_transcripts_unique.txt", header = TRUE, sep = "\t")

# Step 2: Filter out chrY genes
non_chrY_genes <- gene_list[gene_list$chromosome != "chrY", "gene"]

# Step 3: Subset the Seurat object
# Assuming 'male_cells_B_in' is your Seurat object
filtered_seurat_object <- subset(
  male_cells_B_mem,
  features = non_chrY_genes
)


# Load the gene mapping file
gene_mapping <- read.csv("Onek1k_mart_export.txt", header = FALSE, sep = ",", col.names = c("gene_id", "gene_name"))

# Perform DGE between LOY (1) and non-LOY (0) within 'B-in' cells
dge_results <- FindMarkers(
  object = filtered_seurat_object,
  ident.1 = "1",  # LOY cells
  ident.2 = "0",  # Non-LOY cells
  test.use = "MAST",  # Default test
  verbose = FALSE
)

# Ensure rownames (gene IDs) are a column for mapping
dge_results <- cbind(gene_id = rownames(dge_results), dge_results)

# Merge DGE results with gene names
dge_results <- merge(dge_results, gene_mapping, by = "gene_id", all.x = TRUE)

# View top DGE results
head(dge_results)

# Save the results to a CSV file
write.csv(dge_results, file = "DGE_results_CD4et_mem_LOY_with_gene_names.csv", row.names = FALSE)

# Required libraries
library(ggplot2)
library(ggrepel)

# Ensure gene_name is available in the DGE results
dge_results$Significance <- ifelse(
  dge_results$p_val_adj < 0.05 & abs(dge_results$avg_log2FC) > 0.1,
  "Significant",
  "Not Significant"
)

# Create the volcano plot with gene names
volcano_plot <- ggplot(dge_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significance), alpha = 0.6) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  labs(
    title = "Volcano Plot: B-in Cells (LOY vs Non-LOY)",
    x = "Average Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_classic() +  # Use a white background theme
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

# Highlight top genes (optional, based on significance and p-value ranking)
top_genes <- head(dge_results[order(dge_results$p_val_adj), ], 10)  # Top 10 by p-value
volcano_plot <- volcano_plot +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name),  # Use the gene_name column for labels
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = 10
  )

# Display the plot
print(volcano_plot)

# Save the volcano plot
ggsave("volcano_plot_CD4et_LOY.png", volcano_plot, width = 8, height = 6)


#################


# Ensure 'B-in' cells are the focus
Idents(object= male_cells) <- "new_classification"  # Set cluster classification as identity
male_cells_B_mem <- subset(male_cells, idents = "TCL")  # Subset only 'B-in' cells
male_cells_B_mem

# Set identity to 'loy_status' for differential analysis
male_cells$LOY_status <- as.character(male_cells$LOY_status)
Idents(male_cells_B_mem) <- "LOY_status"


# Step 1: Load the gene list
gene_list <- read.table("longest_transcripts_unique.txt", header = TRUE, sep = "\t")

# Step 2: Filter out chrY genes
non_chrY_genes <- gene_list[gene_list$chromosome != "chrY", "gene"]

# Step 3: Subset the Seurat object
# Assuming 'male_cells_B_in' is your Seurat object
filtered_seurat_object <- subset(
  male_cells_B_mem,
  features = non_chrY_genes
)


# Load the gene mapping file
gene_mapping <- read.csv("Onek1k_mart_export.txt", header = FALSE, sep = ",", col.names = c("gene_id", "gene_name"))

# Perform DGE between LOY (1) and non-LOY (0) within 'B-in' cells
dge_results <- FindMarkers(
  object = filtered_seurat_object,
  ident.1 = "1",  # LOY cells
  ident.2 = "0",  # Non-LOY cells
  test.use = "MAST",  # Default test
  verbose = FALSE
)

# Ensure rownames (gene IDs) are a column for mapping
dge_results <- cbind(gene_id = rownames(dge_results), dge_results)

# Merge DGE results with gene names
dge_results <- merge(dge_results, gene_mapping, by = "gene_id", all.x = TRUE)

# View top DGE results
head(dge_results)

# Save the results to a CSV file
write.csv(dge_results, file = "DGE_results_TCL_LOY_with_gene_names.csv", row.names = FALSE)

# Required libraries
library(ggplot2)
library(ggrepel)

# Ensure gene_name is available in the DGE results
dge_results$Significance <- ifelse(
  dge_results$p_val_adj < 0.05 & abs(dge_results$avg_log2FC) > 0.1,
  "Significant",
  "Not Significant"
)

# Create the volcano plot with gene names
volcano_plot <- ggplot(dge_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significance), alpha = 0.6) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  labs(
    title = "Volcano Plot: B-in Cells (LOY vs Non-LOY)",
    x = "Average Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_classic() +  # Use a white background theme
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

# Highlight top genes (optional, based on significance and p-value ranking)
top_genes <- head(dge_results[order(dge_results$p_val_adj), ], 10)  # Top 10 by p-value
volcano_plot <- volcano_plot +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name),  # Use the gene_name column for labels
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = 10
  )

# Display the plot
print(volcano_plot)

# Save the volcano plot
ggsave("volcano_plot_TCL_LOY.png", volcano_plot, width = 8, height = 6)


#################


# Ensure 'B-in' cells are the focus
Idents(object= male_cells) <- "new_classification"  # Set cluster classification as identity
male_cells_B_mem <- subset(male_cells, idents = "Treg")  # Subset only 'B-in' cells
male_cells_B_mem

# Set identity to 'loy_status' for differential analysis
male_cells$LOY_status <- as.character(male_cells$LOY_status)
Idents(male_cells_B_mem) <- "LOY_status"


# Step 1: Load the gene list
gene_list <- read.table("longest_transcripts_unique.txt", header = TRUE, sep = "\t")

# Step 2: Filter out chrY genes
non_chrY_genes <- gene_list[gene_list$chromosome != "chrY", "gene"]

# Step 3: Subset the Seurat object
# Assuming 'male_cells_B_in' is your Seurat object
filtered_seurat_object <- subset(
  male_cells_B_mem,
  features = non_chrY_genes
)


# Load the gene mapping file
gene_mapping <- read.csv("Onek1k_mart_export.txt", header = FALSE, sep = ",", col.names = c("gene_id", "gene_name"))

# Perform DGE between LOY (1) and non-LOY (0) within 'B-in' cells
dge_results <- FindMarkers(
  object = filtered_seurat_object,
  ident.1 = "1",  # LOY cells
  ident.2 = "0",  # Non-LOY cells
  test.use = "MAST",  # Default test
  verbose = FALSE
)

# Ensure rownames (gene IDs) are a column for mapping
dge_results <- cbind(gene_id = rownames(dge_results), dge_results)

# Merge DGE results with gene names
dge_results <- merge(dge_results, gene_mapping, by = "gene_id", all.x = TRUE)

# View top DGE results
head(dge_results)

# Save the results to a CSV file
write.csv(dge_results, file = "DGE_results_Treg_LOY_with_gene_names.csv", row.names = FALSE)

# Required libraries
library(ggplot2)
library(ggrepel)

# Ensure gene_name is available in the DGE results
dge_results$Significance <- ifelse(
  dge_results$p_val_adj < 0.05 & abs(dge_results$avg_log2FC) > 0.1,
  "Significant",
  "Not Significant"
)

# Create the volcano plot with gene names
volcano_plot <- ggplot(dge_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significance), alpha = 0.6) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  labs(
    title = "Volcano Plot: B-in Cells (LOY vs Non-LOY)",
    x = "Average Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_classic() +  # Use a white background theme
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

# Highlight top genes (optional, based on significance and p-value ranking)
top_genes <- head(dge_results[order(dge_results$p_val_adj), ], 10)  # Top 10 by p-value
volcano_plot <- volcano_plot +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name),  # Use the gene_name column for labels
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = 10
  )

# Display the plot
print(volcano_plot)

# Save the volcano plot
ggsave("volcano_plot_Treg_LOY.png", volcano_plot, width = 8, height = 6)

