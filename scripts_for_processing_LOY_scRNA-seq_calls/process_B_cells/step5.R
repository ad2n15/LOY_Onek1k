# Load required libraries
library(Seurat)
library(ggplot2)
library(magrittr)
library(tibble)



male_cells <- readRDS("male_B_cells_191124_2.rds")

# Step 1: Load the PHATE results
phate_results <- read.csv("phate_result.csv", row.names = 1)

# Check the structure of the PHATE results
head(phate_results)

# Step 2: Ensure barcodes in PHATE results match those in Seurat object
# Replace `male_cells` with your actual Seurat object name
if (!all(rownames(phate_results) %in% colnames(male_cells))) {
  stop("Some barcodes in PHATE results do not match the Seurat object. Check for mismatches.")
}

# Optional: Subset PHATE results to match the Seurat object (if needed)
phate_results <- phate_results[colnames(male_cells), ]

# Step 3: Add PHATE results to the Seurat object as a dimensional reduction
phate_matrix <- as.matrix(phate_results) # Convert PHATE results to a matrix
male_cells[["PHATE"]] <- CreateDimReducObject(
  embeddings = phate_matrix,
  key = "PHATE_",
  assay = DefaultAssay(male_cells)
)


# Verify the updated Seurat object
print(male_cells)

saveRDS(male_cells, file = "male_B_cells_191124_3.rds")


# Step 4: Visualize the PHATE results using DimPlot
p1 <- DimPlot(male_cells, reduction = "PHATE", group.by = "seurat_clusters") +

  theme(
    plot.title = element_text(size = 16, face = "bold"), # Bold and increase title size
    axis.title = element_text(size = 14, face = "bold"), # Bold and increase axis titles
    axis.text = element_text(size = 12, face = "bold"), # Bold and increase axis tick labels
    legend.title = element_text(size = 12, face = "bold"), # Bold and increase legend title
    legend.text = element_text(size = 10, face = "bold") # Bold and increase legend text
  )


ggsave(filename = "UMAP_Clustering_Male_B_Cells_surat_clusters_phate.png", plot = p1,
       width = 8, height = 6, dpi = 300)

p1 <- DimPlot(male_cells, reduction = "PHATE", group.by = "new_classification") +

  theme(
    plot.title = element_text(size = 16, face = "bold"), # Bold and increase title size
    axis.title = element_text(size = 14, face = "bold"), # Bold and increase axis titles
    axis.text = element_text(size = 12, face = "bold"), # Bold and increase axis tick labels
    legend.title = element_text(size = 12, face = "bold"), # Bold and increase legend title
    legend.text = element_text(size = 10, face = "bold") # Bold and increase legend text
  )








ggsave(filename = "UMAP_Clustering_Male_B_Cells_new_classification_phate.png", plot = p1,
       width = 8, height = 6, dpi = 300)


# Optional Step 5: Access PHATE embeddings for further analysis
phate_embeddings <- Embeddings(male_cells, reduction = "PHATE")
head(phate_embeddings)

# Example: Use PHATE embeddings for clustering (if applicable)
# Use clustering or other Seurat functions as required

###############

library(slingshot)


em <- Embeddings(male_cells, reduction = "PHATE")
sds <- slingshot(em, clusterLabels = male_cells$new_classification, start.clus = "B-in", end.clus = "B-mem")


#######################################################################################################




library(dplyr)
library(tidyr)
library(Seurat)

# Assuming `sds` is your Slingshot object and `male_cells` is your Seurat object with male cells

# Define the curve index for Slingshot and extract pseudotime
i <- 1
curve_1 <- slingCurves(sds)[[i]]
lambda <- curve_1$lambda %>% as.data.frame() %>% set_names("lambda")
subdata <- AddMetaData(male_cells, lambda, col.name = "lambda")


# Add pseudotime data to the Seurat object
pseudotime <- slingPseudotime(sds, na = FALSE)
subdata <- AddMetaData(subdata, pseudotime, col.name = "pseudotime")
#################
curve_1 <- curve_1$s[curve_1$ord, 1:2]
colnames(curve_1) <- c("PHATE_1", "PHATE_2")

library(ggplot2)
library(viridis)

# Prepare data for plotting
phate_data <- subdata@reductions$PHATE@cell.embeddings %>%
    as.data.frame() %>%
    setNames(c("PHATE_1", "PHATE_2"))
phate_data$lambda <- subdata@meta.data$lambda

# Create the plot
p <- ggplot(phate_data, aes(x = PHATE_2, y = PHATE_1, color = lambda)) +
    geom_point(size = 0.5) +
    geom_path(data = as.data.frame(curve_1), aes(x = PHATE_2, y = PHATE_1),
              linewidth = 1, inherit.aes = FALSE, color = "black") +
    xlab("PHATE 1") +
    ylab("PHATE 2") +
    scale_color_viridis_c(name = "Pseudotime") +
    theme_bw() +
    theme(
        plot.title = element_text(size = 16, face = "bold"), # Bold and larger title
        axis.title = element_text(size = 14, face = "bold"), # Bold axis labels
        axis.text = element_text(size = 12, face = "bold"), # Bold and larger axis tick labels
        legend.title = element_text(size = 12, face = "bold"), # Bold legend title
        legend.text = element_text(size = 10, face = "bold") # Bold legend text
    )


# Save the plot
ggsave("cells_curve1_B_cells.png", p, height = 5.5, width = 8.5)




# Ensure pseudotime and lambda are numeric
#subdata@meta.data$pseudotime <- as.numeric(as.character(subdata@meta.data$pseudotime))
#subdata@meta.data$lambda <- as.numeric(as.character(subdata@meta.data$lambda))

# Create quantile groups based on lambda
subdata@meta.data$quantile <- cut(
  subdata@meta.data$pseudotime,
  breaks = quantile(subdata@meta.data$pseudotime, probs = 0:6/6, na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)


saveRDS(subdata, file = "male_B_cells_191124_4.rds")


)

# Add quantile to phate_data for plotting
phate_data$quantile <- subdata@meta.data$quantile

# Create the plot (color by quantile)
p <- ggplot(phate_data, aes(x = PHATE_2, y = PHATE_1, color = factor(quantile))) +
  geom_point(size = 0.5) +  # Scatter plot with color by quantile
  geom_path(data = as.data.frame(curve_1), aes(x = PHATE_2, y = PHATE_1),
            linewidth = 1, inherit.aes = FALSE, color = "black") +  # Slingshot trajectory (optional)
  xlab("PHATE 1") +
  ylab("PHATE 2") +
  scale_color_viridis_d(name = "Quantile") +  # Discrete color scale for quantiles
  theme_bw() +  # Clean background
  theme(
    plot.title = element_text(size = 16, face = "bold"),  # Bold and larger title
    axis.title = element_text(size = 14, face = "bold"),  # Bold axis labels
    axis.text = element_text(size = 12, face = "bold"),  # Bold and larger axis tick labels
    legend.title = element_text(size = 12, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 10, face = "bold")  # Bold legend text
  ) +
  labs(
    title = "PHATE Plot with Quantiles",
    subtitle = "Colored by Quantile (Based on Pseudotime)"
  )

# Save the plot
ggsave("PHATE_with_quantiles_B_cells.png", p, height = 5.5, width = 8.5)

################
# plot phate with celltype

# Convert new_classification to a factor
phate_data$new_classification <- as.factor(subdata@meta.data$new_classification)

# Extract unique categories
unique_classes <- unique(phate_data$new_classification)
print(unique_classes)  # Check unique values

# Define base colors (starting with Light Blue & Light Red, then using a color palette)
base_colors <- c("#ADD8E6", "#FFA07A")  # Light Blue & Light Red
extra_colors <- colorspace::rainbow_hcl(length(unique_classes) - length(base_colors))  # More colors if needed

# Assign colors dynamically
custom_colors <- setNames(c(base_colors, extra_colors)[seq_along(unique_classes)], unique_classes)

# Create PHATE plot
p <- ggplot(phate_data, aes(x = PHATE_2, y = PHATE_1, color = new_classification)) +
  geom_point(size = 0.5) +  
  geom_path(data = as.data.frame(curve_1), aes(x = PHATE_2, y = PHATE_1),
            linewidth = 1, inherit.aes = FALSE, color = "black") +  
  xlab("PHATE 1") +
  ylab("PHATE 2") +
  scale_color_manual(name = "Cell Type", values = custom_colors) +  
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold")
  ) +
  labs(
    title = "PHATE Plot with Cell Type",
    subtitle = "Colored by Cell Type"
  )

# Save the plot
ggsave("PHATE_with_celltype_B_cells.png", p, height = 5.5, width = 8.5)

###########


# Select relevant columns for the table
# Ensure that the columns 'age', 'seurat_clusters', 'new_classification', and 'LOY_status' are present in metadata
columns_to_include <- c("lambda", "pseudotime", "age", "seurat_clusters", "new_classification", "quantile", "LOY_status", "donor_id", "pool_number")
metadata_table <- subdata@meta.data %>%
  select(all_of(columns_to_include)) %>%
  as.data.frame()

# Print the first few rows of the metadata table to ensure it looks correct
head(metadata_table)



# Load the gene mapping file to map Ensembl IDs to gene names
gene_mapping <- read.csv("Onek1k_mart_export.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(gene_mapping) <- c("ensembl_id", "gene_name")



# Step 1: Load the gene list
gene_list <- read.table("longest_transcripts_unique.txt", header = TRUE, sep = "\t")
# Step 2: Filter out chrY genes
non_chrY_genes <- gene_list[gene_list$chromosome != "chrY", "gene"]
# Step 3: Subset the Seurat object
# Assuming 'male_cells_B_in' is your Seurat object
subdata <- subset(subdata, features = non_chrY_genes)


# Get the list of highly variable genes
hvg_genes <- rownames(subdata@assays$RNA@data)
length(hvg_genes)
# Re-scale the data to include all genes
subdata <- ScaleData(subdata, features = rownames(subdata))

# Extract all genes from the scaled data
scaled_data <- subdata@assays$RNA@scale.data %>%
  as.data.frame() %>%
  t() %>% # Transpose so cells are rows, genes are columns
  as.data.frame()

# Add cell names as a column
scaled_data$cell <- rownames(scaled_data)

# Combine with metadata
metadata_table2 <- cbind(metadata_table, scaled_data[rownames(metadata_table), , drop = FALSE])

# Load the additional metadata (e.g., 10 PCA components, etc.)
data <- read.csv("onek1k_10pca_981.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# Merge metadata with the additional data
merged_data3 <- left_join(metadata_table2, data, by = "donor_id")




#########################################

# Create a dataframe to store the results
results <- data.frame(
  Gene = character(),
  LOY_status_estimate = numeric(),
  LOY_status_pvalue = numeric(),
  stringsAsFactors = FALSE
)



# Filter genes in hvg_genes based on expression in at least 3 quantiles
filtered_genes <- sapply(hvg_genes, function(gene) {
  # Check if the gene is present in the count matrix
  if (gene %in% rownames(subdata@assays$RNA@counts)) {
    # Extract gene expression counts
    gene_expression <- subdata@assays$RNA@counts[gene, ]
    
    # Get the quantile information for each cell
    quantile_info <- subdata@meta.data$quantile
    
    # Identify quantiles where the gene is expressed (non-zero expression)
    expressed_quantiles <- unique(quantile_info[gene_expression > 0])
    
    # Check if the gene is expressed in at least 3 unique quantiles
    return(length(expressed_quantiles) >= 3)
  } else {
    return(FALSE)  # Exclude genes not present in the dataset
  }
})

# Extract the names of filtered genes
filtered_gene_names <- hvg_genes[filtered_genes]

# Print the filtered gene names
print(filtered_gene_names)


library(ggrepel)   # For adjusting text labels in plots
library(lme4)        # For mixed-effects model
# Initialize results for both tests as data.tables
library(data.table)
original_results <- data.table(
  Gene = character(),
  Interaction_estimate = numeric(),
  chisq_stat = numeric(),
  Interaction_pvalue = numeric()
)

quadratic_results <- data.table(
  Gene = character(),
  Interaction_estimate_linear = numeric(),
  Interaction_estimate_quadratic = numeric(),
  chisq_stat = numeric(),
  Interaction_pvalue = numeric()
)

# Iterate over genes
for (gene in filtered_gene_names) {
  if (gene %in% colnames(merged_data3)) {
    # Original test
    formula_null <- as.formula(paste(
      gene, "~ quantile + LOY_status + age.x + pca1 + pca2 + pca3 +
      pca4 + pca5 + pca6 + pca7 + pca8 + pca9 + pca10 + (1 | donor_id)"
    ))
    formula_augmented <- as.formula(paste(
      gene, "~ quantile + LOY_status + age.x + pca1 + pca2 + pca3 +
      pca4 + pca5 + pca6 + pca7 + pca8 + pca9 + pca10 + (1 | donor_id) + quantile*LOY_status"
    ))

    fit0 <- lmer(formula_null, data = merged_data3, REML = FALSE)
    fit1 <- lmer(formula_augmented, data = merged_data3, REML = FALSE)

    anova_result <- anova(fit1, fit0)
    chisq_stat <- anova_result$Chisq[2]
    p_value <- anova_result$`Pr(>Chisq)`[2]
    estimate <- fixef(fit1)["quantile:LOY_status"]

    original_results <- rbindlist(list(original_results, data.table(
      Gene = gene,
      Interaction_estimate = estimate,
      chisq_stat = chisq_stat,
      Interaction_pvalue = p_value
    )))

    # Quadratic test
    merged_data3$quantile_2 <- merged_data3$quantile^2
    formula_null_quad <- as.formula(paste(
      gene, "~ quantile + quantile_2 + LOY_status + age.x + pca1 + pca2 + pca3 +
      pca4 + pca5 + pca6 + pca7 + pca8 + pca9 + pca10 + (1 | donor_id)"
    ))
    formula_augmented_quad <- as.formula(paste(
      gene, "~ quantile + quantile_2 + LOY_status + age.x + pca1 + pca2 + pca3 +
      pca4 + pca5 + pca6 + pca7 + pca8 + pca9 + pca10 + (1 | donor_id) + 
      quantile*LOY_status + quantile_2*LOY_status"
    ))

    fit0_quad <- lmer(formula_null_quad, data = merged_data3, REML = FALSE)
    fit1_quad <- lmer(formula_augmented_quad, data = merged_data3, REML = FALSE)

    anova_result_quad <- anova(fit1_quad, fit0_quad)
    chisq_stat_quad <- anova_result_quad$Chisq[2]
    p_value_quad <- anova_result_quad$`Pr(>Chisq)`[2]
    estimate_linear <- fixef(fit1_quad)["quantile:LOY_status"]
    estimate_quadratic <- fixef(fit1_quad)["quantile_2:LOY_status"]

    quadratic_results <- rbindlist(list(quadratic_results, data.table(
      Gene = gene,
      Interaction_estimate_linear = estimate_linear,
      Interaction_estimate_quadratic = estimate_quadratic,
      chisq_stat = chisq_stat_quad,
      Interaction_pvalue = p_value_quad
    )))
  } else {
    warning(paste("Gene", gene, "not found in the data"))
  }
}

# Combine results
combined_results <- merge(
  original_results, quadratic_results, by = "Gene", all = TRUE, suffixes = c("_original", "_quadratic")
)

# Add gene names using merge
combined_results <- merge(combined_results, gene_mapping, by.x = "Gene", by.y = "ensembl_id", all.x = TRUE)

# FDR correction for p-values
combined_results[, adjusted_pvalue_original := p.adjust(Interaction_pvalue_original, method = "fdr")]
combined_results[, adjusted_pvalue_quadratic := p.adjust(Interaction_pvalue_quadratic, method = "fdr")]

# Log10 transformation for adjusted p-values
combined_results[, log10_adjusted_pvalue_original := -log10(pmax(adjusted_pvalue_original, 1e-300))]
combined_results[, log10_adjusted_pvalue_quadratic := -log10(pmax(adjusted_pvalue_quadratic, 1e-300))]

# Convert LOY_status_estimate to Odds Ratio (OR)
combined_results$Interaction_odds_ratio <- exp(combined_results$Interaction_estimate)
# Convert the Odds Ratio to log2 scale
combined_results$log2_Interaction_odds_ratio <- log2(combined_results$Interaction_odds_ratio)

# Convert LOY_status_estimate to Odds Ratio (OR)
combined_results$Interaction_quad_odds_ratio <- exp(combined_results$Interaction_estimate_quadratic)
# Convert the Odds Ratio to log2 scale
combined_results$log2_Interaction_quad_odds_ratio <- log2(combined_results$Interaction_quad_odds_ratio)


# Load necessary package
library(ggrepel)

# Load the gene mapping file
gene_mapping <- read.csv("Onek1k_mart_export.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(gene_mapping) <- c("ensembl_id", "gene_name")

# Add gene names to the results dataframe
combined_results$GeneName <- sapply(combined_results$Gene, function(x) {
  # Find the corresponding gene name from the mapping
  gene_name <- gene_mapping$gene_name[gene_mapping$ensembl_id == x]
  return(ifelse(length(gene_name) > 0, gene_name, NA))  # return NA if no gene name is found
})


# View the results after the conversion
results = combined_results
head(results)


# Save combined results
fwrite(combined_results, "Combined_LOY_status_results_with_gene_names_B_cells.csv")



# Extract top 10 genes with the highest absolute effect size
top_10_genes <- results[order(abs(results$Interaction_estimate), decreasing = TRUE), ][1:10, ]

# Volcano plot with log2(odds_ratio)
volcano_plot <- ggplot(results, aes(x = log2_Interaction_odds_ratio, y = log10_adjusted_pvalue_original)) +
  geom_point(aes(color = adjusted_pvalue_original < 0.05), size = 3, alpha = 0.7) +  # Color based on significance
  scale_color_manual(values = c("gray", "red"), name = "Significance", labels = c("Not Significant", "Significant")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.8) +  # Updated line width
  
  # Use ggrepel to automatically adjust text position for top 10 genes
  geom_text_repel(data = top_10_genes,
                  aes(label = GeneName),
                  size = 3, 
                  box.padding = 0.2,  # Space around the text
                  max.overlaps = 10,  # Limit the number of overlaps
                  color = "black",   # Set text color
                  force = 2,         # Increase force to move text further
                  nudge_y = 0.1,     # Slight nudge in the y direction
                  nudge_x = 0.2) +   # Slight nudge in the x direction
  labs(
    title = "Volcano Plot (log2(OR))",
    x = "Log2(Odds Ratio)",
    y = "-log10(Adjusted p-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    text = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA),  # Panel background white
    plot.background = element_rect(fill = "white", color = NA)   # Plot background white
  )

# Save the plot as a PNG file
ggsave("volcano_plot_with_top10_genes_repelled_B_cells_linear.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

# Print the top 10 genes with their names and effect sizes
print(top_10_genes)


# Extract top 10 genes with the highest absolute effect size
top_10_genes <- results[order(abs(results$Interaction_estimate_quadratic), decreasing = TRUE), ][1:10, ]

# Volcano plot with log2(odds_ratio)
volcano_plot <- ggplot(results, aes(x = log2_Interaction_quad_odds_ratio, y = log10_adjusted_pvalue_quadratic)) +
  geom_point(aes(color = adjusted_pvalue_quadratic < 0.05), size = 3, alpha = 0.7) +  # Color based on significance
  scale_color_manual(values = c("gray", "red"), name = "Significance", labels = c("Not Significant", "Significant")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.8) +  # Updated line width

  # Use ggrepel to automatically adjust text position for top 10 genes
  geom_text_repel(data = top_10_genes,
                  aes(label = GeneName),
                  size = 3,
                  box.padding = 0.2,  # Space around the text
                  max.overlaps = 10,  # Limit the number of overlaps
                  color = "black",   # Set text color
                  force = 2,         # Increase force to move text further
                  nudge_y = 0.1,     # Slight nudge in the y direction
                  nudge_x = 0.2) +   # Slight nudge in the x direction
  labs(
    title = "Volcano Plot (log2(OR))",
    x = "Log2(Odds Ratio)",
    y = "-log10(Adjusted p-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    text = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA),  # Panel background white
    plot.background = element_rect(fill = "white", color = NA)   # Plot background white
  )

# Save the plot as a PNG file
ggsave("volcano_plot_with_top10_genes_repelled_B_cells_quadratic.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)



############################

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# 1. Summarize data for new_classification
classification_summary <- subdata@meta.data %>%
  group_by(new_classification, LOY_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(new_classification) %>%
  mutate(percentage = (count / sum(count)) * 100)

# 2. Create the plot for new_classification
classification_plot <- ggplot(classification_summary, aes(x = new_classification, y = count, fill = as.factor(LOY_status))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Counts and Percentages of LOY_status in each new_classification", 
       y = "Cell Counts", x = "New Classification", fill = "LOY_status") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the classification plot
ggsave("classification_LOY_status_plot_B_cells.png", classification_plot, width = 8, height = 6, dpi = 300)

# 3. Summarize data for quantile
quantile_summary <- subdata@meta.data %>%
  group_by(quantile, LOY_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(quantile) %>%
  mutate(percentage = (count / sum(count)) * 100)

# 4. Create the plot for quantile
quantile_plot <- ggplot(quantile_summary, aes(x = quantile, y = count, fill = as.factor(LOY_status))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Counts and Percentages of LOY_status in each Quantile", 
       y = "Cell Counts", x = "Quantile", fill = "LOY_status") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the quantile plot
ggsave("quantile_LOY_status_plot_B_Cells.png", quantile_plot, width = 8, height = 6, dpi = 300)

#############


# List of genes to analyze
genes_of_interest <- c("ENSG00000234745", "ENSG00000229807", "ENSG00000100721") # Replace with your gene list

# Loop through each gene in the list
for (gene in genes_of_interest) {

  # Ensure the gene exists in the data
  if (gene %in% rownames(subdata@assays$RNA@counts)) {
    
    # Extract expression data for the gene
    subdata@meta.data$gene_expression <- FetchData(subdata, vars = gene)[, 1]
    
    # Classify cells as expressing or not
    subdata@meta.data$expressing <- ifelse(subdata@meta.data$gene_expression > 0, "Expressing", "Non-Expressing")
    
    # Summarize data by quantile and LOY_status
    gene_expression_summary <- subdata@meta.data %>%
      group_by(quantile, LOY_status, expressing) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(quantile, LOY_status) %>%
      mutate(percentage = (count / sum(count)) * 100) %>%
      filter(expressing == "Expressing") # Focus on cells that express the gene
    
    # Print the summary table for the gene
    cat("Summary for gene:", gene, "\n")
    print(gene_expression_summary)
    
    # Check if the gene is in the scaled data
    if (gene %in% rownames(subdata@assays$RNA@scale.data)) {
      
      # Extract the gene expression data for the gene of interest
      gene_expression <- subdata@assays$RNA@scale.data[gene, ]
      
      # Add expression data to metadata for easier plotting
      subdata@meta.data$gene_expression <- gene_expression
      
      # Ensure 'quantile' is treated as a factor
      subdata@meta.data$quantile <- factor(subdata@meta.data$quantile, levels = 1:6)
      
      # Create the plot with 12 boxplots (6 for LOY_status=1, 6 for LOY_status=0)
      expression_plot <- ggplot(subdata@meta.data, aes(x = quantile, y = gene_expression, fill = as.factor(LOY_status))) +
        geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.7) +  # Add boxplot and remove outliers
        geom_smooth(aes(group = LOY_status, color = as.factor(LOY_status)),  # Add fit lines
                    method = "lm", se = FALSE, size = 1.2) +  # Use linear fit, no confidence interval
        facet_wrap(~ LOY_status, ncol = 2) +  # Create two facets for LOY_status 0 and 1
        labs(
          title = paste("Expression of", gene, "Across Quantiles by LOY Status"),
          x = "Quantile",
          y = "Gene Expression",
          fill = "LOY Status",
          color = "LOY Status"  # Legend for fit lines
        ) +
        theme_minimal(base_size = 14) +
        theme(
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA), # Ensures the plot has a white background
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          legend.position = "top"
        ) +
        scale_fill_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e")) +  # Customize fill colors
        scale_color_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e"))   # Customize fit line colors
      
      # Save the plot using ggsave()
      file_name <- paste0("expression_across_quantiles_with_fit_lines_B_", gene, ".png")
      ggsave(file_name, plot = expression_plot, width = 10, height = 8, dpi = 300, bg = "white") # bg ensures a white background
      
      cat("Plot saved for gene:", gene, "\n")
    } else {
      cat("The gene", gene, "is not found in the scaled data.\n")
    }
  } else {
    cat("The gene", gene, "is not found in the counts data.\n")
  }
}

