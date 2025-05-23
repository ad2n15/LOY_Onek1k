# Load required libraries
library(Seurat)
library(ggplot2)
library(magrittr)
library(tibble)



male_cells <- readRDS("male_CD4T_cells_191124_2.rds")

# Step 1: Load the PHATE results
phate_results <- read.csv("phate_result_CD4T.csv", row.names = 1)

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

saveRDS(male_cells, file = "male_CD4T_cells_191124_3.rds")


# Step 4: Visualize the PHATE results using DimPlot
p1 <- DimPlot(male_cells, reduction = "PHATE", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +

  theme(
    plot.title = element_text(size = 16, face = "bold"), # Bold and increase title size
    axis.title = element_text(size = 14, face = "bold"), # Bold and increase axis titles
    axis.text = element_text(size = 12, face = "bold"), # Bold and increase axis tick labels
    legend.title = element_text(size = 12, face = "bold"), # Bold and increase legend title
    legend.text = element_text(size = 10, face = "bold") # Bold and increase legend text
  )

ggsave(filename = "UMAP_Clustering_Male_CD4T_Cells_surat_clusters_phate.png", plot = p1,
       width = 8, height = 6, dpi = 300)


p1 <- DimPlot(male_cells, reduction = "PHATE", group.by = "new_classification") +

  theme(
    plot.title = element_text(size = 16, face = "bold"), # Bold and increase title size
    axis.title = element_text(size = 14, face = "bold"), # Bold and increase axis titles
    axis.text = element_text(size = 12, face = "bold"), # Bold and increase axis tick labels
    legend.title = element_text(size = 12, face = "bold"), # Bold and increase legend title
    legend.text = element_text(size = 10, face = "bold") # Bold and increase legend text
  )

ggsave(filename = "UMAP_Clustering_Male_CD4T_Cells_new_classification_phate.png", plot = p1,
       width = 8, height = 6, dpi = 300)


# Optional Step 5: Access PHATE embeddings for further analysis
phate_embeddings <- Embeddings(male_cells, reduction = "PHATE")
head(phate_embeddings)

# Example: Use PHATE embeddings for clustering (if applicable)
# Use clustering or other Seurat functions as required

###############

library(slingshot)



male_cells <- subset(male_cells, seurat_clusters != 3 & seurat_clusters != 9)

em <- Embeddings(male_cells, reduction = "PHATE")
#sds <- slingshot(em, clusterLabels = male_cells$seurat_clusters, end.clus = c(3,8,9))


sds <- slingshot(em, clusterLabels = male_cells$new_classification, start.clus = "CD4nc", end.clus = "Treg")


########
i <- 1
curve_1 <- slingCurves(sds)[[i]]
pt <- curve_1$lambda %>% as.data.frame() %>% set_names("pt")
subdata <- AddMetaData(male_cells, pt, col.name="pt")
curve_1 <- curve_1$s[curve_1$ord, 1:2]
colnames(curve_1) <- c("PHATE_1", "PHATE_2")

library(ggplot2)
library(viridis)

# Prepare data for plotting
phate_data <- subdata@reductions$PHATE@cell.embeddings %>%
    as.data.frame() %>%
    setNames(c("PHATE_1", "PHATE_2"))
phate_data$pt <- subdata@meta.data$pt

# Create the plot
p <- ggplot(phate_data, aes(x = PHATE_2, y = PHATE_1, color = pt)) +
    geom_point(size = 0.5) +
    geom_path(data = as.data.frame(curve_1), aes(x = PHATE_2, y = PHATE_1),
              linewidth = 1, inherit.aes = FALSE, color = "black") +
    xlab("PHATE 1") +
    ylab("PHATE 2") +
    scale_color_viridis_c(name = "Pseudotime") +
    theme_bw() + # White background
    theme(
        plot.title = element_text(size = 16, face = "bold"), # Bold and larger title
        axis.title = element_text(size = 14, face = "bold"), # Bold axis labels
        axis.text = element_text(size = 12, face = "bold"), # Bold and larger axis tick labels
        legend.title = element_text(size = 12, face = "bold"), # Bold legend title
        legend.text = element_text(size = 10, face = "bold") # Bold legend text
    )

# Save the plot
ggsave("cells_curve1_CD4T.png", p, height = 5.5, width = 8.5)




#######################################################################################################




library(dplyr)
library(tidyr)
library(Seurat)

# Assuming `sds` is your Slingshot object and `male_cells` is your Seurat object with male cells

# Define the curve index for Slingshot and extract pseudotime
i <- 1
curve_2 <- slingCurves(sds)[[i]]
lambda <- curve_2$lambda %>% as.data.frame() %>% set_names("lambda")

# Add pseudotime data to the Seurat object
pseudotime <- slingPseudotime(sds)[[i]]
subdata <- AddMetaData(male_cells, pseudotime, col.name = "pseudotime")

# Add lambda to the metadata
subdata <- AddMetaData(subdata, lambda, col.name = "lambda")

# Ensure pseudotime and lambda are numeric
subdata@meta.data$pseudotime <- as.numeric(as.character(subdata@meta.data$pseudotime))
subdata@meta.data$lambda <- as.numeric(as.character(subdata@meta.data$lambda))

# Create quantile groups based on lambda
subdata@meta.data$quantile <- cut(
  subdata@meta.data$lambda,
  breaks = quantile(subdata@meta.data$lambda, probs = 0:6/6, na.rm = TRUE),
  include.lowest = TRUE,
  labels = FALSE
)

saveRDS(subdata, file = "male_CD4T_cells_191124_4.rds")

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
ggsave("PHATE_with_quantiles_CD4T_cells.png", p, height = 5.5, width = 8.5)

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
ggsave("PHATE_with_celltype_CD4T_cells.png", p, height = 5.5, width = 8.5)

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



library(data.table)
library(lme4)
library(parallel)
library(ggplot2)
library(ggrepel)

# Initialize results tables
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

# Parallelize the computation
process_gene <- function(gene, data, use_quadratic = FALSE) {
  if (gene %in% colnames(data)) {
    tryCatch({
      if (!use_quadratic) {
        # Original test
        formula_null <- as.formula(paste(
          gene, "~ quantile + LOY_status + age.x + pca1 + pca2 + pca3 +
          pca4 + pca5 + pca6 + pca7 + pca8 + pca9 + pca10 + (1 | donor_id)"
        ))
        formula_augmented <- as.formula(paste(
          gene, "~ quantile + LOY_status + age.x + pca1 + pca2 + pca3 +
          pca4 + pca5 + pca6 + pca7 + pca8 + pca9 + pca10 + (1 | donor_id) + quantile * LOY_status"
        ))

        fit0 <- lmer(formula_null, data = data, REML = FALSE)
        fit1 <- lmer(formula_augmented, data = data, REML = FALSE)

        anova_result <- anova(fit1, fit0)
        list(
          Gene = gene,
          Interaction_estimate = fixef(fit1)["quantile:LOY_status"],
          chisq_stat = anova_result$Chisq[2],
          Interaction_pvalue = anova_result$`Pr(>Chisq)`[2]
        )
      } else {
        # Quadratic test
        data$quantile_2 <- data$quantile^2
        formula_null_quad <- as.formula(paste(
          gene, "~ quantile + quantile_2 + LOY_status + age.x + pca1 + pca2 + pca3 +
          pca4 + pca5 + pca6 + pca7 + pca8 + pca9 + pca10 + (1 | donor_id)"
        ))
        formula_augmented_quad <- as.formula(paste(
          gene, "~ quantile + quantile_2 + LOY_status + age.x + pca1 + pca2 + pca3 +
          pca4 + pca5 + pca6 + pca7 + pca8 + pca9 + pca10 + (1 | donor_id) +
          quantile * LOY_status + quantile_2 * LOY_status"
        ))

        fit0_quad <- lmer(formula_null_quad, data = data, REML = FALSE)
        fit1_quad <- lmer(formula_augmented_quad, data = data, REML = FALSE)

        anova_result_quad <- anova(fit1_quad, fit0_quad)
        list(
          Gene = gene,
          Interaction_estimate_linear = fixef(fit1_quad)["quantile:LOY_status"],
          Interaction_estimate_quadratic = fixef(fit1_quad)["quantile_2:LOY_status"],
          chisq_stat = anova_result_quad$Chisq[2],
          Interaction_pvalue = anova_result_quad$`Pr(>Chisq)`[2]
        )
      }
    }, error = function(e) NULL)
  } else {
    NULL
  }
}

# Use parallel processing for efficiency
num_cores <- detectCores() - 1
original_results_list <- mclapply(
  filtered_gene_names,
  process_gene,
  data = merged_data3,
  use_quadratic = FALSE,
  mc.cores = num_cores
)
original_results <- rbindlist(original_results_list, fill = TRUE)

quadratic_results_list <- mclapply(
  filtered_gene_names,
  process_gene,
  data = merged_data3,
  use_quadratic = TRUE,
  mc.cores = num_cores
)
quadratic_results <- rbindlist(quadratic_results_list, fill = TRUE)

# Combine results
combined_results <- merge(
  original_results, quadratic_results, by = "Gene", all = TRUE, suffixes = c("_original", "_quadratic")
)

# Add gene names using merge
gene_mapping <- fread("Onek1k_mart_export.txt", header = FALSE)
setnames(gene_mapping, c("ensembl_id", "gene_name"))
combined_results <- merge(combined_results, gene_mapping, by.x = "Gene", by.y = "ensembl_id", all.x = TRUE)

# FDR correction and odds ratio calculation
combined_results[, adjusted_pvalue_original := p.adjust(Interaction_pvalue_original, method = "fdr")]
combined_results[, adjusted_pvalue_quadratic := p.adjust(Interaction_pvalue_quadratic, method = "fdr")]
combined_results[, log10_adjusted_pvalue_original := -log10(pmax(adjusted_pvalue_original, 1e-300))]
combined_results[, log10_adjusted_pvalue_quadratic := -log10(pmax(adjusted_pvalue_quadratic, 1e-300))]
combined_results[, Interaction_odds_ratio := exp(Interaction_estimate)]
combined_results[, log2_Interaction_odds_ratio := log2(Interaction_odds_ratio)]
combined_results[, Interaction_quad_odds_ratio := exp(Interaction_estimate_quadratic)]
combined_results[, log2_Interaction_quad_odds_ratio := log2(Interaction_quad_odds_ratio)]

# Save results
fwrite(combined_results, "Combined_LOY_status_results_with_gene_names_CD4T_cells.csv")

# Visualization
top_10_genes <- combined_results[order(abs(Interaction_estimate), decreasing = TRUE)][1:10]
volcano_plot <- ggplot(combined_results, aes(x = log2_Interaction_odds_ratio, y = log10_adjusted_pvalue_original)) +
  geom_point(aes(color = adjusted_pvalue_original < 0.05), size = 3, alpha = 0.7) +
  scale_color_manual(values = c("gray", "red"), name = "Significance") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.8) +
  geom_text_repel(data = top_10_genes, aes(label = gene_name), size = 3, box.padding = 0.2) +
  labs(title = "Volcano Plot (log2(OR))", x = "Log2(Odds Ratio)", y = "-log10(Adjusted p-value)") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("volcano_plot_with_top10_genes_repelled_CD4T_cells_linear.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

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
ggsave("classification_LOY_status_plot_CD4T_cells.png", classification_plot, width = 8, height = 6, dpi = 300)

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
ggsave("quantile_LOY_status_plot_CD4T_Cells.png", quantile_plot, width = 8, height = 6, dpi = 300)

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
      file_name <- paste0("expression_across_quantiles_with_fit_lines_CD4T_", gene, ".png")
      ggsave(file_name, plot = expression_plot, width = 10, height = 8, dpi = 300, bg = "white") # bg ensures a white background
      
      cat("Plot saved for gene:", gene, "\n")
    } else {
      cat("The gene", gene, "is not found in the scaled data.\n")
    }
  } else {
    cat("The gene", gene, "is not found in the counts data.\n")
  }
}

