library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

# List of Seurat file paths for different cell types
seurat_file_paths <- c(
  "male_B_cells_191124_2.rds",
  "male_CD8T_cells_191124_2.rds",
  "male_CD4T_cells_191124_2.rds",
  "male_NK_cells_191124_2.rds",
  "male_mono_cells_191124_2.rds"
)

# Load gene list and gene mapping (do this once)
gene_list <- read.table("longest_transcripts_unique.txt", header = TRUE, sep = "\t")
non_chrY_genes <- gene_list[gene_list$chromosome != "chrY", "gene"]
gene_mapping <- read.csv("Onek1k_mart_export.txt", header = FALSE, sep = ",", col.names = c("gene_id", "gene_name"))

# XIST Ensemble ID
xist_ensemble_id <- "ENSG00000229807"

# Loop through each Seurat file
for (seurat_file in seurat_file_paths) {

  # Load the Seurat object
  seurat_object <- readRDS(seurat_file)

  # Check if XIST gene is present in the dataset
  if (!(xist_ensemble_id %in% rownames(seurat_object))) {
    cat("Warning: XIST gene not found in dataset:", seurat_file, "\n")
    next
  }

  # Identify cells with XIST expression
  xist_expression <- FetchData(seurat_object, vars = xist_ensemble_id, slot = "counts")
  xist_cells <- rownames(xist_expression)[xist_expression[, xist_ensemble_id] > 0]

  # Print total number of cells before excluding XIST-expressing cells
  cat("\nBefore excluding cells expressing XIST for", seurat_file, "\n")
  print(table(seurat_object$LOY_status))

  # Exclude cells expressing XIST
  seurat_object <- subset(seurat_object, cells = setdiff(colnames(seurat_object), xist_cells))

  # Remove XIST from features
  cat("\n--- Removing XIST from Features ---\n")
  seurat_object <- subset(seurat_object, features = setdiff(rownames(seurat_object), xist_ensemble_id))

  # Re-scale the data to eliminate residual effects
  cat("\n--- Re-scaling Data After Feature Removal ---\n")
  seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))

  # Print total number of cells after excluding XIST-expressing cells
  cat("\nAfter excluding cells expressing XIST for", seurat_file, "\n")
  print(table(seurat_object$LOY_status))

  # Get unique classifications from $new_classification
  unique_classifications <- unique(seurat_object$new_classification)

  # Loop over each unique classification
  for (classification in unique_classifications) {
    cat("\nProcessing classification:", classification, "in", seurat_file, "\n")

    # Subset Seurat object to cells with the current classification
    classification_subset <- subset(seurat_object, subset = new_classification == classification)

    # Set identity to LOY_status for differential analysis
    classification_subset$LOY_status <- as.character(classification_subset$LOY_status)
    Idents(classification_subset) <- "LOY_status"

    # Filter the genes (non-chrY genes)
    filtered_subset <- subset(classification_subset, features = non_chrY_genes)

    # Perform DGE between LOY_status = 1 (LOY cells) vs LOY_status = 0 (non-LOY cells)
    dge_results <- FindMarkers(
      object = filtered_subset,
      ident.1 = "1",  # LOY_status = 1 (LOY cells)
      ident.2 = "0",  # LOY_status = 0 (non-LOY cells)
      test.use = "MAST",
      verbose = FALSE
    )

    # Add gene_id for mapping
    dge_results <- cbind(gene_id = rownames(dge_results), dge_results)

    # Merge with gene names
    dge_results <- merge(dge_results, gene_mapping, by = "gene_id", all.x = TRUE)

    # Save results to a CSV file
    output_file <- paste0("DGE_results_without_XIST_", gsub(".rds", "", basename(seurat_file)), "_", classification, "_LOY_status.csv")
    write.csv(dge_results, file = output_file, row.names = FALSE)

    # Create volcano plot
    dge_results$Significance <- ifelse(
      dge_results$p_val_adj < 0.05 & abs(dge_results$avg_log2FC) > 0.1,
      "Significant",
      "Not Significant"
    )

    volcano_plot <- ggplot(dge_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      geom_point(aes(color = Significance), alpha = 0.6) +
      scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
      labs(
        title = paste("Volcano Plot:", classification, "-", gsub(".rds", "", basename(seurat_file))),
        x = "Average Log2 Fold Change",
        y = "-Log10 Adjusted P-Value"
      ) +
      theme_classic() +
      theme(
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )

    # Highlight top genes
    top_genes <- head(dge_results[order(dge_results$p_val_adj), ], 10)
    volcano_plot <- volcano_plot +
      geom_text_repel(
        data = top_genes,
        aes(label = gene_name),
        size = 3,
        box.padding = 0.3,
        point.padding = 0.3,
        max.overlaps = 10
      )

    # Save and display the plot
    volcano_plot_file <- paste0("volcano_plot_without_XIST_", gsub(".rds", "", basename(seurat_file)), "_", classification, "_LOY_status.png")
    ggsave(volcano_plot_file, volcano_plot, width = 8, height = 6)
    print(volcano_plot)
  }
}

