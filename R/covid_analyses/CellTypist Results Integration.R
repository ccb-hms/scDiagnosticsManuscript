# ----------------------------------------------
# COVID-19 - CellTypist Results Integration
# ----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)

# Souce files
source("R/auxiliary/EnvironmentSetupCellTypist.R")

# ----------------------------------------------

# Ensure Python environment is set up
cat("Ensuring Python environment is available...\n")
environmentSetupCellTypist()

# Configuration
ORIGINAL_COVID_FILE <- "data/covid_data.rds"
ANNOTATED_QUERY_FILE <- "data/annotated_query_data.h5ad"
FINAL_OUTPUT_FILE <- "data/covid_data_annotated.rds"

cat("Loading original COVID-19 SCE object...\n")
covid_data_original <- readRDS(ORIGINAL_COVID_FILE)

cat("Loading CellTypist annotated results...\n")
adata_annotated <- reticulate::import("scanpy")$read_h5ad(ANNOTATED_QUERY_FILE)

cat("Converting annotated AnnData to SCE (annotations only)...\n")
sce_annotations_only <- AnnData2SCE(adata_annotated)

cat("Original COVID data dimensions:", dim(covid_data_original), "\n")
cat("Annotations data dimensions:", dim(sce_annotations_only), "\n")

# Check cell order matches
cell_names_original <- colnames(covid_data_original)
cell_names_annotated <- colnames(sce_annotations_only)

if (!identical(cell_names_original, cell_names_annotated)) {
  warning("Cell order doesn't match exactly, attempting to align...")
  
  # Find common cells
  common_cells <- intersect(cell_names_original, cell_names_annotated)
  cat("Common cells found:", length(common_cells), "\n")
  cat("Original cells:", length(cell_names_original), "\n")
  cat("Annotated cells:", length(cell_names_annotated), "\n")
  
  if(length(common_cells) == 0) {
    stop("No common cells found between original and annotated data!")
  }
  
  # Align both objects to common cells
  covid_data_original <- covid_data_original[, common_cells]
  sce_annotations_only <- sce_annotations_only[, common_cells]
  
  cat("After alignment - both datasets have", ncol(covid_data_original), "cells\n")
}

# Extract annotation columns from CellTypist
annotation_coldata <- colData(sce_annotations_only)
cat("CellTypist annotation columns:\n")
print(colnames(annotation_coldata))

# Add annotations to original COVID data
cat("Adding CellTypist annotations to original COVID data...\n")
original_coldata <- colData(covid_data_original)
combined_coldata <- cbind(original_coldata, annotation_coldata)

# Update the original SCE object with annotations
colData(covid_data_original) <- combined_coldata

cat("Final colData columns:\n")
print(colnames(colData(covid_data_original)))

# Summary of annotations
cat("\n=== CellTypist Annotation Summary ===\n")

if ("predicted_labels" %in% colnames(colData(covid_data_original))) {
  cat("Predicted cell types:\n")
  pred_table <- table(colData(covid_data_original)$predicted_labels)
  print(pred_table)
}

if ("majority_voting" %in% colnames(colData(covid_data_original))) {
  cat("\nMajority voting cell types:\n")
  maj_table <- table(colData(covid_data_original)$majority_voting)
  print(maj_table)
}

if ("conf_score" %in% colnames(colData(covid_data_original))) {
  conf_scores <- colData(covid_data_original)$conf_score
  cat("\nConfidence score summary:\n")
  print(summary(conf_scores))
  
  low_conf_count <- sum(conf_scores < 0.5, na.rm = TRUE)
  low_conf_pct <- round(100 * low_conf_count / length(conf_scores), 1)
  cat("Low confidence cells (<0.5):", low_conf_count, 
      paste0("(", low_conf_pct, "%)\n"))
}

# Save final annotated object
cat("Saving final annotated SCE object...\n")
saveRDS(covid_data_original, FINAL_OUTPUT_FILE)

cat("\n=== Integration Complete ===\n")
cat("Final annotated object saved as:", FINAL_OUTPUT_FILE, "\n")
cat("Final object dimensions:", dim(covid_data_original), "\n")
cat("Total colData columns:", ncol(colData(covid_data_original)), "\n")