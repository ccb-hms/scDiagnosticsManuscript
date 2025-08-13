# ----------------------------------------------
# COVID-19 - CellTypist Results Integration
# ----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)

# Source files
source("R/auxiliary/environmentSetupCellTypist.R")

# ----------------------------------------------

# Ensure Python environment is set up
cat("Ensuring Python environment is available...\n")
environmentSetupCellTypist()

# Configuration
ORIGINAL_COVID_FILE <- "data/covid/covid_data.rds"
ANNOTATED_QUERY_FILE <- "data/covid/annotated_query_data.h5ad"

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

# Identify which CellTypist annotation columns are actually present
annotation_coldata <- colData(sce_annotations_only)
available_cols <- colnames(annotation_coldata)

# Define possible CellTypist columns and check which ones exist
possible_celltypist_cols <- c("predicted_labels", "majority_voting", "conf_score")
present_celltypist_cols <- intersect(possible_celltypist_cols, available_cols)

cat("Available annotation columns in data:\n")
print(available_cols)

if(length(present_celltypist_cols) == 0) {
  stop("No CellTypist annotation columns found! Expected at least 'predicted_labels'.")
}

cat("CellTypist columns found:\n")
print(present_celltypist_cols)

# Extract only the CellTypist columns that are present
celltypist_annotations <- annotation_coldata[, present_celltypist_cols, drop = FALSE]

cat("Adding CellTypist annotation columns:\n")
print(colnames(celltypist_annotations))

# Add CellTypist annotations to original COVID data
cat("Adding CellTypist annotations to original COVID data...\n")
original_coldata <- colData(covid_data_original)
combined_coldata <- cbind(original_coldata, celltypist_annotations)

# Update the original SCE object with annotations
colData(covid_data_original) <- combined_coldata

cat("Final colData columns:\n")
print(colnames(colData(covid_data_original)))

# Summary of annotations (conditional on what's available)
cat("\n=== CellTypist Annotation Summary ===\n")

# Always check for predicted_labels (this should always be present)
if ("predicted_labels" %in% colnames(colData(covid_data_original))) {
  cat("Predicted cell types (top 10):\n")
  pred_table <- sort(table(colData(covid_data_original)$predicted_labels), decreasing = TRUE)
  print(head(pred_table, 10))
} else {
  cat("⚠️  Warning: 'predicted_labels' column not found!\n")
}

# Check for majority_voting (may or may not be present)
if ("majority_voting" %in% colnames(colData(covid_data_original))) {
  cat("\nMajority voting cell types (top 10):\n")
  maj_table <- sort(table(colData(covid_data_original)$majority_voting), decreasing = TRUE)
  print(head(maj_table, 10))
  
  # Compare predicted_labels vs majority_voting if both are present
  if ("predicted_labels" %in% colnames(colData(covid_data_original))) {
    pred_labels <- colData(covid_data_original)$predicted_labels
    maj_labels <- colData(covid_data_original)$majority_voting
    agreement <- sum(pred_labels == maj_labels, na.rm = TRUE) / length(pred_labels)
    cat("Agreement between predicted_labels and majority_voting:", 
        round(100 * agreement, 1), "%\n")
  }
} else {
  cat("\nMajority voting: Not performed (majority_voting=False in Python script)\n")
}

# Check for confidence scores
if ("conf_score" %in% colnames(colData(covid_data_original))) {
  conf_scores <- colData(covid_data_original)$conf_score
  cat("\nConfidence score summary:\n")
  print(summary(conf_scores))
  
  low_conf_count <- sum(conf_scores < 0.5, na.rm = TRUE)
  low_conf_pct <- round(100 * low_conf_count / length(conf_scores), 1)
  cat("Low confidence cells (<0.5):", low_conf_count, 
      paste0("(", low_conf_pct, "%)\n"))
  
  high_conf_count <- sum(conf_scores >= 0.8, na.rm = TRUE)
  high_conf_pct <- round(100 * high_conf_count / length(conf_scores), 1)
  cat("High confidence cells (≥0.8):", high_conf_count, 
      paste0("(", high_conf_pct, "%)\n"))
} else {
  cat("\nConfidence scores: Not available\n")
}

# SAVE BACK TO ORIGINAL FILE (overwriting it)
cat("Saving annotated data back to original file...\n")
saveRDS(covid_data_original, ORIGINAL_COVID_FILE)

cat("\n=== Integration Complete ===\n")
cat("CellTypist annotations added to:", ORIGINAL_COVID_FILE, "\n")
cat("Final object dimensions:", dim(covid_data_original), "\n")
cat("Total colData columns:", ncol(colData(covid_data_original)), "\n")

# Summary of what was added
cat("\nCellTypist annotations successfully integrated:\n")
for(col in present_celltypist_cols) {
  cat("✓ ", col, "\n", sep = "")
}