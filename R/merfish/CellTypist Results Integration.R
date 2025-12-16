# ----------------------------------------------------------
# MERFISH Mouse Colon IBD - CellTypist Results Integration
# ----------------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)

# Source auxiliary functions
source("R/auxiliary/environmentSetupCellTypist.R")
source("R/auxiliary/addMergedCellTypes.R")

# ----------------------------------------------------------

# ________________________
# Setup and Data Loading
# ________________________

message("--- Initializing CellTypist Integration (MERFISH) ---")

# Ensure Python environment is set up
message("Ensuring Python environment is available...")
environmentSetupCellTypist()

# Define file paths
# Note: Integrating back into DSS9 query data
ORIGINAL_QUERY_FILE <- "data/merfish/dss9_data.rds"
ANNOTATED_QUERY_FILE <- "data/merfish/annotated_query_data.h5ad"

# Load data
message("\nLoading original DSS9 SCE object...")
dss9_data_original <- readRDS(ORIGINAL_QUERY_FILE)

message("Loading CellTypist annotated results from Python...")
adata_annotated <- reticulate::import("scanpy")$read_h5ad(ANNOTATED_QUERY_FILE)
sce_annotations_only <- AnnData2SCE(adata_annotated)

message("Original DSS9 data dimensions: ", ncol(dss9_data_original), " cells x ", nrow(dss9_data_original), " genes")
message("Annotation data dimensions: ", ncol(sce_annotations_only), " cells x ", nrow(sce_annotations_only), " genes")

# _______________________________
# Data Validation and Alignment
# _______________________________

message("\n--- Validating Cell Alignment ---")

cell_names_original <- colnames(dss9_data_original)
cell_names_annotated <- colnames(sce_annotations_only)

if (!identical(cell_names_original, cell_names_annotated)) {
  warning("Cell order doesn't match exactly, attempting to align...")
  
  common_cells <- intersect(cell_names_original, cell_names_annotated)
  message("Common cells found: ", length(common_cells))
  
  if(length(common_cells) == 0) {
    stop("No common cells found between original and annotated data!")
  }
  
  # Align both objects to the common set of cells in the correct order
  dss9_data_original <- dss9_data_original[, common_cells]
  sce_annotations_only <- sce_annotations_only[, common_cells]
  
  message("Alignment complete. Both datasets now have ", ncol(dss9_data_original), " cells.")
} else {
    message("✓ Cell names and order match perfectly.")
}

# ____________________________________
# Annotation Integration and Renaming
# ____________________________________

message("\n--- Integrating and Renaming CellTypist Annotations ---")

# 1. Identify which columns are present from the Python script
annotation_coldata <- colData(sce_annotations_only)
possible_python_cols <- c("predicted_labels", "majority_voting", "conf_score")
present_python_cols <- intersect(possible_python_cols, colnames(annotation_coldata))

if(length(present_python_cols) == 0) {
  stop("No CellTypist annotation columns found! Expected at least 'predicted_labels'.")
}
message("Found raw Python columns: ", paste(present_python_cols, collapse = ", "))

# 2. Extract only these columns
celltypist_annotations <- annotation_coldata[, present_python_cols, drop = FALSE]

# 3. Rename the columns to add the 'celltypist_' prefix for clarity
name_map <- c(
    "predicted_labels"  = "celltypist_predicted_labels",
    "majority_voting"   = "celltypist_majority_voting",
    "conf_score"        = "celltypist_conf_score"
)
current_names <- colnames(celltypist_annotations)
colnames(celltypist_annotations) <- name_map[current_names]
message("Renamed columns to: ", paste(colnames(celltypist_annotations), collapse = ", "))

# 4. Add the renamed columns to the original SCE object
colData(dss9_data_original) <- cbind(colData(dss9_data_original), celltypist_annotations)
message("✓ Prefixed annotation columns successfully added.")

# __________________________________
# Create Merged Annotation Columns
# __________________________________

message("\n--- Creating Merged Versions of Annotations ---")

# Use the new, prefixed column names as input for the merging function
if ("celltypist_predicted_labels" %in% colnames(colData(dss9_data_original))) {
    dss9_data_original <- addMergedCellTypes(
        sce_object = dss9_data_original,
        input_col_name = "celltypist_predicted_labels", 
        dataset = "MERFISH"
    )
}

# Defensively merge majority_voting labels only if they exist
if ("celltypist_majority_voting" %in% colnames(colData(dss9_data_original))) {
    dss9_data_original <- addMergedCellTypes(
        sce_object = dss9_data_original,
        input_col_name = "celltypist_majority_voting", 
        dataset = "MERFISH"
    )
}
message("✓ Merged cell type columns created.")

# __________________________
# Final Summary and Saving
# __________________________

message("\n--- Final Summary and Validation ---")

# Check for the primary merged predictions
if ("celltypist_predicted_labels_merged" %in% colnames(colData(dss9_data_original))) {
  cat("\nMerged predicted cell types (top 10):\n")
  print(head(sort(table(dss9_data_original$celltypist_predicted_labels_merged), decreasing = TRUE), 10))
} else {
  cat("\nWarning: Merged predicted labels column not found!\n")
}

# Check for confidence scores
if ("celltypist_conf_score" %in% colnames(colData(dss9_data_original))) {
  conf_scores <- dss9_data_original$celltypist_conf_score
  cat("\nConfidence score summary:\n")
  print(summary(conf_scores))
} else {
  cat("\nConfidence scores: Not available.\n")
}

# Save the final, fully annotated object
message("\nSaving fully annotated data back to original file...")
saveRDS(dss9_data_original, ORIGINAL_QUERY_FILE)

message("\n✓✓✓ Integration Complete! ✓✓✓")
message("CellTypist annotations have been successfully added to: ", ORIGINAL_QUERY_FILE)

# Final confirmation of what was added
final_celltypist_cols <- grep("celltypist", colnames(colData(dss9_data_original)), value = TRUE)
message("\nFinal CellTypist-related columns in object:")
for(col in final_celltypist_cols) {
  message("- ", col)
}