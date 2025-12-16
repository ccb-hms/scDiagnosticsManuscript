# -------------------------------------------------------------
# MERFISH Mouse Colon IBD - scVI/scArches Results Integration
# -------------------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(dplyr)

# Load source file
source("R/auxiliary/addMergedCellTypes.R")

# -------------------------------------------------------------

# Configuration
ORIGINAL_QUERY_FILE <- "data/merfish/dss9_data.rds"
ORIGINAL_REF_FILE <- "data/merfish/healthy_data.rds"
QUERY_PRED_CSV_FILE <- "data/merfish/scvi_predictions.csv"
REF_UMAP_CSV_FILE <- "data/merfish/scvi_reference_umap.csv"

# Main Script
message("--- Integrating scVI/scArches Annotation Results (MERFISH) ---")

# Load Data
message("Loading original R objects...")
dss9_data <- readRDS(ORIGINAL_QUERY_FILE)
healthy_data <- readRDS(ORIGINAL_REF_FILE)

message("Loading predictions and UMAP coordinates from Python CSV outputs...")
query_df <- read.csv(QUERY_PRED_CSV_FILE, row.names = 1, check.names = FALSE)
ref_umap_df <- read.csv(REF_UMAP_CSV_FILE, row.names = 1, check.names = FALSE)

# _________________________________________
# Process and Integrate DSS9 (Query) Data 
# _________________________________________

message("\n--- Processing Query (DSS9) Data ---")

# Align Query Data
message("Validating and aligning query cell barcodes...")
original_query_barcodes <- colnames(dss9_data)
# Reorder the python output to match the original R object perfectly
query_df <- query_df[original_query_barcodes, , drop = FALSE]

if(any(is.na(rownames(query_df)))) {
    stop("Error: Query alignment failed. Not all barcodes from the original SCE were found.")
}
message("✓ Query cell barcodes aligned.")

# Integrate Query Annotations and UMAP
message("Integrating query annotations and UMAP coordinates...")
dss9_data$scvi_prediction <- query_df$scvi_prediction
reducedDim(dss9_data, "UMAP_scVI") <- as.matrix(query_df[, c("UMAP_scVI_1", "UMAP_scVI_2")])

# Create Merged Annotation Column for Query
dss9_data <- addMergedCellTypes(
    sce_object = dss9_data,
    input_col_name = "scvi_prediction", 
    dataset = "MERFISH"
)
message("✓ Query data integration complete.")

# ______________________________________________
# Process and Integrate Healthy (Reference) Data 
# ______________________________________________

message("\n--- Processing Reference (Healthy) Data ---")

# Align Reference Data
message("Validating and aligning reference cell barcodes...")
original_ref_barcodes <- colnames(healthy_data)
ref_umap_df <- ref_umap_df[original_ref_barcodes, , drop = FALSE]

if(any(is.na(rownames(ref_umap_df)))) {
    stop("Error: Reference alignment failed. Not all barcodes from the original SCE were found.")
}
message("✓ Reference cell barcodes aligned.")

# Integrate Reference UMAP
message("Integrating reference UMAP coordinates...")
reducedDim(healthy_data, "UMAP_scVI") <- as.matrix(ref_umap_df)
message("✓ Reference data integration complete.")

# _______________________
# Final Save and Cleanup
# _______________________

message("\n--- Saving Final Annotated Objects ---")

# Save Both Final Objects
saveRDS(dss9_data, ORIGINAL_QUERY_FILE)
saveRDS(healthy_data, ORIGINAL_REF_FILE)

message("\n✓✓✓ scVI Integration Complete for both datasets! ✓✓✓")
message("Updated file: ", ORIGINAL_QUERY_FILE)
message("Updated file: ", ORIGINAL_REF_FILE)

rm(list=ls())
gc()