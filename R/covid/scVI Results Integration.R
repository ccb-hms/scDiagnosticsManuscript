# ----------------------------------------------
# COVID-19 - scVI/scArches Results Integration
# ----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(dplyr)

# Load source file
source("R/auxiliary/addMergedCellTypes.R")

# ----------------------------------------------

# Configuration
ORIGINAL_COVID_FILE <- "data/covid/covid_data_sce.rds"
ORIGINAL_NORMAL_FILE <- "data/covid/normal_data_sce.rds"
QUERY_PRED_CSV_FILE <- "data/covid/scvi_predictions.csv"
REF_UMAP_CSV_FILE <- "data/covid/scvi_reference_umap.csv"

# Main Script
message("--- Integrating scVI/scArches Annotation Results ---")

# Load Data
message("Loading original R objects...")
covid_data <- readRDS(ORIGINAL_COVID_FILE)
normal_data <- readRDS(ORIGINAL_NORMAL_FILE)

message("Loading predictions and UMAP coordinates from Python CSV outputs...")
query_df <- read.csv(QUERY_PRED_CSV_FILE, row.names = 1, check.names = FALSE)
ref_umap_df <- read.csv(REF_UMAP_CSV_FILE, row.names = 1, check.names = FALSE)

# _________________________________________
# Process and Integrate COVID (Query) Data 
# _________________________________________

message("\n--- Processing Query (COVID) Data ---")

# Align Query Data
message("Validating and aligning query cell barcodes...")
original_query_barcodes <- colnames(covid_data)
# Reorder the python output to match the original R object perfectly
query_df <- query_df[original_query_barcodes, , drop = FALSE]

if(any(is.na(rownames(query_df)))) {
    stop("Error: Query alignment failed. Not all barcodes from the original SCE were found.")
}
message("✓ Query cell barcodes aligned.")

# Integrate Query Annotations and UMAP
message("Integrating query annotations and UMAP coordinates...")
covid_data$scvi_prediction <- query_df$scvi_prediction
reducedDim(covid_data, "UMAP_scVI") <- as.matrix(query_df[, c("UMAP_scVI_1", "UMAP_scVI_2")])

# Create Merged Annotation Column for Query
covid_data <- addMergedCellTypes(
    sce_object = covid_data,
    input_col_name = "scvi_prediction"
)
message("✓ Query data integration complete.")

# ______________________________________________
# Process and Integrate Normal (Reference) Data 
# ______________________________________________

message("\n--- Processing Reference (Normal) Data ---")

# Align Reference Data
message("Validating and aligning reference cell barcodes...")
original_ref_barcodes <- colnames(normal_data)
ref_umap_df <- ref_umap_df[original_ref_barcodes, , drop = FALSE]

if(any(is.na(rownames(ref_umap_df)))) {
    stop("Error: Reference alignment failed. Not all barcodes from the original SCE were found.")
}
message("✓ Reference cell barcodes aligned.")

# Integrate Reference UMAP
message("Integrating reference UMAP coordinates...")
reducedDim(normal_data, "UMAP_scVI") <- as.matrix(ref_umap_df)
message("✓ Reference data integration complete.")

# _______________________
# Final Save and Cleanup
# _______________________

message("\n--- Saving Final Annotated Objects ---")

# Save Both Final Objects
saveRDS(covid_data, ORIGINAL_COVID_FILE)
saveRDS(normal_data, ORIGINAL_NORMAL_FILE)

message("\n✓✓✓ scVI Integration Complete for both datasets! ✓✓✓")
message("Updated file: ", ORIGINAL_COVID_FILE)
message("Updated file: ", ORIGINAL_NORMAL_FILE)

rm(list=ls())
gc()