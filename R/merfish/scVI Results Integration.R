# -----------------------------------------------
# MERFISH - scVI/scArches Results Integration
# -----------------------------------------------

library(SingleCellExperiment)
library(dplyr)

source("R/auxiliary/addMergedCellTypes.R")

# Configuration (MERFISH-specific)
ORIGINAL_DAY0_FILE <- "data/merfish/healthy_data.rds"
ORIGINAL_DAY9_FILE <- "data/merfish/dss9_data.rds"
QUERY_PRED_CSV_FILE <- "data/merfish/scvi_predictions.csv"
REF_UMAP_CSV_FILE <- "data/merfish/scvi_reference_umap.csv"

message("--- Integrating scVI/scArches Annotation Results (MERFISH) ---")

# Load Data
message("Loading original R objects...")
merfish_day0 <- readRDS(ORIGINAL_DAY0_FILE)
merfish_day9 <- readRDS(ORIGINAL_DAY9_FILE)

message("Loading predictions, confidence scores, and UMAP coordinates from Python CSV outputs...")
query_df <- read.csv(QUERY_PRED_CSV_FILE, row.names = 1, check.names = FALSE)
ref_umap_df <- read.csv(REF_UMAP_CSV_FILE, row.names = 1, check.names = FALSE)

# _________________________________________
# Process and Integrate Day 9 (Query) Data 
# _________________________________________

message("\n--- Processing Query (Day 9) Data ---")

message("Validating and aligning query cell barcodes...")
original_query_barcodes <- colnames(merfish_day9)
query_df <- query_df[original_query_barcodes, , drop = FALSE]

if(any(is.na(rownames(query_df)))) {
    stop("Error: Query alignment failed.")
}
message("✓ Query cell barcodes aligned.")

message("Integrating query annotations, confidence scores, and UMAP coordinates...")
merfish_day9$scvi_prediction <- query_df$scvi_prediction
merfish_day9$scvi_confidence <- query_df$scvi_confidence
reducedDim(merfish_day9, "UMAP_scVI") <- as.matrix(query_df[, c("UMAP_scVI_1", "UMAP_scVI_2")])

merfish_day9 <- addMergedCellTypes(
    sce_object = merfish_day9,
    input_col_name = "scvi_prediction", 
    dataset = "MERFISH"
)
message("✓ Query data integration complete.")

# ______________________________________________
# Process and Integrate Day 0 (Reference) Data 
# ______________________________________________

message("\n--- Processing Reference (Day 0) Data ---")

message("Validating and aligning reference cell barcodes...")
original_ref_barcodes <- colnames(merfish_day0)
ref_umap_df <- ref_umap_df[original_ref_barcodes, , drop = FALSE]

if(any(is.na(rownames(ref_umap_df)))) {
    stop("Error: Reference alignment failed.")
}
message("✓ Reference cell barcodes aligned.")

message("Integrating reference UMAP coordinates...")
reducedDim(merfish_day0, "UMAP_scVI") <- as.matrix(ref_umap_df)
message("✓ Reference data integration complete.")

# _______________________
# Final Save and Cleanup
# _______________________

message("\n--- Saving Final Annotated Objects ---")

saveRDS(merfish_day9, ORIGINAL_DAY9_FILE)
saveRDS(merfish_day0, ORIGINAL_DAY0_FILE)

message("\n✓✓✓ scVI/scArches Integration Complete for MERFISH! ✓✓✓")
message("Updated file: ", ORIGINAL_DAY9_FILE)
message("Updated file: ", ORIGINAL_DAY0_FILE)

rm(list=ls())
gc()