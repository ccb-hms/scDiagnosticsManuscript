# ---------------------------------------------------------
# MERFISH Mouse Colon IBD - scVI/scArches Data Preparation
# ---------------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(zellkonverter)

# ---------------------------------------------------------

# Configuration 
REFERENCE_SCE_FILE <- "data/merfish/healthy_data.rds"
QUERY_SCE_FILE <- "data/merfish/dss9_data.rds"
ANNDATA_REF_FILE <- "data/merfish/scvi_reference_data.h5ad"
ANNDATA_QUERY_FILE <- "data/merfish/scvi_query_data.h5ad"

# Main Script 
message("Preparing MERFISH data for scVI/scArches ")

# Load SCE objects
message("Loading SCE objects...")
ref_sce <- readRDS(REFERENCE_SCE_FILE)
query_sce <- readRDS(QUERY_SCE_FILE)

message("Reference (Healthy) dimensions: ", nrow(ref_sce), " genes x ", ncol(ref_sce), " cells")
message("Query (DSS9) dimensions: ", nrow(query_sce), " genes x ", ncol(query_sce), " cells")

# Function to clean and prepare an SCE for scVI
prepare_for_scvi <- function(sce, is_reference = TRUE) {
    # 1. scVI requires raw counts. Verify the 'counts' assay exists.
    if (!"counts" %in% assayNames(sce)) {
        stop("A 'counts' assay is required for scVI but was not found.")
    }

    # 2. Keep only essential colData columns.
    # scVI needs 'sample_id' for batch correction.
    # The reference also needs 'tier2' for annotation training.
    
    # NOTE: Ensure 'sample_id' exists in your object. If not, and you don't have batches, 
    # you may need to create a dummy column: sce$sample_id <- "batch1"
    essential_cols <- "sample_id" 
    
    if (is_reference) {
        essential_cols <- c(essential_cols, "tier2")
    }

    missing_cols <- setdiff(essential_cols, colnames(colData(sce)))
    if (length(missing_cols) > 0) {
        stop("Missing required colData columns: ", paste(missing_cols, collapse = ", "))
    }
    colData(sce) <- colData(sce)[, essential_cols, drop = FALSE]

    # 3. Create the AnnData object using raw counts and clean up other slots.
    sce <- SingleCellExperiment(
        assays = list(counts = assay(sce, "counts")),
        colData = colData(sce),
        rowData = rowData(sce) # Keep gene names
    )
    adata <- SCE2AnnData(sce, X_name = "counts")

    return(adata)
}

# Prepare reference and query data
message("\nPreparing reference data...")
adata_ref <- prepare_for_scvi(ref_sce, is_reference = TRUE)

message("Preparing query data...")
adata_query <- prepare_for_scvi(query_sce, is_reference = FALSE)

# Save AnnData objects
message("\nSaving AnnData files for Python...")
adata_ref$write_h5ad(ANNDATA_REF_FILE)
adata_query$write_h5ad(ANNDATA_QUERY_FILE)

message("✓ Data preparation for scVI/scArches complete.")
message("Files created:")
message("- ", ANNDATA_REF_FILE)
message("- ", ANNDATA_QUERY_FILE)

# Clean up memory
rm(list = ls())
gc()