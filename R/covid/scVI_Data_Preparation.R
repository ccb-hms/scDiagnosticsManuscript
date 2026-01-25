# --------------------------------------------------
# COVID-19 PBMC - scVI/scArches Data Preparation
# --------------------------------------------------

library(SingleCellExperiment)
library(zellkonverter)
library(scran)  

# Configuration 
REFERENCE_SCE_FILE <- "data/covid/normal_data_sce.rds"
QUERY_SCE_FILE <- "data/covid/covid_data_sce.rds"
ANNDATA_REF_FILE <- "data/covid/scvi_reference_data.h5ad"
ANNDATA_QUERY_FILE <- "data/covid/scvi_query_data.h5ad"
N_HVGS <- 2000  

message("Preparing data for scVI/scArches with HVG filtering")

# Load SCE objects
message("Loading SCE objects...")
ref_sce <- readRDS(REFERENCE_SCE_FILE)
query_sce <- readRDS(QUERY_SCE_FILE)

message("Reference dimensions: ", nrow(ref_sce), " genes x ", ncol(ref_sce), " cells")
message("Query dimensions: ", nrow(query_sce), " genes x ", ncol(query_sce), " cells")

# _______________________________
# Compute HVGs for both datasets
# _______________________________

message("\nComputing highly variable genes...")

# Use logcounts if available, otherwise compute from counts
if (!"logcounts" %in% assayNames(ref_sce)) {
    library(scuttle)
    ref_sce <- logNormCounts(ref_sce)
    query_sce <- logNormCounts(query_sce)
}

# Compute HVGs for reference
message("  Computing HVGs for reference...")
dec_ref <- modelGeneVar(ref_sce)
hvg_ref <- getTopHVGs(dec_ref, n = N_HVGS)

# Compute HVGs for query
message("  Computing HVGs for query...")
dec_query <- modelGeneVar(query_sce)
hvg_query <- getTopHVGs(dec_query, n = N_HVGS)

# Take union of HVGs
hvg_union <- union(hvg_ref, hvg_query)
message("  Reference HVGs: ", length(hvg_ref))
message("  Query HVGs: ", length(hvg_query))
message("  Union HVGs: ", length(hvg_union))

# Subset both datasets to union of HVGs
ref_sce <- ref_sce[hvg_union, ]
query_sce <- query_sce[hvg_union, ]

message("\nAfter HVG filtering:")
message("  Reference: ", nrow(ref_sce), " genes x ", ncol(ref_sce), " cells")
message("  Query: ", nrow(query_sce), " genes x ", ncol(query_sce), " cells")

# __________________________
# Prepare for scVI/scArches
# __________________________

prepare_for_scvi <- function(sce, is_reference = TRUE) {
    if (!"counts" %in% assayNames(sce)) {
        stop("A 'counts' assay is required for scVI but was not found.")
    }

    essential_cols <- "sample_id"
    if (is_reference) {
        essential_cols <- c(essential_cols, "author_cell_type")
    }

    missing_cols <- setdiff(essential_cols, colnames(colData(sce)))
    if (length(missing_cols) > 0) {
        stop("Missing required colData columns: ", paste(missing_cols, collapse = ", "))
    }
    colData(sce) <- colData(sce)[, essential_cols, drop = FALSE]

    sce <- SingleCellExperiment(
        assays = list(counts = assay(sce, "counts")),
        colData = colData(sce),
        rowData = rowData(sce)
    )
    adata <- SCE2AnnData(sce, X_name = "counts")

    return(adata)
}

# Prepare and save
message("\nPreparing reference data...")
adata_ref <- prepare_for_scvi(ref_sce, is_reference = TRUE)

message("Preparing query data...")
adata_query <- prepare_for_scvi(query_sce, is_reference = FALSE)

message("\nSaving AnnData files for Python...")
adata_ref$write_h5ad(ANNDATA_REF_FILE)
adata_query$write_h5ad(ANNDATA_QUERY_FILE)

message("✓ Data preparation for scVI/scArches complete.")
message("Files created:")
message("- ", ANNDATA_REF_FILE)
message("- ", ANNDATA_QUERY_FILE)

rm(list = ls())
gc()