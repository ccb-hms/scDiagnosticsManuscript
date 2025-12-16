# ---------------------------------------------------------------
# MERFISH Mouse Colon IBD - Data Preparation (Healthy vs. DSS9)
# ---------------------------------------------------------------

# Load libraries
library(MerfishData)
library(SpatialExperiment)
library(scran)
library(scater)
library(BiocParallel)
library(Matrix) 

# Source auxiliary functions
source("R/auxiliary/addReferencePCA.R")
source("R/auxiliary/addMergedCellTypes.R")

# ---------------------------------------------------------------

# ______________________________________
# Data Loading and Condition Splitting
# ______________________________________

# Load MERFISH mouse colon IBD data
spe <- MouseColonIbdCadinu2024()

# Filter for the slice of interest
spe_s2 <- spe[, spe$slice_id == "2"]
message("Filtered for slice 2, containing ", ncol(spe_s2), " cells.")

# Split into two separate objects
spe_healthy_full <- spe_s2[, spe_s2$sample_type == "Healthy"]
spe_dss9_full <- spe_s2[, spe_s2$sample_type == "DSS9"]

# Select a single mouse ID for each dataset to avoid batch effects
sce_healthy <- spe_healthy_full[, spe_healthy_full$mouse_id == "082421_D0_m6"]
sce_dss9    <- spe_dss9_full[, spe_dss9_full$mouse_id == "062221_D9_m3"]

message("Filtered for single mouse per condition. Full cell counts: ",
        "Healthy = ", ncol(sce_healthy), ", ",
        "DSS9 = ", ncol(sce_dss9))


# _________________________
# Per-Cell Quality Control
# _________________________

# Identify mitochondrial genes
is_mito <- grepl("^mt-", rownames(sce_healthy), ignore.case = TRUE)
message("Found ", sum(is_mito), " mitochondrial genes in the panel.")

# Calculate QC metrics for both objects
sce_healthy <- addPerCellQC(sce_healthy, subsets = list(Mito = is_mito))
sce_dss9    <- addPerCellQC(sce_dss9, subsets = list(Mito = is_mito))

# Identify outliers using adaptive thresholds
qc_lib_size_h   <- isOutlier(sce_healthy$sum, log = TRUE, type = "lower")
qc_n_features_h <- isOutlier(sce_healthy$detected, log = TRUE, type = "lower")
qc_mito_h       <- isOutlier(sce_healthy$subsets_Mito_percent, type = "higher")
discard_h <- qc_lib_size_h | qc_n_features_h | qc_mito_h
sce_healthy <- sce_healthy[, !discard_h]

qc_lib_size_d9   <- isOutlier(sce_dss9$sum, log = TRUE, type = "lower")
qc_n_features_d9 <- isOutlier(sce_dss9$detected, log = TRUE, type = "lower")
qc_mito_d9       <- isOutlier(sce_dss9$subsets_Mito_percent, type = "higher")
discard_d9 <- qc_lib_size_d9 | qc_n_features_d9 | qc_mito_d9
sce_dss9 <- sce_dss9[, !discard_d9]

message("After QC, cells remaining: ",
        "Healthy = ", ncol(sce_healthy), " (removed ", sum(discard_h), "), ",
        "DSS9 = ", ncol(sce_dss9), " (removed ", sum(discard_d9), ")")


# _____________________________________________
# Biological Filtering (Healthy Reference Only)
# _____________________________________________

# Goal: Remove inflammation-associated cells (IAE, IAF, IASMC) from the Healthy dataset
# so the reference represents a truly "healthy" state.

inflamed_patterns <- "^(IAE|IAF|IASMC)"
is_inflamed <- grepl(inflamed_patterns, sce_healthy$tier2)

if (sum(is_inflamed) > 0) {
    message(sprintf("Removing %d inflammation-associated cells (IAE, IAF, IASMC) from Healthy dataset.", sum(is_inflamed)))
    sce_healthy <- sce_healthy[, !is_inflamed]
} else {
    message("No inflammation-associated cells found in Healthy dataset.")
}

message("Healthy reference count after biological filtering: ", ncol(sce_healthy))


# _______________
# Normalization
# _______________

# Re-calculate logcounts on the filtered data to ensure correct library size factors
sce_healthy <- logNormCounts(sce_healthy, assay.type = "counts")
sce_dss9    <- logNormCounts(sce_dss9, assay.type = "counts")
message("Re-calculated logcounts on filtered data.")


# _________________________
# Dimensionality Reduction
# _________________________

# Use the custom function to compute HVGs and run PCA on the reference object.
sce_healthy <- addReferencePCA(
    ref_sce = sce_healthy,
    query_sce = sce_dss9,
    ref_name = "Healthy_vs_DSS9"
)

# Apply PCA to the query object using the same diagnostic HVGs for comparability
message("\nApplying the same PCA space to the query (DSS9) object...")
diagnostic_hvgs <- metadata(sce_healthy)$diagnostic_hvgs
sce_dss9 <- runPCA(sce_dss9, subset_row = diagnostic_hvgs, ncomponents = 50)
message("PCA computed for both objects.")


# _____________________________
# Cell Type Annotation Merging
# _____________________________

message("Merging fine-grained cell type annotations into broader categories...")

# Use function to add the merged column
sce_healthy <- addMergedCellTypes(sce_object = sce_healthy, input_col_name = "tier2",
                                  dataset = "MERFISH")
sce_dss9 <- addMergedCellTypes(sce_object = sce_dss9, input_col_name = "tier2",
                               dataset = "MERFISH")

message("Cell type merging complete for both datasets.")
message("Merged cell type distribution in reference (normal):")
print(sort(table(sce_healthy$tier2_merged), decreasing = TRUE))


# _____________________________
# Final Formatting and Saving
# _____________________________

# Define a function to apply final formatting steps
formatAndFinalizeSCE <- function(sce_object, cols_to_keep) {

    # 1. Subset colData to keep only necessary columns
    # Ensure 'sample_id' or 'mouse_id' is kept if needed for scVI batch correction
    # Adding "mouse_id" to the keep list just in case.
    cols_to_keep <- c(cols_to_keep, "mouse_id", "sample_id")
    cols_to_keep <- intersect(cols_to_keep, colnames(colData(sce_object)))
    colData(sce_object) <- colData(sce_object)[, cols_to_keep, drop = FALSE]
    
    # 2. Handle Assays (Crucial for scVI)
    # scVI requires raw counts. We must NOT delete them.
    # Convert both counts and logcounts to in-memory sparse matrices.
    
    if ("counts" %in% assayNames(sce_object)) {
        assay(sce_object, "counts") <- as(assay(sce_object, "counts"), "CsparseMatrix")
    }
    
    if ("logcounts" %in% assayNames(sce_object)) {
        assay(sce_object, "logcounts") <- as(assay(sce_object, "logcounts"), "CsparseMatrix")
    }
    
    # 3. Clean up
    # This removes the nested 'Blank' experiment that holds on-disk data
    altExps(sce_object) <- NULL
  
    # 4. Add cell names
    colnames(sce_object) <- paste0("cell", 1:ncol(sce_object))
    
    return(sce_object)
}

# Apply final formatting to both objects
# Note: Adding "sample_id" if it exists, as the scVI script looks for it
if(!"sample_id" %in% colnames(colData(sce_healthy))) sce_healthy$sample_id <- "Healthy_Ref"
if(!"sample_id" %in% colnames(colData(sce_dss9))) sce_dss9$sample_id <- "DSS9_Query"

sce_healthy_final <- formatAndFinalizeSCE(sce_healthy, cols_to_keep = c("tier2", "tier2_merged"))
sce_dss9_final    <- formatAndFinalizeSCE(sce_dss9, cols_to_keep = c("tier2", "tier2_merged"))

# Create directory if it doesn't exist
if (!dir.exists("data/merfish")) {
  dir.create("data/merfish", recursive = TRUE)
}

# Save the final, fully in-memory objects
saveRDS(sce_healthy_final, file = "data/merfish/healthy_data.rds")
saveRDS(sce_dss9_final, file = "data/merfish/dss9_data.rds")
message("Saved final, formatted objects to 'data/merfish/'")


# _______________
# Final Summary
# _______________

message("\n--- Pipeline Complete! ---")
message("Final 'Healthy' dataset contains ", ncol(sce_healthy_final), " cells.")
message("Final 'DSS9' dataset contains ", ncol(sce_dss9_final), " cells.")

# Clean up memory
rm(list = ls())
gc()