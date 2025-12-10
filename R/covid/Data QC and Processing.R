# ----------------------------------------------------------------------------
# COVID-19 PBMC - Data QC, Processing, and Splitting
# ----------------------------------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scater)
library(scran)
library(biomaRt)
library(Matrix)

# Source auxiliary functions
source("R/auxiliary/convertEnsemblToSymbols.R")
source("R/auxiliary/addReferencePCA.R")
source("R/auxiliary/addMergedCellTypes.R")

# ----------------------------------------------------------------------------

# _____________
# Data Loading
# _____________

message("Loading cleaned SCE object from 'data/covid/covid_data_clean.rds'...")
sce <- readRDS("data/covid/covid_data_clean.rds")

message(sprintf("Loaded data with %d genes x %d cells.", nrow(sce), ncol(sce)))

# Setting the seed
set.seed(0)

# ____________________________________________
# Assay Renaming and Pseudo-Count Generation
# ____________________________________________

message("Original data is log-normalized. Renaming primary assay to 'logcounts'.")

# This block makes the script robust, correctly identifying the log-transformed
# data whether it's named 'X' or was mislabeled 'counts' in the prior step.
if ("X" %in% assayNames(sce)) {
    assayNames(sce)[assayNames(sce) == "X"] <- "logcounts"
} else if ("counts" %in% assayNames(sce) && !"logcounts" %in% assayNames(sce) && max(assay(sce, "counts")) < 50) {
    message("Found assay misnamed 'counts'; renaming to 'logcounts'.")
    assayNames(sce)[assayNames(sce) == "counts"] <- "logcounts"
}

message("Generating a pseudo-counts assay by reversing the log-transformation...")
# Reverse the log-normalization (assumed to be log2(counts + 1))
logcounts_matrix <- assay(sce, "logcounts")
pseudo_counts_matrix <- as(2^logcounts_matrix - 1, "CsparseMatrix")

# Ensure no negative values from floating point inaccuracies
pseudo_counts_matrix[pseudo_counts_matrix < 0] <- 0

# Add the new 'counts' assay to the object
assay(sce, "counts") <- pseudo_counts_matrix

# Reorder assays for convention (counts first)
assay_order <- intersect(c("counts", "logcounts"), assayNames(sce))
assays(sce) <- assays(sce)[assay_order]

message("Successfully generated 'counts' assay. Available assays:")
print(assayNames(sce))

# __________________________
# Quality Control Filtering
# __________________________

message("Applying QC filters to cells and genes...")

# 1. Filter cells with low total UMI counts (from original colData)
keep_cells <- sce$total_counts >= 1000
message(sprintf("Filtering cells: Keeping %d of %d cells with >= 1000 total counts.", sum(keep_cells), ncol(sce)))

# 2. Filter lowly expressed genes using the newly generated 'counts' matrix
keep_genes <- rowSums(assay(sce, "counts")[, keep_cells] > 0) >= 10
message(sprintf("Filtering genes: Keeping %d of %d genes expressed in >= 10 cells.", sum(keep_genes), nrow(sce)))

# Apply filters
sce_filtered <- sce[keep_genes, keep_cells]
message(sprintf("Dimensions after QC: %d genes x %d cells.", nrow(sce_filtered), ncol(sce_filtered)))

# __________________________________
# Ensembl to Gene Symbol Conversion
# __________________________________

message("Converting Ensembl gene IDs to HGNC symbols...")
# This function uses biomaRt and handles duplicates/missing symbols.
sce_symbols <- convertEnsemblToSymbols(sce_filtered)
message("Gene ID conversion complete.")

# __________________
# Dataset Splitting
# __________________

message("Splitting data into reference (normal) and query (COVID-19) sets.")

ref_sce <- sce_symbols[, sce_symbols$disease == "normal"]
query_sce <- sce_symbols[, sce_symbols$disease == "COVID-19"]

ref_sce$dataset_type <- "reference"
query_sce$dataset_type <- "query"

message(sprintf("Reference (normal) dataset: %d genes x %d cells", nrow(ref_sce), ncol(ref_sce)))
message(sprintf("Query (COVID-19) dataset: %d genes x %d cells", nrow(query_sce), ncol(query_sce)))

# Verify that both final objects have the required assays
message("Assays in reference object: ", paste(assayNames(ref_sce), collapse=", "))
message("Assays in query object: ", paste(assayNames(query_sce), collapse=", "))

# ___________________________________________
# Diagnostic-Aware HVG Selection and PCA
# ___________________________________________

message("Computing PCA using a diagnostic-aware union of HVGs...")

# This function finds HVGs in both reference and query, takes their union,
# and runs PCA on the reference using this combined feature set.
# This captures both stable and disease-specific biological signals.
ref_sce <- addReferencePCA(ref_sce, query_sce, ref_name = "Healthy Reference")

# Extract the union HVGs computed in the function
diagnostic_hvgs <- metadata(ref_sce)$diagnostic_hvgs
message(sprintf("Using %d diagnostic HVGs for PCA on both datasets.", length(diagnostic_hvgs)))

# Run PCA on the query data using the *same* set of HVGs
# modelGeneVar and runPCA both default to using 'logcounts' if available.
query_sce <- runPCA(query_sce, subset_row = diagnostic_hvgs, ncomponents = 50)

# Store HVG metadata in the query object for consistency
metadata(query_sce)$diagnostic_hvgs <- diagnostic_hvgs
metadata(query_sce)$ref_hvgs <- metadata(ref_sce)$ref_hvgs
metadata(query_sce)$query_hvgs <- metadata(ref_sce)$query_hvgs

message("PCA computed on both datasets using a shared, diagnostic-aware feature space.")

# _____________________________
# Cell Type Annotation Merging
# _____________________________

message("Merging fine-grained cell type annotations into broader categories...")

# Use function to add the merged column
ref_sce <- addMergedCellTypes(sce_object = ref_sce, input_col_name = "author_cell_type")
query_sce <- addMergedCellTypes(sce_object = query_sce, input_col_name = "author_cell_type")

message("Cell type merging complete for both datasets.")
message("Merged cell type distribution in reference (normal):")
print(sort(table(ref_sce$author_cell_type_merged), decreasing = TRUE))

# ________________________
# Save Processed Datasets
# ________________________

message("Saving final processed reference and query datasets.")

saveRDS(ref_sce, "data/covid/normal_data_sce.rds")
saveRDS(query_sce, "data/covid/covid_data_sce.rds")

message("Saved reference SCE to: data/covid/normal_data_sce.rds")
message("Saved query SCE to: data/covid/covid_data_sce.rds")

# _______________
# Memory Cleanup
# _______________

message("Cleaning up memory...")
rm(list = ls())
gc()

message("Script finished.")