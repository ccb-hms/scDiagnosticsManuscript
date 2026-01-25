# ----------------------------------------------------------------------------
# COVID-19 PBMC - Data Cleaning and Subsetting
# ----------------------------------------------------------------------------

# Load libraries
library(zellkonverter)
library(SingleCellExperiment)
library(HDF5Array)
library(Matrix)
library(dplyr)

# ----------------------------------------------------------------------------

# _____________________________
# File Paths and Initial Load
# _____________________________

input_path <- "data/covid/covid_data.h5ad"
output_path <- "data/covid/covid_data_clean.rds"

message("Loading raw h5ad file: ", input_path)
original_size_gb <- file.info(input_path)$size / (1024^3)
message(sprintf("Original file size: %.2f GB", original_size_gb))

# Load using HDF5-backed mode for memory efficiency
sce_full <- readH5AD(input_path, use_hdf5 = TRUE, reader = "R")
message(sprintf("Initial dimensions: %d genes x %d cells", nrow(sce_full), ncol(sce_full)))

# Setting the seed
set.seed(0)

# _________________
# Metadata Removal
# _________________

message("Creating a clean SCE object by removing pre-existing metadata...")

# Re-instantiate the SCE object to remove all metadata, reducedDims, altExps
# The assay is kept as a DelayedMatrix for now.
sce_cleaned <- SingleCellExperiment(
    assays = list(X = assay(sce_full, "X")),
    colData = colData(sce_full),
    rowData = rowData(sce_full)
)

message("Removed metadata, neighbors, PCA, and other pre-computed slots.")

# __________________________________
# Filtering for Relevant Conditions
# __________________________________

message("Filtering cells to retain only healthy and COVID-19 samples.")

# Filter out "respiratory system disorder" cases
keep_disease <- sce_cleaned$disease %in% c("normal", "COVID-19")
sce_filtered <- sce_cleaned[, keep_disease]

# Filter out LPS-treated samples from the healthy controls
is_healthy_control <- sce_filtered$disease == "normal" & sce_filtered$Status_on_day_collection_summary == "Healthy"
is_covid_patient <- sce_filtered$disease == "COVID-19"
keep_final <- is_healthy_control | is_covid_patient

sce_filtered <- sce_filtered[, keep_final]

message(sprintf("Retained %d cells after filtering for condition.", ncol(sce_filtered)))
message("Disease distribution after filtering:")
print(table(sce_filtered$disease))

# ________________________________
# Stratified Cell Sampling
# ________________________________

message("Performing stratified sampling by sample and condition...")

# Identify samples with a sufficient number of cells (>= 500)
sample_counts <- table(sce_filtered$sample_id)
samples_to_keep <- names(sample_counts)[sample_counts >= 500]
message(sprintf("Identified %d samples with >= 500 cells.", length(samples_to_keep)))

# Stratify samples by disease status
disease_status <- sce_filtered$disease[match(samples_to_keep, sce_filtered$sample_id)]
normal_samples <- samples_to_keep[disease_status == "normal"]
covid_samples <- samples_to_keep[disease_status == "COVID-19"]

message(sprintf("Found %d healthy control samples and %d COVID-19 samples.", length(normal_samples), length(covid_samples)))

# Sample up to 1000 cells from each normal (reference) sample
normal_indices <- unlist(sapply(normal_samples, function(id) {
    sample_indices <- which(sce_filtered$sample_id == id)
    n_to_sample <- min(length(sample_indices), 1000)
    sample(sample_indices, n_to_sample)
}))

# Sample exactly 500 cells from each COVID (query) sample
covid_indices <- unlist(sapply(covid_samples, function(id) {
    sample_indices <- which(sce_filtered$sample_id == id)
    sample(sample_indices, 500)
}))

# Combine and subset the SCE object
sce_sampled <- sce_filtered[, c(normal_indices, covid_indices)]

message(sprintf("Total cells after sampling: %d", ncol(sce_sampled)))
message(sprintf("  - Healthy cells: %d", sum(sce_sampled$disease == "normal")))
message(sprintf("  - COVID-19 cells: %d", sum(sce_sampled$disease == "COVID-19")))
message(sprintf("Sampled dimensions: %d genes x %d cells", nrow(sce_sampled), ncol(sce_sampled)))

# _________________________
# Sparse Matrix Conversion
# _________________________

message("Converting assay from DelayedMatrix to in-memory CsparseMatrix...")

# Realize the sampled data into memory
realized_matrix <- realize(assay(sce_sampled, "X"))

# Convert to the final sparse matrix format
expr_matrix <- as(realized_matrix, "CsparseMatrix")

# Create the final clean SCE object
sce_clean <- SingleCellExperiment(
    assays = list(counts = expr_matrix), # Name assay 'counts'
    colData = colData(sce_sampled),
    rowData = rowData(sce_sampled)
)

message(sprintf("Final matrix class: %s", class(assay(sce_clean, "counts"))))
message(sprintf("Final object memory size: %.2f GB", as.numeric(object.size(sce_clean)) / (1024^3)))

# ____________
# File Saving
# ____________

message("Saving cleaned and sampled SCE object to: ", output_path)
saveRDS(sce_clean, output_path, compress = TRUE)

saved_size_gb <- file.info(output_path)$size / (1024^3)
message(sprintf("Saved file size: %.2f GB", saved_size_gb))
message(sprintf("Total size reduction of %.1f%%.", (1 - saved_size_gb / original_size_gb) * 100))

# _______________
# Memory Cleanup
# _______________

message("Cleaning up memory...")
rm(list = ls())
gc()

message("Script finished.")