# -----------------------------------------------
# COVID-19 - Data Cleaning and Metadata Removal
# -----------------------------------------------

# Load libraries
library(zellkonverter)
library(SingleCellExperiment)
library(HDF5Array)
library(Matrix)
library(dplyr)

# -----------------------------------------------

cat(paste(rep("=", 45), collapse = ""), "\n")
cat("COVID-19 DATA CLEANING - REMOVING METADATA\n")
cat(paste(rep("=", 45), collapse = ""), "\n")

# _____________________________
# File Information and Loading
# _____________________________

input_path <- "data/covid/covid_data.h5ad"
output_path <- "data/covid/covid_data_clean.rds"

cat("\nLoading original file:", input_path, "\n")
original_size <- file.info(input_path)$size / (1024^3)
cat(sprintf("Original file size: %.2f GB\n", original_size))

# Load the original file
sce_full <- readH5AD(input_path, use_hdf5 = TRUE, reader = "R")
cat(sprintf("Data dimensions: %d cells × %d genes\n", ncol(sce_full), nrow(sce_full)))

cat(sprintf("Memory usage before cleaning: %.2f GB\n", 
    as.numeric(object.size(sce_full)) / (1024^3)))

# ______________
# Data Cleaning
# ______________

cat("\nCleaning metadata but keeping DelayedMatrix format...\n")

# Create SCE object without heavy metadata but keep DelayedMatrix
sce_cleaned <- SingleCellExperiment(
    assays = list(X = assay(sce_full, "X")),  # Keep as DelayedMatrix for now
    colData = colData(sce_full),  # Keep all colData
    rowData = rowData(sce_full)   # Keep all rowData
)

# Clear the heavy metadata components
metadata(sce_cleaned) <- list()
reducedDims(sce_cleaned) <- list()
altExps(sce_cleaned) <- list()

cat("✓ Removed all metadata (neighbors, antibody data, PCA, etc.)\n")
cat("✓ Kept DelayedMatrix format for efficient sampling\n")

# ________________________
# Cell Sampling by Sample
# ________________________

cat("\nFiltering and sampling cells by condition...\n")

# First, filter out respiratory system disorder cases
cat("Filtering out respiratory system disorder cases...\n")
keep_indices <- sce_cleaned$disease %in% c("normal", "COVID-19")
sce_filtered <- sce_cleaned[, keep_indices]

# Second, filter out LPS-treated samples from healthy controls
cat("Filtering out LPS-treated samples from healthy controls...\n")
healthy_samples <- (sce_filtered$disease == "normal" & 
                   sce_filtered$Status_on_day_collection_summary == "Healthy")
covid_samples <- sce_filtered$disease == "COVID-19"
keep_clean <- healthy_samples | covid_samples

sce_filtered <- sce_filtered[, keep_clean]

cat(sprintf("Filtered data: %d cells retained\n", ncol(sce_filtered)))
cat("Disease distribution after filtering:\n")
print(table(sce_filtered$disease))
cat("Status distribution for normal samples:\n")
normal_status <- table(sce_filtered$Status_on_day_collection_summary[sce_filtered$disease == "normal"])
print(normal_status)

# Check sample distribution after filtering
sample_counts <- table(sce_filtered$sample_id)
cat(sprintf("Total samples after filtering: %d\n", length(sample_counts)))
cat(sprintf("Cells per sample (range): %d - %d\n", min(sample_counts), max(sample_counts)))

# Filter samples with at least 500 cells
samples_with_500plus <- names(sample_counts)[sample_counts >= 500]
cat(sprintf("Samples with ≥500 cells: %d\n", length(samples_with_500plus)))

# Check disease status for each qualifying sample
disease_status <- sce_filtered$disease[match(samples_with_500plus, sce_filtered$sample_id)]
normal_samples <- samples_with_500plus[disease_status == "normal"]
covid_samples <- samples_with_500plus[disease_status == "COVID-19"]

cat(sprintf("  - Normal samples: %d\n", length(normal_samples)))
cat(sprintf("  - COVID samples: %d\n", length(covid_samples)))

# Sample cells with different strategies by condition
set.seed(0)  # For reproducibility
sampled_indices <- c()

# Sample from normal samples (1000 cells or all if <1000)
for(sample in normal_samples) {
    sample_indices <- which(sce_filtered$sample_id == sample)
    n_cells <- length(sample_indices)
    n_sample <- ifelse(n_cells >= 1000, 1000, n_cells)  # Take all if <1000, otherwise 1000
    sampled_indices <- c(sampled_indices, sample(sample_indices, n_sample))
}

# Sample from COVID samples (500 cells each)
for(sample in covid_samples) {
    sample_indices <- which(sce_filtered$sample_id == sample)
    sampled_indices <- c(sampled_indices, sample(sample_indices, 500))
}

# Subset the SCE object
sce_sampled <- sce_filtered[, sampled_indices]

# Report sampling results
normal_cells <- sum(sce_sampled$disease == "normal")
covid_cells <- sum(sce_sampled$disease == "COVID-19")

cat(sprintf("✓ Sampled from %d normal samples (up to 1000 cells each)\n", length(normal_samples)))
cat(sprintf("✓ Sampled from %d COVID samples (500 cells each)\n", length(covid_samples)))
cat(sprintf("✓ Total normal cells: %d\n", normal_cells))
cat(sprintf("✓ Total COVID cells: %d\n", covid_cells))
cat(sprintf("✓ Total sampled cells: %d\n", ncol(sce_sampled)))
cat(sprintf("New dimensions: %d cells × %d genes\n", ncol(sce_sampled), nrow(sce_sampled)))

# _________________________
# Sparse Matrix Conversion
# _________________________

cat("\nConverting to sparse matrix...\n")
cat("Step 1: Realizing matrix (this may take a moment)...\n")
realized_matrix <- realize(assay(sce_sampled, "X"))
cat(sprintf("Realized matrix class: %s\n", class(realized_matrix)))

cat("Step 2: Converting to CsparseMatrix...\n")
expr_matrix <- as(realized_matrix, "CsparseMatrix")
cat(sprintf("Final matrix class: %s\n", class(expr_matrix)))
cat(sprintf("Sparse matrix size: %.2f GB\n", as.numeric(object.size(expr_matrix)) / (1024^3)))

# Create final SCE object with sparse matrix
sce_clean <- SingleCellExperiment(
    assays = list(X = expr_matrix),
    colData = colData(sce_sampled),
    rowData = rowData(sce_sampled)
)

cat("✓ Converted to memory-efficient sparse matrix format\n")
cat(sprintf("Final memory usage: %.2f GB\n", as.numeric(object.size(sce_clean)) / (1024^3)))

# ____________
# File Saving
# ____________

cat("\nSaving cleaned SCE object as RDS...\n")

# Save as RDS file
saveRDS(sce_clean, output_path, compress = TRUE)

# Check the saved file size
saved_size <- file.info(output_path)$size / (1024^3)

cat("✓ File saved successfully\n")
cat(sprintf("✓ Saved file size: %.2f GB\n", saved_size))
cat(sprintf("✓ Space reduction: %.2f GB (%.1f%% of original)\n", 
    original_size - saved_size, 
    (saved_size/original_size)*100))

cat(sprintf("\nCleaned data saved to: %s\n", output_path))
cat("Original h5ad file preserved unchanged.\n")

cat("\n", paste(rep("=", 45), collapse = ""), "\n")
cat("DATA CLEANING COMPLETED SUCCESSFULLY\n")
cat(paste(rep("=", 45), collapse = ""), "\n")

# _______________
# Memory Cleanup
# _______________

cat("\nCleaning up memory...\n")

# Remove all objects from environment
rm(list = ls())

# Force garbage collection
gc()

cat("✓ Memory cleared\n")

# -----------------------------------------------