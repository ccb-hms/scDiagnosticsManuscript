# ----------------------------------------------------
# COVID-19 - Data QC, Processing & Dataset Splitting
# ----------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scater)
library(scran)
library(biomaRt)
library(Matrix)

# Source files
source("R/auxiliary/convertEnsemblToSymbols.R")
source("R/auxiliary/addReferencePCA.R")
source("R/auxiliary/mergeCellTypes.R")

# ----------------------------------------------------

cat(paste(rep("=", 50), collapse = ""), "\n")
cat("COVID-19 DATA QC, PROCESSING & DATASET SPLITTING\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# _____________
# Data Loading
# _____________

cat("\nLoading cleaned COVID data...\n")
sce_clean <- readRDS("data/covid/covid_data_clean.rds")

cat(sprintf("Initial data dimensions: %d genes × %d cells\n", nrow(sce_clean), ncol(sce_clean)))
cat("Disease distribution:\n")
print(table(sce_clean$disease))

# __________________________
# Quality Control Filtering
# __________________________

cat("\nApplying quality control filters...\n")

# 1. Remove cells with <1000 mRNA molecules
cat("Filtering cells with <1000 total counts...\n")
keep_cells_counts <- sce_clean$total_counts >= 1000
cat(sprintf("Keeping %d/%d cells (%.1f%%)\n", 
    sum(keep_cells_counts), length(keep_cells_counts), 
    100*sum(keep_cells_counts)/length(keep_cells_counts)))

# 2. Remove genes expressed in <10 cells (optional but typical)
cat("Filtering lowly expressed genes...\n")
keep_genes <- rowSums(assay(sce_clean, "X") > 0) >= 10
cat(sprintf("Keeping %d/%d genes (%.1f%%)\n", 
    sum(keep_genes), length(keep_genes), 
    100*sum(keep_genes)/length(keep_genes)))

# Apply filters
sce_filtered <- sce_clean[keep_genes, keep_cells_counts]

cat(sprintf("After QC filtering: %d genes × %d cells\n", nrow(sce_filtered), ncol(sce_filtered)))

# ______________________
# Additional QC Metrics
# ______________________

cat("\nCalculating additional QC metrics...\n")

# Calculate mitochondrial gene percentage (if not already present)
if(!"pct_counts_mt" %in% colnames(colData(sce_filtered))) {
    mito_genes <- grep("^MT-|^mt-", rownames(sce_filtered), value = TRUE)
    if(length(mito_genes) > 0) {
        sce_filtered <- addPerCellQC(sce_filtered, subsets = list(Mito = mito_genes))
        cat(sprintf("Found %d mitochondrial genes\n", length(mito_genes)))
    } else {
        cat("No mitochondrial genes found with MT- prefix\n")
    }
} else {
    cat("Mitochondrial percentages already calculated\n")
}

# ________________________________
# Convert Ensembl to Gene Symbols
# ________________________________

cat("\nConverting Ensembl IDs to gene symbols...\n")
sce_symbols <- convertEnsemblToSymbols(sce_filtered)

# __________________
# Dataset Splitting
# __________________

cat("\nSplitting into normal and COVID datasets...\n")

# Split datasets
normal_indices <- sce_symbols$disease == "normal"
covid_indices <- sce_symbols$disease == "COVID-19"

normal_data_sce <- sce_symbols[, normal_indices]
covid_data_sce <- sce_symbols[, covid_indices]

cat(sprintf("Normal dataset: %d genes × %d cells\n", nrow(normal_data_sce), ncol(normal_data_sce)))
cat(sprintf("COVID dataset: %d genes × %d cells\n", nrow(covid_data_sce), ncol(covid_data_sce)))

# _____________________________
# Assay Naming and Verification
# _____________________________

cat("\nChecking and updating assay names...\n")

# Check data ranges to determine data type
normal_range <- range(assay(normal_data_sce, "X"))
covid_range <- range(assay(covid_data_sce, "X"))

cat(sprintf("Normal data range: %.3f to %.3f\n", normal_range[1], normal_range[2]))
cat(sprintf("COVID data range: %.3f to %.3f\n", covid_range[1], covid_range[2]))

# Based on ranges (0-8.7, 0-9.1), this appears to be log-normalized data
cat("Data appears to be log-normalized based on ranges\n")

# Rename assay to logcounts
assayNames(normal_data_sce) <- "logcounts"
assayNames(covid_data_sce) <- "logcounts"

cat("✓ Renamed 'X' assay to 'logcounts'\n")

# _____________________________
# HVG and PCA on Both Datasets
# _____________________________

cat("\nComputing HVGs and PCA using diagnostic-optimized approach...\n")

# Data is already log-normalized, so skip logNormCounts step
cat("Using existing log-normalized data for HVG detection\n")

# Use union HVG approach for diagnostic-optimized PCA
normal_data_sce <- addReferencePCA(normal_data_sce, covid_data_sce, "Normal Reference")

# Extract diagnostic HVGs for COVID data PCA
diagnostic_hvgs <- metadata(normal_data_sce)$diagnostic_hvgs
cat(sprintf("Using %d diagnostic HVGs for COVID data PCA\n", length(diagnostic_hvgs)))

# Run PCA on COVID data using same diagnostic HVGs
covid_data_sce <- runPCA(covid_data_sce, subset_row = diagnostic_hvgs, ncomponents = 50)

cat("✓ PCA computed on both normal and COVID data\n")

# Store diagnostic HVG information in COVID data as well
metadata(covid_data_sce)$diagnostic_hvgs <- diagnostic_hvgs
metadata(covid_data_sce)$ref_hvgs <- metadata(normal_data_sce)$ref_hvgs
metadata(covid_data_sce)$query_hvgs <- metadata(normal_data_sce)$query_hvgs

cat("✓ Diagnostic HVG information stored in both datasets\n")

# Summary of HVG selection
cat("\nHVG Selection Summary:\n")
cat(sprintf("  Reference HVGs: %d\n", length(metadata(normal_data_sce)$ref_hvgs)))
cat(sprintf("  Query HVGs: %d\n", length(metadata(normal_data_sce)$query_hvgs))) 
cat(sprintf("  Diagnostic HVGs (union): %d\n", length(diagnostic_hvgs)))
cat(sprintf("  PCA components: %d\n", ncol(reducedDim(normal_data_sce, "PCA"))))

# _____________________________
# Cell Type Merging
# _____________________________

cat("\nMerging cell types for both datasets...\n")

# Apply cell type merging to both datasets
normal_data_sce$author_cell_type_merged <- mergeCellTypes(normal_data_sce)
covid_data_sce$author_cell_type_merged <- mergeCellTypes(covid_data_sce)

# Show results for normal dataset
cat("\nOriginal vs merged cell types (Normal dataset):\n")
original_normal_types <- table(normal_data_sce$author_cell_type)
merged_normal_types <- table(normal_data_sce$author_cell_type_merged)

cat("Original cell type count:", length(unique(normal_data_sce$author_cell_type)), "\n")
cat("Merged cell type count:", length(unique(normal_data_sce$author_cell_type_merged)), "\n")

cat("\nMerged cell type distribution (Normal):\n")
print(sort(merged_normal_types, decreasing = TRUE))

# Show results for COVID dataset
cat("\nMerged cell type distribution (COVID):\n")
merged_covid_types <- table(covid_data_sce$author_cell_type_merged)
print(sort(merged_covid_types, decreasing = TRUE))

cat("\nCell type merging summary:\n")
cat(sprintf("  Normal dataset: %d -> %d cell types\n", 
    length(unique(normal_data_sce$author_cell_type)), 
    length(unique(normal_data_sce$author_cell_type_merged))))
cat(sprintf("  COVID dataset: %d -> %d cell types\n", 
    length(unique(covid_data_sce$author_cell_type)), 
    length(unique(covid_data_sce$author_cell_type_merged))))

cat("✓ Cell type merging completed for both datasets\n")

# _________________
# Final Processing
# _________________

# Add dataset information
normal_data_sce$dataset_type <- "reference"
covid_data_sce$dataset_type <- "query"

cat("✓ Datasets processed and ready\n")

# ________________________
# Save Processed Datasets
# ________________________

cat("\nSaving processed datasets...\n")

saveRDS(normal_data_sce, "data/covid/normal_data_sce.rds")
saveRDS(covid_data_sce, "data/covid/covid_data_sce.rds")

cat("✓ Normal data saved to: data/covid/normal_data_sce.rds\n")
cat("✓ COVID data saved to: data/covid/covid_data_sce.rds\n")

# ___________________
# Summary Statistics
# ___________________

cat(paste(rep("=", 50), collapse = ""), "\n")
cat("PROCESSING SUMMARY\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

cat(sprintf("Reference (Normal) Dataset: %d cells × %d genes\n", 
    ncol(normal_data_sce), nrow(normal_data_sce)))
cat(sprintf("Query (COVID) Dataset: %d cells × %d genes\n", 
    ncol(covid_data_sce), nrow(covid_data_sce)))
cat(sprintf("Highly Variable Genes: %d\n", length(diagnostic_hvgs)))
cat(sprintf("PCA Components: %d\n", ncol(reducedDim(normal_data_sce, "PCA"))))

cat("\n✓ Data processing completed successfully!\n")

# _______________
# Memory Cleanup
# _______________

cat("\nCleaning up memory...\n")
rm(list = ls())
gc()
cat("✓ Memory cleared\n")

# -----------------------------------------------
