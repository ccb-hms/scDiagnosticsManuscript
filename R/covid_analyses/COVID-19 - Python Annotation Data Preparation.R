# ----------------------------------------------
# COVID-19 - Python Annotation Data Preparation
# ----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)

# ----------------------------------------------

# Setup Python environment
cat("Setting up Python environment...\n")
source("R/auxiliary/EnvironmentSetupCellTypist.R")

# Configuration
REFERENCE_FILE <- "data/normal_data.rds"
QUERY_FILE <- "data/covid_data.rds"
CELL_TYPE_COLUMN <- "cell_type"  

# Output files
ANNDATA_REF_FILE <- "data/reference_data.h5ad"
ANNDATA_QUERY_FILE <- "data/query_data.h5ad"

cat("Loading SCE objects...\n")
normal_data <- readRDS(REFERENCE_FILE)
covid_data <- readRDS(QUERY_FILE)

cat("Reference data dimensions:", dim(normal_data), "\n")
cat("Query data dimensions:", dim(covid_data), "\n")

# Check for counts assay
if(!"counts" %in% names(assays(normal_data))) {
  cat("Available assays in reference:", names(assays(normal_data)), "\n")
  stop("'counts' assay not found in reference data!")
}

if(!"counts" %in% names(assays(covid_data))) {
  cat("Available assays in query:", names(assays(covid_data)), "\n")
  stop("'counts' assay not found in query data!")
}

cat("✓ Both datasets have 'counts' assay - perfect for CellTypist!\n")

# Check if cell type column exists in reference
if (!CELL_TYPE_COLUMN %in% colnames(colData(normal_data))) {
  stop("Cell type column '", CELL_TYPE_COLUMN, "' not found in reference data!")
}

cat("Cell types in reference:\n")
print(table(colData(normal_data)[[CELL_TYPE_COLUMN]]))

# Prepare reference data - keep only cell type annotation
cat("Preparing reference data - keeping only cell type annotation...\n")
ref_coldata_minimal <- colData(normal_data)[, CELL_TYPE_COLUMN, drop = FALSE]
colData(normal_data) <- ref_coldata_minimal

# Prepare query data - completely clear colData (we'll reload original later)
cat("Preparing query data - clearing colData for efficiency...\n")
colData(covid_data) <- DataFrame(row.names = colnames(covid_data))

cat("Converting SCE objects to AnnData using raw counts...\n")

# Convert to AnnData using counts
adata_ref <- SCE2AnnData(normal_data, X_name = "counts")
adata_query <- SCE2AnnData(covid_data, X_name = "counts")

# Save AnnData objects
cat("Saving AnnData files...\n")
adata_ref$write_h5ad(ANNDATA_REF_FILE)
adata_query$write_h5ad(ANNDATA_QUERY_FILE)

cat("Data preparation complete!\n")
cat("Files created:\n")
cat("- ", ANNDATA_REF_FILE, "\n")
cat("- ", ANNDATA_QUERY_FILE, "\n")
cat("- Used assay: counts (optimal for CellTypist)\n")
cat("✓ No need to modify Python script - it will work perfectly with raw counts!\n")

cat("\nNote: Original covid_data.rds will be reloaded during result merging.\n")