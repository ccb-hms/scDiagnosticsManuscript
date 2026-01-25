# --------------------------------------------------------------
# MERFISH Mouse Colon IBD - Python Annotation Data Preparation
# --------------------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)

# Source files 
source("R/auxiliary/environmentSetupCellTypist.R")

# --------------------------------------------------------------

# Setting the seed
set.seed(0)

# Setup Python environment
cat("Setting up Python environment...\n")
environmentSetupCellTypist()

# Configuration
REFERENCE_FILE <- "data/merfish/healthy_data.rds"
QUERY_FILE <- "data/merfish/dss9_data.rds"
CELL_TYPE_COLUMN <- "tier2"  

# Output files
ANNDATA_REF_FILE <- "data/merfish/reference_data.h5ad"
ANNDATA_QUERY_FILE <- "data/merfish/query_data.h5ad"

cat("Loading SCE objects...\n")
normal_data <- readRDS(REFERENCE_FILE)
dss9_data <- readRDS(QUERY_FILE)

cat("Reference (Healthy) dimensions:", dim(normal_data), "\n")
cat("Query (DSS9) dimensions:", dim(dss9_data), "\n")

# Convert logcounts to approximate counts and remove logcounts
cleanSCEObject <- function(sce_obj, obj_name) {
  if(!"counts" %in% names(assays(sce_obj))) {
    if("logcounts" %in% names(assays(sce_obj))) {
      cat("Converting logcounts to approximate counts for", obj_name, "...\n")
      
      logcounts_data <- logcounts(sce_obj)
      
      # Detect log base and convert back to counts
      max_val <- max(logcounts_data)
      if(max_val > 20) {
        cat("Detected natural log normalization\n")
        approx_counts <- exp(logcounts_data) - 1
      } else {
        cat("Detected log2 normalization\n")
        approx_counts <- 2^(logcounts_data) - 1
      }
      
      # Clean up negative values and round
      approx_counts[approx_counts < 0] <- 0
      approx_counts <- round(approx_counts)
      
      # Add as counts assay
      assay(sce_obj, "counts") <- approx_counts
      cat("✓ Added 'counts' assay to", obj_name, "\n")
    } else {
      stop("Neither 'counts' nor 'logcounts' assay found in ", obj_name, "!")
    }
  }
  
  # Remove ALL other assays (keep only counts)
  assay_names <- names(assays(sce_obj))
  for(assay_name in assay_names) {
    if(assay_name != "counts") {
      assays(sce_obj)[[assay_name]] <- NULL
      cat("✓ Removed '", assay_name, "' assay from ", obj_name, "\n", sep = "")
    }
  }
  
  return(sce_obj)
}

# Convert both objects and clean up assays
normal_data <- cleanSCEObject(normal_data, "reference")
dss9_data <- cleanSCEObject(dss9_data, "query")

# Check if cell type column exists in reference
if (!CELL_TYPE_COLUMN %in% colnames(colData(normal_data))) {
  available_cols <- colnames(colData(normal_data))
  cat("Available colData columns:", paste(available_cols, collapse = ", "), "\n")
  stop("Cell type column '", CELL_TYPE_COLUMN, "' not found in reference data!")
}

cat("Cell types in reference:\n")
print(table(colData(normal_data)[[CELL_TYPE_COLUMN]]))

# AGGRESSIVE colData cleaning
cat("Cleaning colData for maximum efficiency...\n")

# Reference: Keep ONLY cell_type column
cat("Reference before cleaning:", ncol(colData(normal_data)), "colData columns\n")
ref_coldata_minimal <- colData(normal_data)[, CELL_TYPE_COLUMN, drop = FALSE]
colData(normal_data) <- ref_coldata_minimal
cat("Reference after cleaning:", ncol(colData(normal_data)), "colData columns\n")

# Query: Remove ALL colData (CellTypist doesn't need any metadata for prediction)
cat("Query before cleaning:", ncol(colData(dss9_data)), "colData columns\n")
colData(dss9_data) <- DataFrame(row.names = colnames(dss9_data))
cat("Query after cleaning:", ncol(colData(dss9_data)), "colData columns\n")

# Clean rowData
rowData(normal_data) <- DataFrame(row.names = rownames(normal_data))
rowData(dss9_data) <- DataFrame(row.names = rownames(dss9_data))

# Remove any reducedDims
reducedDims(normal_data) <- list()
reducedDims(dss9_data) <- list()

cat("Converting ultra-lean SCE objects to AnnData...\n")

# Convert to AnnData
adata_ref <- SCE2AnnData(normal_data, X_name = "counts")
adata_query <- SCE2AnnData(dss9_data, X_name = "counts")

# Save AnnData objects
cat("Saving optimized AnnData files...\n")
adata_ref$write_h5ad(ANNDATA_REF_FILE)
adata_query$write_h5ad(ANNDATA_QUERY_FILE)

# Show file sizes
ref_size_mb <- round(file.size(ANNDATA_REF_FILE) / 1024^2, 1)
query_size_mb <- round(file.size(ANNDATA_QUERY_FILE) / 1024^2, 1)

cat("Data preparation complete!\n")
cat("Files created:\n")
cat("- ", ANNDATA_REF_FILE, " (", ref_size_mb, " MB)\n", sep = "")
cat("- ", ANNDATA_QUERY_FILE, " (", query_size_mb, " MB)\n", sep = "")

# Clean up memory
rm(normal_data, dss9_data, adata_ref, adata_query)
gc()

cat("✓ Ultra-lean data preparation complete!\n")