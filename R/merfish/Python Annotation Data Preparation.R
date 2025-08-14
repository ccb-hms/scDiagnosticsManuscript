# -------------------------------------------------------------
# MERFISH Mouse Colon IBD - Python Annotation Data Preparation
# -------------------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)

# Source files 
source("R/auxiliary/environmentSetupCellTypist.R")

# -------------------------------------------------------------

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
healthy_data <- readRDS(REFERENCE_FILE)
dss9_data <- readRDS(QUERY_FILE)

cat("Reference data dimensions:", dim(healthy_data), "\n")
cat("Query data dimensions:", dim(dss9_data), "\n")

# Check available assays and colData
cat("Available assays in reference:", names(assays(healthy_data)), "\n")
cat("Available assays in query:", names(assays(dss9_data)), "\n")
cat("colData columns in reference:", ncol(colData(healthy_data)), "\n")
cat("colData columns in query:", ncol(colData(dss9_data)), "\n")

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
healthy_data <- cleanSCEObject(healthy_data, "reference")
dss9_data <- cleanSCEObject(dss9_data, "query")

# Check if cell type column exists in reference
if (!CELL_TYPE_COLUMN %in% colnames(colData(healthy_data))) {
  available_cols <- colnames(colData(healthy_data))
  cat("Available colData columns:", paste(available_cols, collapse = ", "), "\n")
  stop("Cell type column '", CELL_TYPE_COLUMN, "' not found in reference data!")
}

cat("Cell types in reference:\n")
print(table(colData(healthy_data)[[CELL_TYPE_COLUMN]]))

# AGGRESSIVE colData cleaning
cat("Cleaning colData for maximum efficiency...\n")

# Reference: Keep ONLY cell_type column
cat("Reference before cleaning:", ncol(colData(healthy_data)), "colData columns\n")
ref_coldata_minimal <- colData(healthy_data)[, CELL_TYPE_COLUMN, drop = FALSE]
colData(healthy_data) <- ref_coldata_minimal
cat("Reference after cleaning:", ncol(colData(healthy_data)), "colData columns\n")

# Query: Remove ALL colData (CellTypist doesn't need any metadata for prediction)
cat("Query before cleaning:", ncol(colData(dss9_data)), "colData columns\n")
colData(dss9_data) <- DataFrame(row.names = colnames(dss9_data))
cat("Query after cleaning:", ncol(colData(dss9_data)), "colData columns\n")

# Also clean up rowData (gene metadata) - CellTypist only needs gene names
cat("Cleaning rowData (gene metadata)...\n")
cat("Reference rowData before:", ncol(rowData(healthy_data)), "columns\n")
cat("Query rowData before:", ncol(rowData(dss9_data)), "columns\n")

# Keep minimal rowData (just the essentials)
rowData(healthy_data) <- DataFrame(row.names = rownames(healthy_data))
rowData(dss9_data) <- DataFrame(row.names = rownames(dss9_data))

cat("Reference rowData after:", ncol(rowData(healthy_data)), "columns\n")
cat("Query rowData after:", ncol(rowData(dss9_data)), "columns\n")

# Remove any reducedDims (PCA, UMAP, etc.) - not needed for CellTypist
reducedDims(healthy_data) <- list()
reducedDims(dss9_data) <- list()
cat("✓ Removed all reducedDims (not needed for CellTypist)\n")

cat("Converting ultra-lean SCE objects to AnnData...\n")

# Convert to AnnData (now super clean and minimal)
adata_ref <- SCE2AnnData(healthy_data, X_name = "counts")
adata_query <- SCE2AnnData(dss9_data, X_name = "counts")

# Save AnnData objects (now much smaller!)
cat("Saving optimized AnnData files...\n")
adata_ref$write_h5ad(ANNDATA_REF_FILE)
adata_query$write_h5ad(ANNDATA_QUERY_FILE)

# Show file sizes
ref_size_mb <- round(file.size(ANNDATA_REF_FILE) / 1024^2, 1)
query_size_mb <- round(file.size(ANNDATA_QUERY_FILE) / 1024^2, 1)

cat("Data preparation complete!\n")
cat("Files created (ultra-optimized):\n")
cat("- ", ANNDATA_REF_FILE, " (", ref_size_mb, " MB)\n", sep = "")
cat("- ", ANNDATA_QUERY_FILE, " (", query_size_mb, " MB)\n", sep = "")
cat("- Reference: counts + cell_type only\n")
cat("- Query: counts only (no metadata)\n")

# Clean up memory
cat("Cleaning up memory...\n")
rm(healthy_data, dss9_data, adata_ref, adata_query)
gc()

cat("✓ Ultra-lean data preparation complete!\n")