# ---------------------------------------------------
# MERFISH Mouse Colon IBD - Create Azimuth Reference
# ---------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(Azimuth)
library(Seurat)
library(Matrix)

# ----------------------------------------------

# _______________________________________
# Input Parameters for Azimuth Reference
# _______________________________________

# Read in reference data (Healthy MERFISH)
# Note: Ensure this path matches where you saved the processed healthy data
input_sce_name <- "healthy_data"  
healthy_data_path <- "data/merfish/healthy_data.rds"

if(file.exists(healthy_data_path)) {
    healthy_data <- readRDS(healthy_data_path)
} else {
    stop("Healthy reference data not found at: ", healthy_data_path)
}

# Output directory
output_directory <- "data/merfish/Azimuth/custom_azimuth_reference"

# Cell type column to use (Specified by user)
celltype_column <- "tier2"  

# Setting the seed
set.seed(0)

# ______________________________________
# Computational Optimization Parameters 
# ______________________________________

# SCTransform: Use standard parameters for reproducibility
n_variable_features <- 3000  # Keep default for full gene coverage
sct_clip_range <- c(-10, 10)  # Standard clipping for numerical stability

# PCA: Keep standard dimensionality for downstream analyses
n_pca_dims <- 50  # Standard for single-cell analysis

# UMAP: Minimal modifications for computational efficiency with large datasets
umap_dims <- 1:30  # Use top 30 PCs (standard practice)
umap_n_neighbors <- 30  # Keep default for proper local structure
umap_min_dist <- 0.3  # Keep default for proper embedding
umap_n_epochs <- 300  # Reduced from 500 for efficiency
umap_metric <- "cosine"  # Cosine distance appropriate for high-dimensional data
umap_learning_rate <- 1.0  # Keep default

# AzimuthReference: Conservative modifications for large dataset compatibility
azimuth_dims <- 1:50  # Use default
azimuth_k_param <- 100  # Higher than default (31) to ensure compatibility with RunAzimuth

cat("=== Publication-Ready Azimuth Reference Creation (MERFISH) ===\n")

# _______________________
# Load and Prepare Data
# _______________________

cat("=== Loading and Preparing Data ===\n")
input_sce <- healthy_data
cat("Processing dataset:", ncol(input_sce), "cells x", nrow(input_sce), "genes\n")

# Create output directory
if(!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)

# ____________________
# Process SCE Object
# ____________________

cat("=== Processing SCE Object ===\n")

# Check cell_type column
if(!celltype_column %in% colnames(colData(input_sce))) {
    stop("Cell type column '", celltype_column, "' not found in colData.")
}

# Convert logcounts to counts if needed
if(!"counts" %in% assayNames(input_sce)) {
    cat("Converting logcounts to counts (required for Azimuth)...\n")
    logcounts_data <- logcounts(input_sce)
    
    # Simple reversal of log-normalization (assuming log2(x+1))
    # Note: If size factors were used, this is an approximation suitable for Azimuth's SCT step
    max_val <- max(logcounts_data)
    if(max_val > 20) {
        approx_counts <- exp(logcounts_data) - 1
    } else {
        approx_counts <- 2^(logcounts_data) - 1
    }
    
    approx_counts[approx_counts < 0] <- 0
    approx_counts <- round(approx_counts)
    assay(input_sce, "counts") <- approx_counts
}

# _____________________
# Create Seurat Object
# _____________________

cat("=== Creating Seurat Object ===\n")

counts_data <- assay(input_sce, "counts")
logcounts_data <- logcounts(input_sce)
metadata <- as.data.frame(colData(input_sce))
metadata$cell_type <- factor(metadata[[celltype_column]])

cat("Cell types identified:", length(unique(metadata$cell_type)), "types\n")

normal_seurat <- CreateSeuratObject(
    counts = counts_data,
    meta.data = metadata,
    project = "healthy_merfish_reference"
)

normal_seurat <- SetAssayData(normal_seurat, 
                             layer = "data", 
                             new.data = logcounts_data)

# ____________________________
# Standard Processing Pipeline
# ____________________________

cat("=== Processing Pipeline ===\n")

# SCTransform
cat("Running SCTransform...\n")
start_time <- Sys.time()
normal_seurat <- SCTransform(
    normal_seurat, 
    variable.features.n = n_variable_features, 
    clip.range = sct_clip_range,
    verbose = FALSE
)
cat("SCTransform completed.\n")

# PCA
cat("Running PCA...\n")
normal_seurat <- RunPCA(
    normal_seurat, 
    npcs = n_pca_dims, 
    verbose = FALSE
)

# UMAP
cat("Running UMAP...\n")
normal_seurat <- RunUMAP(
    normal_seurat, 
    dims = umap_dims, 
    n.neighbors = umap_n_neighbors, 
    min.dist = umap_min_dist, 
    n.epochs = umap_n_epochs, 
    metric = umap_metric, 
    learning.rate = umap_learning_rate,
    return.model = TRUE, 
    verbose = FALSE
)

# ____________________________
# Azimuth Reference Creation
# ____________________________

cat("=== Creating Azimuth Reference ===\n")

# Create plot metadata for reference
plot_metadata_df <- data.frame(
    cell_type = normal_seurat$cell_type,
    row.names = colnames(normal_seurat)
)

reference <- AzimuthReference(
    object = normal_seurat,
    refUMAP = "umap",
    refDR = "pca", 
    refAssay = "SCT",
    dims = azimuth_dims,
    k.param = azimuth_k_param, 
    plotref = "umap",
    plot.metadata = plot_metadata_df,
    ori.index = NULL,
    colormap = NULL,
    assays = NULL,
    metadata = c("cell_type"),
    reference.version = paste0("merfish_reference_", ncol(normal_seurat), "cells_v1.0"),
    verbose = FALSE
)

# _______________
# Save Reference
# _______________

cat("=== Saving Reference ===\n")

SaveAzimuthReference(
    object = reference,
    folder = paste0(output_directory, "/")
)

# ____________________
# Results and Cleanup
# ____________________

reference_info <- list(
    path = output_directory,
    n_cells = ncol(normal_seurat),
    cell_types = levels(normal_seurat$cell_type),
    creation_date = Sys.time(),
    version = paste0("merfish_reference_", ncol(normal_seurat), "cells_v1.0")
)

# Save reference info
saveRDS(reference_info, file.path(output_directory, "reference_info.rds"))

# Cleanup
rm(healthy_data, normal_seurat, reference, counts_data, logcounts_data, metadata, input_sce)
gc()

cat("\n✓ Publication-ready MERFISH Azimuth reference created successfully!\n")