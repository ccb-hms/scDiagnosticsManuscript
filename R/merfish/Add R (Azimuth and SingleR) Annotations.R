# -----------------------------------------------------------------
# MERFISH Mouse Colon IBD - Add R (Azimuth & SingleR) Annotations
# -----------------------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(BiocParallel)
library(parallel)

# Source files
source("R/auxiliary/performSingleRWithSubsampling.R")
source("R/auxiliary/performAzimuthAnnotation.R")

# -----------------------------------------------------------------

# Read all processed datasets
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds") 

# Setting the seed
set.seed(0)

# More stable SnowParam setup for Windows
if(.Platform$OS.type == "windows") {
    bpparam <- SnowParam(
        workers = min(detectCores() - 4, 6)              
    )
} else {
    bpparam <- MulticoreParam(workers = min(detectCores() - 4, 6))
}

# Create output directories
if(!dir.exists("data/merfish/Azimuth")) dir.create("data/merfish/Azimuth", recursive = TRUE)
if(!dir.exists("data/merfish/SingleR")) dir.create("data/merfish/SingleR", recursive = TRUE)

# ____________________
# Azimuth Annotations
# ____________________

# Set path to your custom reference (created by the reference creation script)
custom_reference_path <- "data/merfish/Azimuth/custom_azimuth_reference"

# Verify the reference exists
if(!dir.exists(custom_reference_path)) {
    stop("Custom reference not found at: ", custom_reference_path, 
         "\nPlease run the Azimuth reference creation script first!")
}

# Check if reference files exist
ref_files <- c(file.path(custom_reference_path, "ref.Rds"),
               file.path(custom_reference_path, "idx.annoy"))

if(!all(file.exists(ref_files))) {
    missing_files <- ref_files[!file.exists(ref_files)]
    stop("Missing reference files: ", paste(missing_files, collapse = ", "),
         "\nPlease run the Azimuth reference creation script first!")
}

cat("Using custom reference from:", custom_reference_path, "\n")

# Load reference info if available
ref_info_path <- file.path(custom_reference_path, "reference_info.rds")
if(file.exists(ref_info_path)) {
    ref_info <- readRDS(ref_info_path)
    cat("Reference info:\n")
    cat("- Cells:", ref_info$n_cells, "\n")
    cat("- Cell types:", length(ref_info$cell_types), "\n")
    cat("- Created:", ref_info$creation_date, "\n")
}

# Choose reference to use
reference_to_use <- custom_reference_path

# Option 2: Use standard PBMC reference instead
# reference_to_use <- "pbmcref"

# Option 3: Use path to downloaded reference
# reference_to_use <- "/path/to/downloaded/reference"

# Apply Azimuth annotations
cat("Running Azimuth annotation with custom reference...\n")
azimuth_merfish_output <- performAzimuthAnnotation(
    sce_obj = dss9_data,
    dataset_name = "merfish", 
    reference_path = reference_to_use
)
dss9_data <- azimuth_merfish_output$sce_annotated

cat("✓ Azimuth annotation complete!\n")

# _____________________
# SingleR Annotations
# _____________________

# SingleR: Healthy -> DSS9
singler_healthy_output <- performSingleRWithSubsampling(
    ref_sce = healthy_data,
    query_sce = dss9_data,
    ref_name = "healthy",
    annotation_col = "cell_type",
    max_cells_ref = NULL,
    bpparam = bpparam
)
dss9_data$singler_annotations <- singler_healthy_output$annotations
dss9_data$singler_scores <- singler_healthy_output$scores

# Set annotations as character data
dss9_data$singler_annotations <- as.character(dss9_data$singler_annotations)

# __________________
# SAVE ALL RESULTS
# __________________

# Save annotated SCE objects
saveRDS(dss9_data, "data/merfish/dss9_data.rds") 

# Save SingleR outputs
saveRDS(singler_healthy_output, "data/merfish/SingleR/singler_healthy_output.rds")

# Save Azimuth outputs
saveRDS(azimuth_merfish_output$azimuth_results, "data/merfish/Azimuth/azimuth_merfish_results.rds")

# Save Seurat objects with Azimuth annotations
saveRDS(azimuth_merfish_output$seurat_object, "data/merfish/Azimuth/azimuth_merfish_seurat.rds")

# Save custom reference info
reference_info <- list(
    reference_path = reference_to_use,
    reference_type = if(reference_to_use == custom_reference_path) "custom_healthy" else "standard",
    creation_date = Sys.Date()
)
saveRDS(reference_info, "data/merfish/Azimuth/reference_info.rds")

# Clean up
rm(healthy_data, dss9_data, 
   singler_healthy_output,
   azimuth_merfish_output,
   custom_reference_path, reference_to_use, reference_info)
