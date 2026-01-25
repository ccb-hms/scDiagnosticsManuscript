# ----------------------------------------------------------------
# MERFISH Mouse Colon IBD - Add R (Azimuth & SingleR) Annotations
# ----------------------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(BiocParallel)
library(parallel)

# Source files (Ensure these paths are correct relative to your project root)
source("R/auxiliary/performSingleRWithSubsampling.R")
source("R/auxiliary/performAzimuthAnnotation.R")
source("R/auxiliary/addMergedCellTypes.R")

# -----------------------------------------------------

# Read processed datasets
# Using MERFISH paths
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

# Create output directories for MERFISH
if(!dir.exists("data/merfish/Azimuth")) dir.create("data/merfish/Azimuth", recursive = TRUE)
if(!dir.exists("data/merfish/SingleR")) dir.create("data/merfish/SingleR", recursive = TRUE)

# ____________________
# Azimuth Annotations
# ____________________

# Set path to your custom reference
custom_reference_path <- "data/merfish/Azimuth/custom_azimuth_reference"

# Verify the reference exists
if(!dir.exists(custom_reference_path)) {
    stop("Custom reference not found at: ", custom_reference_path, 
         "\nPlease run the Azimuth reference creation script first!")
}

cat("Using custom reference from:", custom_reference_path, "\n")

# Apply Azimuth annotations
cat("Running Azimuth annotation with custom reference...\n")
azimuth_dss9_output <- performAzimuthAnnotation(
    sce_obj = dss9_data,
    dataset_name = "dss9", 
    reference_path = custom_reference_path
)

# Extract annotated object
dss9_data <- azimuth_dss9_output$sce_annotated

cat("✓ Azimuth annotation complete!\n")

# _____________________
# SingleR Annotations
# _____________________

# SingleR: Healthy -> DSS9
# Note: Using 'tier2' as the annotation column from healthy_data
singler_output <- performSingleRWithSubsampling(
    ref_sce = healthy_data,
    query_sce = dss9_data,
    ref_name = "healthy_merfish",
    annotation_col = "tier2",
    max_cells_ref = NULL,
    bpparam = bpparam
)

dss9_data$singler_annotations <- singler_output$annotations
dss9_data$singler_scores <- singler_output$scores

# Set annotations as character data
dss9_data$singler_annotations <- as.character(dss9_data$singler_annotations)

# __________________________________
# Create Merged Annotation Columns
# __________________________________

message("--- Creating Merged Versions of Annotations ---")

# Use the flexible function to add merged columns directly to the SCE object
dss9_data <- addMergedCellTypes(sce_object = dss9_data, input_col_name = "singler_annotations", 
                                dataset = "MERFISH")
dss9_data <- addMergedCellTypes(sce_object = dss9_data, input_col_name = "azimuth_celltype_l1", 
                                dataset = "MERFISH")
dss9_data <- addMergedCellTypes(sce_object = dss9_data, input_col_name = "azimuth_celltype_l2", 
                                dataset = "MERFISH")

message("Merged annotation columns created successfully.")
message("Example of merged SingleR annotations:")
print(head(table(dss9_data$singler_annotations_merged)))

# __________________
# Save All Results
# __________________

# Save annotated SCE objects
# Saving back to MERFISH data directory
saveRDS(dss9_data, "data/merfish/dss9_data.rds") 

# Save SingleR outputs
saveRDS(singler_output, "data/merfish/SingleR/singler_dss9_output.rds")

# Save Azimuth outputs
saveRDS(azimuth_dss9_output$azimuth_results, "data/merfish/Azimuth/azimuth_dss9_results.rds")

# Save Seurat objects with Azimuth annotations
saveRDS(azimuth_dss9_output$seurat_object, "data/merfish/Azimuth/azimuth_dss9_seurat.rds")

# Save custom reference info
reference_info <- list(
    reference_path = custom_reference_path,
    reference_type = "custom_healthy_merfish",
    creation_date = Sys.Date()
)
saveRDS(reference_info, "data/merfish/Azimuth/reference_info.rds")

# Clean up
rm(healthy_data, dss9_data, 
   singler_output,
   azimuth_dss9_output,
   custom_reference_path, reference_info)

gc()
message("Pipeline finished. Annotated MERFISH objects saved.")