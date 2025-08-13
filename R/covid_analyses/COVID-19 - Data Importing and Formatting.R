# -----------------------------------------------
# COVID-19 - Data Importing and Formatting
# -----------------------------------------------

# Load libraries
library(zellkonverter)
library(SingleCellExperiment)
library(HDF5Array)
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)

# Source files
source("R/auxiliary/convertEnsemblToSymbols.R")

# -----------------------------------------------

# _____________________________
# Data Loading and Initial QC
# _____________________________

# Load with HDF5 backend (data stays on disk)
sce_full <- readH5AD("data/covid_data.h5ad", use_hdf5 = TRUE, reader = "R")
message("Initial data loaded with ", ncol(sce_full), " cells and ", nrow(sce_full), " features.")

# Apply basic QC filters
qc_pass <- sce_full$n_genes >= 500 & 
    sce_full$n_genes <= 8000 & 
    sce_full$pct_counts_mt <= 20
sce_full <- sce_full[, qc_pass]
message("After standard QC: ", ncol(sce_full), " cells remain.")


# _____________________
# Biological Filtering
# _____________________

# Remove specific experimental conditions and diseases
sce_full <- sce_full[, !(sce_full$Status_on_day_collection_summary %in% 
    c("LPS_10hours", "LPS_90mins"))]
sce_full <- sce_full[, sce_full$disease != "respiratory system disorder"]
message("After biological filtering: ", ncol(sce_full), " cells remain.")


# _____________________________________________
# Cell Type Eligibility and Balanced Sampling
# _____________________________________________

set.seed(0)

# Get cell metadata for easier manipulation
cell_data <- colData(sce_full) |>
    as.data.frame() |>
    tibble::rownames_to_column("cell_id")

# A cell type is eligible if it has >= 2500 cells in BOTH normal and COVID-19 groups.
MIN_CELLS_PER_GROUP <- 2500

cell_counts_by_disease <- cell_data |>
    group_by(disease, cell_type) |>
    summarise(n = n(), .groups = "drop")

eligible_types_table <- cell_counts_by_disease |>
    pivot_wider(names_from = disease, values_from = n, values_fill = 0) |>
    filter(normal >= MIN_CELLS_PER_GROUP, `COVID-19` >= MIN_CELLS_PER_GROUP)

eligible_types <- eligible_types_table$cell_type

message("Found ", length(eligible_types), " cell types with at least ", 
        MIN_CELLS_PER_GROUP, " cells in both normal and COVID-19 groups.")
message("Eligible types: ", paste(eligible_types, collapse = ", "))

SAMPLES_PER_TYPE <- 2500

# Sample for the 'normal' group
sampled_normal_ids <- lapply(eligible_types, function(ct) {
    ids <- cell_data$cell_id[cell_data$disease == "normal" & cell_data$cell_type == ct]
    sample(ids, SAMPLES_PER_TYPE)
}) |> unlist()

# Sample for the 'COVID-19' group
sampled_covid_ids <- lapply(eligible_types, function(ct) {
    ids <- cell_data$cell_id[cell_data$disease == "COVID-19" & cell_data$cell_type == ct]
    sample(ids, SAMPLES_PER_TYPE)
}) |> unlist()

# Create the final subset objects
sce_normal <- sce_full[, sampled_normal_ids]
sce_covid <- sce_full[, sampled_covid_ids]

message("\nCreated 'normal' dataset with ", ncol(sce_normal), " cells.")
message("Created 'COVID-19' dataset with ", ncol(sce_covid), " cells.")


# ____________________________
# Final Formatting and Saving
# ____________________________

# Define a function to apply final formatting to avoid repeating code
formatSCE <- function(sce_object) {
    # Keep only selected metadata columns
    selected_cols <- c("sample_id", "disease", "cell_type", "Site", "Days_from_onset",
                       "sex", "Smoker", "donor_id", "development_stage", "Status_on_day_collection_summary")
    colData(sce_object) <- colData(sce_object)[, selected_cols]
    
    # Remove existing embeddings
    reducedDims(sce_object) <- list()

    # Convert assay to a memory-efficient sparse matrix
    assay(sce_object) <- as(realize(assay(sce_object)), "CsparseMatrix")
    names(assays(sce_object)) <- "logcounts"

    # Clean up rowData and convert gene IDs to symbols
    rowData(sce_object) <- rowData(sce_object)[, "feature_name", drop = FALSE]
    sce_object <- convertEnsemblToSymbols(sce_object)
    
    return(sce_object)
}

# Apply formatting to both objects
sce_normal_final <- formatSCE(sce_normal)
sce_covid_final <- formatSCE(sce_covid)

# Save the final objects
saveRDS(sce_normal_final, "data/normal_data.rds")
saveRDS(sce_covid_final, "data/covid_data.rds")
message("Saved final objects to 'data/normal_data.rds' and 'data/covid_data.rds'")


# _______________
# Final Summary
# _______________

message("\n--- Pipeline Complete! ---")
message("Final 'normal' dataset cell type distribution (perfectly balanced):")
print(table(sce_normal_final$cell_type))

message("\nFinal 'COVID-19' dataset cell type distribution (perfectly balanced):")
print(table(sce_covid_final$cell_type))

# Clean up memory
rm(sce_full, sce_normal, sce_covid, sce_normal_final, sce_covid_final, cell_data)
gc()