# -----------------------------------------------
# COVID-19 PBMC - Cell Type Merging
# -----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(dplyr)

# -----------------------------------------------

# Create a mapping function to merge cell types
addMergedCellTypes <- function(sce_object, input_col_name) {
    
    # Input Validation
    if (!input_col_name %in% colnames(colData(sce_object))) {
        stop(sprintf("Error: Column '%s' not found in the colData of the SCE object.", input_col_name))
    }
    
    # Get Labels and Apply Merging Logic
    labels_vector <- sce_object[[input_col_name]]
    merged_types <- as.character(labels_vector) # Start with a copy
    
    # Define mappings
    cd14_mono_subtypes <- c("CD14_mono", "CD83_CD14_mono")
    cd16_mono_subtypes <- c("CD16_mono", "C1_CD16_mono")
    cd4_subtypes <- c("CD4.CM", "CD4.EM", "CD4.IL22", "CD4.Naive", "CD4.Prolif", "CD4.Tfh", "CD4.Th1", "CD4.Th17", "CD4.Th2")
    cd8_subtypes <- c("CD8.EM", "CD8.Naive", "CD8.Prolif", "CD8.TE")
    b_subtypes <- c("B_exhausted", "B_immature", "B_naive", "B_non-switched_memory", "B_switched_memory", "B_malignant")
    plasma_subtypes <- c("Plasma_cell_IgA", "Plasma_cell_IgG", "Plasma_cell_IgM")
    hspc_subtypes <- c("HSC_CD38neg", "HSC_CD38pos", "HSC_erythroid", "HSC_MK", "HSC_myeloid", "HSC_prolif")
    ilc_subtypes <- c("ILC1_3", "ILC2")

    # Apply mappings
    merged_types[merged_types %in% cd14_mono_subtypes] <- "CD14 mono"
    merged_types[merged_types %in% cd16_mono_subtypes] <- "CD16 mono"
    merged_types[merged_types %in% cd4_subtypes] <- "CD4 T"
    merged_types[merged_types %in% cd8_subtypes] <- "CD8 T"
    merged_types[merged_types %in% b_subtypes] <- "B cell"
    merged_types[merged_types %in% plasma_subtypes] <- "Plasma cell"
    merged_types[merged_types %in% hspc_subtypes] <- "HSPC"
    merged_types[merged_types %in% ilc_subtypes] <- "ILC"

    # Add New Column and Return
    output_col_name <- paste0(input_col_name, "_merged")
    sce_object[[output_col_name]] <- merged_types
    
    message(sprintf("✓ Added merged cell types to new column: '%s'", output_col_name))
    
    return(sce_object)
}
