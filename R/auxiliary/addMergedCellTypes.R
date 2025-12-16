# -----------------------------------------------
# Multi-Dataset Cell Type Merging Function
# -----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(dplyr)

# -----------------------------------------------

# Create a mapping function to merge cell types
# Arguments:
#   sce_object: The SingleCellExperiment object
#   input_col_name: The column containing original annotations
#   dataset: Either "COVID" (default) or "MERFISH"

addMergedCellTypes <- function(sce_object, input_col_name, dataset = "COVID") {
    
    # Input Validation
    if (!input_col_name %in% colnames(colData(sce_object))) {
        stop(sprintf("Error: Column '%s' not found in the colData of the SCE object.", input_col_name))
    }
    
    if (!dataset %in% c("COVID", "MERFISH")) {
        stop("Error: 'dataset' argument must be either 'COVID' or 'MERFISH'.")
    }
    
    # Get Labels
    labels_vector <- as.character(sce_object[[input_col_name]])
    merged_types <- labels_vector # Initialize
    
    if (dataset == "COVID") {
        
        # Define mappings (Exact Matching)
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
        
    } 

    else if (dataset == "MERFISH") {
        
        # Apply mappings (Pattern Matching/Regex via case_when)
        # Note: Order matters in case_when (first match wins)
        merged_types <- dplyr::case_when(
            # FIBROBLAST LINEAGE
            grepl("^IAF", labels_vector) ~ "Inflamed Fibroblast",
            grepl("^Fibro", labels_vector) ~ "Fibroblast",
            
            # STROMAL OTHERS
            grepl("FRC|Pericyte", labels_vector, ignore.case = TRUE) ~ "Pericyte/FRC",
            
            # SMOOTH MUSCLE
            grepl("^IASMC", labels_vector) ~ "Inflamed SMC",
            grepl("SMC", labels_vector) ~ "Smooth Muscle",
            
            # EPITHELIAL & STEM
            grepl("^IAE", labels_vector) ~ "Inflamed Epithelial",
            grepl("Stem|TA", labels_vector) ~ "Stem/TA", 
            grepl("Colonocytes|Goblet|Epithelial|M cells|Repair associated", labels_vector) ~ "Epithelial",
            
            # IMMUNE
            grepl("Neutrophil", labels_vector, ignore.case = TRUE) ~ "Neutrophil",
            grepl("B cell|T|Macrophage|Monocyte|DC|ILC2|Plasma|Mast", labels_vector) ~ "Other Immune",
            
            # OTHERS
            grepl("Endothelial|EC", labels_vector) ~ "Endothelial",
            grepl("Glia|Neuron", labels_vector) ~ "Enteric Nervous",
            grepl("ICC", labels_vector) ~ "ICC",
            grepl("Adipose", labels_vector) ~ "Adipose",
            
            TRUE ~ "Other"
        )
    }

    # Add New Column and Return
    output_col_name <- paste0(input_col_name, "_merged")
    sce_object[[output_col_name]] <- merged_types
    
    message(sprintf("✓ Added merged cell types (Dataset: %s) to new column: '%s'", dataset, output_col_name))
    
    return(sce_object)
}