# -----------------------------------------------
# COVID-19 - Cell Type Merging 
# -----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(dplyr)

# -----------------------------------------------

# Create a mapping function to merge cell types
mergeCellTypes <- function(sce_object) {
    
    current_types <- as.character(sce_object$author_cell_type)
    merged_types <- current_types
    
    # Merge CD14+ monocyte subtypes -> "CD14 mono"
    cd14_mono_subtypes <- c("CD14_mono", "CD83_CD14_mono")
    merged_types[current_types %in% cd14_mono_subtypes] <- "CD14 mono"
    
    # Merge CD16+ monocyte subtypes -> "CD16 mono" 
    cd16_mono_subtypes <- c("CD16_mono", "C1_CD16_mono")
    merged_types[current_types %in% cd16_mono_subtypes] <- "CD16 mono"
    
    # CD4 T cells - merge all CD4 subtypes
    cd4_subtypes <- c("CD4.CM", "CD4.EM", "CD4.IL22", "CD4.Naive", 
                      "CD4.Prolif", "CD4.Tfh", "CD4.Th1", "CD4.Th17", "CD4.Th2")
    merged_types[current_types %in% cd4_subtypes] <- "CD4 T"
    
    # CD8 T cells - merge all CD8 subtypes  
    cd8_subtypes <- c("CD8.EM", "CD8.Naive", "CD8.Prolif", "CD8.TE")
    merged_types[current_types %in% cd8_subtypes] <- "CD8 T"
    
    # B cells - merge all B cell subtypes
    b_subtypes <- c("B_exhausted", "B_immature", "B_naive", 
                    "B_non-switched_memory", "B_switched_memory", "B_malignant")
    merged_types[current_types %in% b_subtypes] <- "B cell"
    
    # Plasma cells - merge all plasma cell subtypes
    plasma_subtypes <- c("Plasma_cell_IgA", "Plasma_cell_IgG", "Plasma_cell_IgM")
    merged_types[current_types %in% plasma_subtypes] <- "Plasma cell"
    
    # HSPC - merge all hematopoietic stem/progenitor cells
    hspc_subtypes <- c("HSC_CD38neg", "HSC_CD38pos", "HSC_erythroid", 
                       "HSC_MK", "HSC_myeloid", "HSC_prolif")
    merged_types[current_types %in% hspc_subtypes] <- "HSPC"
    
    return(merged_types)
}

# -----------------------------------------------
