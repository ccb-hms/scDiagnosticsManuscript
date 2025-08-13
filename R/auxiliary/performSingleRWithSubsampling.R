# -------------------------------------------------------
# Perform SingleR Cell Type Annotation With Subsampling
# -------------------------------------------------------

# Required libraries
library(SingleR)
library(SingleCellExperiment)
library(BiocParallel)
library(dplyr)
library(tibble)

# -------------------------------------------------------

# Combined function for subsampling (with rare cell type preservation) and SingleR annotation
performSingleRWithSubsampling <- function(
    ref_sce,
    query_sce, 
    ref_name,
    annotation_col = "cell_type",
    max_cells_ref = NULL,  # NULL means no subsampling
    rare_threshold = 1000, 
    bpparam = SerialParam()) {
    
    cat("=== SingleR Annotation:", ref_name, "===\n")
    
    # Check if annotation column exists
    if(!annotation_col %in% colnames(colData(ref_sce))) {
        stop("Annotation column '", annotation_col, "' not found in reference data colData")
    }
    
    # If max_cells_ref is NULL or >= total cells, use all cells (no subsampling)
    if(is.null(max_cells_ref) || max_cells_ref >= ncol(ref_sce)) {
        if(is.null(max_cells_ref)) {
            cat("Using all", ncol(ref_sce), "cells from", ref_name, "(no subsampling - max_cells_ref is NULL)\n")
        } else {
            cat("Using all", ncol(ref_sce), "cells from", ref_name, "(max_cells_ref >= total cells)\n")
        }
        
        ref_sce_sub <- ref_sce
        
        subsampling_info <- list(
            original_cells = ncol(ref_sce),
            final_cells = ncol(ref_sce),
            rare_types_preserved = NULL,
            rare_threshold = rare_threshold,
            subsampling_performed = FALSE,
            max_cells_ref = max_cells_ref
        )
        
    } else {
        # Smart subsampling with rare cell type preservation
        cat("Smart subsampling", ref_name, "from", ncol(ref_sce), "to", max_cells_ref, "cells\n")
        
        set.seed(42)  # For reproducibility
        
        # Get cell metadata
        cell_data <- colData(ref_sce) |>
            as.data.frame() |>
            tibble::rownames_to_column("cell_id")
        
        # Count cell types and identify rare vs common
        cell_type_counts <- table(cell_data[[annotation_col]])
        rare_types <- names(cell_type_counts)[cell_type_counts <= rare_threshold]
        common_types <- names(cell_type_counts)[cell_type_counts > rare_threshold]
        
        cat("Found", length(rare_types), "rare cell types (≤", rare_threshold, "cells):", 
            paste(rare_types, collapse = ", "), "\n")
        cat("Found", length(common_types), "common cell types (>", rare_threshold, "cells)\n")
        
        sampled_cells <- c()
        
        # Step 1: Keep ALL cells from rare types
        for(ct in rare_types) {
            ct_cells <- cell_data[cell_data[[annotation_col]] == ct, "cell_id"]
            sampled_cells <- c(sampled_cells, ct_cells)
            cat("Preserving all", length(ct_cells), "cells from rare type:", ct, "\n")
        }
        
        rare_count <- length(sampled_cells)
        remaining_budget <- max_cells_ref - rare_count
        
        cat("Rare cells preserved:", rare_count, "| Remaining budget:", remaining_budget, "\n")
        
        # Step 2: Distribute remaining budget among common types
        if(remaining_budget > 0 && length(common_types) > 0) {
            # Proportional sampling based on original cell type sizes
            common_type_counts <- cell_type_counts[common_types]
            total_common_cells <- sum(common_type_counts)
            
            # Calculate proportional allocation
            proportions <- common_type_counts / total_common_cells
            cells_per_common_type <- round(proportions * remaining_budget)
            
            # Handle rounding discrepancies
            total_allocated <- sum(cells_per_common_type)
            if(total_allocated != remaining_budget) {
                # Adjust the largest type by the difference
                largest_type <- names(which.max(cells_per_common_type))
                cells_per_common_type[largest_type] <- cells_per_common_type[largest_type] + 
                    (remaining_budget - total_allocated)
            }
            
            # Sample from each common type
            for(ct in common_types) {
                ct_cells <- cell_data[cell_data[[annotation_col]] == ct, "cell_id"]
                n_sample <- min(cells_per_common_type[ct], length(ct_cells))
                
                if(n_sample > 0) {
                    sampled_ct <- sample(ct_cells, n_sample)
                    sampled_cells <- c(sampled_cells, sampled_ct)
                    cat("Sampled", n_sample, "/", length(ct_cells), "cells from common type:", ct, "\n")
                }
            }
        } else if(remaining_budget <= 0) {
            cat("Warning: Rare cell types exceed budget! Using", rare_count, "cells total.\n")
        }
        
        # Create subsampled reference
        ref_sce_sub <- ref_sce[, sampled_cells]
        
        # Verify final composition
        final_counts <- table(colData(ref_sce_sub)[[annotation_col]])
        cat("Final reference composition:\n")
        for(ct in names(final_counts)) {
            original_count <- cell_type_counts[ct]
            preserved_pct <- round(100 * final_counts[ct] / original_count, 1)
            cat("  ", ct, ":", final_counts[ct], "/", original_count, 
                "(", preserved_pct, "%)\n")
        }
        
        subsampling_info <- list(
            original_cells = ncol(ref_sce),
            final_cells = ncol(ref_sce_sub),
            rare_types_preserved = rare_types,
            rare_threshold = rare_threshold,
            subsampling_performed = TRUE,
            max_cells_ref = max_cells_ref
        )
    }
    
    cat("Final reference size:", ncol(ref_sce_sub), "cells with", 
        length(unique(colData(ref_sce_sub)[[annotation_col]])), "cell types\n")
    
    # Perform SingleR annotation
    cat("Running SingleR annotation:", ref_name, "-> query\n")
    
    pred <- SingleR(
        test = logcounts(query_sce),
        ref = logcounts(ref_sce_sub), 
        labels = colData(ref_sce_sub)[[annotation_col]],
        de.method = "wilcox",
        de.n = 50,
        BPPARAM = bpparam            
    )
    
    cat("Completed annotation with", ref_name, "\n")
    cat("Predicted", length(unique(pred$labels)), "unique cell types\n\n")
    
    # Return both annotations and scores, plus subsampling info
    return(list(
        annotations = pred$labels,
        scores = pred$scores,
        subsampling_info = subsampling_info
    ))
}