# ------------------------------------------------------
# Perform Azimuth Annotation With 
# ------------------------------------------------------

# Load libraries
library(Azimuth)
library(SingleCellExperiment)
library(Seurat)

# ------------------------------------------------------
 
# Function to perform Azimuth annotation with custom reference
performAzimuthAnnotation <- function(
    sce_obj, 
    dataset_name, 
    reference_path = "pbmcref") {
    
    cat("=== Running Azimuth annotation for", dataset_name, "===\n")
    
    # Ensure counts assay exists (convert from logcounts if needed)
    if(!"counts" %in% assayNames(sce_obj)) {
        cat("Converting logcounts to counts assay for", dataset_name, "...\n")
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
        
        approx_counts[approx_counts < 0] <- 0
        approx_counts <- round(approx_counts)
        assay(sce_obj, "counts") <- approx_counts
    } else {
        cat("Counts assay already exists for", dataset_name, "\n")
    }
    
    # Create Seurat object MANUALLY
    cat("Creating Seurat object manually...\n")
    counts_data <- assay(sce_obj, "counts")
    logcounts_data <- logcounts(sce_obj)
    metadata <- as.data.frame(colData(sce_obj))
    
    seurat_obj <- CreateSeuratObject(
        counts = counts_data,
        meta.data = metadata,
        project = paste0(dataset_name, "_query"),
        assay = "RNA"
    )
    
    seurat_obj <- SetAssayData(seurat_obj, 
                              assay = "RNA",
                              layer = "data", 
                              new.data = logcounts_data)
    DefaultAssay(seurat_obj) <- "RNA"
    
    # Run Azimuth
    cat("Running Azimuth with reference:", reference_path, "...\n")
    seurat_annotated <- RunAzimuth(seurat_obj, reference = reference_path)
    
    # DEBUG: Check what columns actually exist and their types
    cat("Debugging column types after RunAzimuth...\n")
    available_cols <- colnames(seurat_annotated[[]])
    predicted_cols <- grep("predicted", available_cols, value = TRUE)
    cat("Available predicted columns:", paste(predicted_cols, collapse = ", "), "\n")
    
    # Check types of key columns
    for(col in predicted_cols) {
        val <- seurat_annotated[[col]]
        cat("Column", col, "- Class:", class(val), "- Type:", typeof(val), 
            "- Length:", length(val), "\n")
    }
    
    # SAFE EXTRACTION FUNCTION
    safe_extract <- function(obj, col_name, default = NULL) {
        if(col_name %in% colnames(obj[[]])) {
            val <- obj[[col_name, drop = TRUE]]
            cat("Extracting", col_name, "- Class:", class(val), "\n")
            
            # Handle different types safely
            if(is.list(val)) {
                # If it's a list, try to unlist or take first element
                if(length(val) == 1 && is.vector(val[[1]])) {
                    return(val[[1]])
                } else {
                    return(unlist(val))
                }
            } else {
                return(val)
            }
        } else {
            return(default)
        }
    }
    
    # SAFE NUMERIC CONVERSION
    safe_numeric <- function(x, default_val = 0.5) {
        if(is.null(x)) {
            return(rep(default_val, ncol(seurat_annotated)))
        }
        
        tryCatch({
            if(is.list(x)) {
                x <- unlist(x)
            }
            as.numeric(x)
        }, error = function(e) {
            cat("Warning: Could not convert to numeric, using default values\n")
            rep(default_val, ncol(seurat_annotated))
        })
    }
    
    # Extract results safely
    azimuth_results <- list(
        dataset_name = dataset_name,
        n_cells = ncol(sce_obj),
        annotation_date = Sys.Date(),
        azimuth_reference = reference_path
    )
    
    # Try different possible column names and extract safely
    if("predicted.celltype.l1" %in% available_cols) {
        azimuth_results$predicted.celltype.l1 <- safe_extract(seurat_annotated, "predicted.celltype.l1")
        azimuth_results$predicted.celltype.l1.score <- safe_extract(seurat_annotated, "predicted.celltype.l1.score")
    } else if("predicted.cell_type" %in% available_cols) {
        azimuth_results$predicted.celltype.l1 <- safe_extract(seurat_annotated, "predicted.cell_type")
        # Look for score column
        score_cols <- grep("score", predicted_cols, value = TRUE)
        if(length(score_cols) > 0) {
            azimuth_results$predicted.celltype.l1.score <- safe_extract(seurat_annotated, score_cols[1])
        }
    } else {
        # Use first predicted column
        if(length(predicted_cols) > 0) {
            azimuth_results$predicted.celltype.l1 <- safe_extract(seurat_annotated, predicted_cols[1])
            # Try to find corresponding score
            score_col <- paste0(predicted_cols[1], ".score")
            azimuth_results$predicted.celltype.l1.score <- safe_extract(seurat_annotated, score_col)
        }
    }
    
    # For L2, try similar approach
    if("predicted.celltype.l2" %in% available_cols) {
        azimuth_results$predicted.celltype.l2 <- safe_extract(seurat_annotated, "predicted.celltype.l2")
        azimuth_results$predicted.celltype.l2.score <- safe_extract(seurat_annotated, "predicted.celltype.l2.score")
    } else {
        # Use L1 as L2 if L2 doesn't exist
        azimuth_results$predicted.celltype.l2 <- azimuth_results$predicted.celltype.l1
        azimuth_results$predicted.celltype.l2.score <- azimuth_results$predicted.celltype.l1.score
    }
    
    # Extract mapping score
    azimuth_results$mapping.score <- safe_extract(seurat_annotated, "mapping.score", 
                                                  default = rep(0.5, ncol(seurat_annotated)))
    
    # Add UMAP coordinates
    if("ref.umap" %in% names(seurat_annotated@reductions)) {
        azimuth_results$umap.1 <- Embeddings(seurat_annotated, "ref.umap")[,1]
        azimuth_results$umap.2 <- Embeddings(seurat_annotated, "ref.umap")[,2]
    }
    
    # SAFE ADDITION TO SCE colData
    cat("Adding results to SCE colData...\n")
    
    tryCatch({
        colData(sce_obj)$azimuth_celltype_l1 <- as.character(azimuth_results$predicted.celltype.l1)
        colData(sce_obj)$azimuth_celltype_l2 <- as.character(azimuth_results$predicted.celltype.l2)
        colData(sce_obj)$azimuth_score_l1 <- safe_numeric(azimuth_results$predicted.celltype.l1.score)
        colData(sce_obj)$azimuth_score_l2 <- safe_numeric(azimuth_results$predicted.celltype.l2.score)
        colData(sce_obj)$azimuth_mapping_score <- safe_numeric(azimuth_results$mapping.score)
        
        # Add UMAP if available
        if(!is.null(azimuth_results$umap.1)) {
            colData(sce_obj)$azimuth_umap1 <- azimuth_results$umap.1
            colData(sce_obj)$azimuth_umap2 <- azimuth_results$umap.2
        }
        
        cat("✓ Successfully added all results to SCE\n")
        
    }, error = function(e) {
        cat("Error adding results to SCE:", e$message, "\n")
        stop("Failed to add results to SCE: ", e$message)
    })
    
    cat("✓ Azimuth annotation complete for", dataset_name, "\n")
    cat("L1 cell types:", length(unique(azimuth_results$predicted.celltype.l1)), "\n")
    cat("L2 cell types:", length(unique(azimuth_results$predicted.celltype.l2)), "\n")
    
    return(list(
        sce_annotated = sce_obj,
        azimuth_results = azimuth_results,
        seurat_object = seurat_annotated
    ))
}