# -------------------------------------------------------------------
# Add Reference SCE Object Using Union HVGs (Reference & Query HVGs)
# -------------------------------------------------------------------

# Required libraries
library(scran)
library(scater)
library(SingleCellExperiment)

# -------------------------------------------------------------------

# Union HVGs for PCA
addReferencePCA <- function(ref_sce, query_sce, ref_name) {
    cat("Computing diagnostic-optimized PCA for", ref_name, "\n")
    
    # Get HVGs from both datasets
    ref_dec <- suppressWarnings(modelGeneVar(ref_sce))
    ref_hvgs <- getTopHVGs(ref_dec, n = 2000)
    
    query_dec <- suppressWarnings(modelGeneVar(query_sce))  
    query_hvgs <- getTopHVGs(query_dec, n = 2000)
    
    # Union of HVGs (captures both reference and query variation)
    union_hvgs <- union(ref_hvgs, query_hvgs)
    
    # Keep only genes present in both datasets
    common_genes <- intersect(rownames(ref_sce), rownames(query_sce))
    diagnostic_hvgs <- intersect(union_hvgs, common_genes)
    
    cat("Reference HVGs:", length(ref_hvgs), "\n")
    cat("Query HVGs:", length(query_hvgs), "\n") 
    cat("Diagnostic HVGs (union):", length(diagnostic_hvgs), "\n")
    
    # Compute PCA on reference using union HVGs
    ref_sce <- runPCA(ref_sce, subset_row = diagnostic_hvgs, ncomponents = 50)
    
    # Store HVG info for reference
    metadata(ref_sce)$ref_hvgs <- ref_hvgs
    metadata(ref_sce)$query_hvgs <- query_hvgs  
    metadata(ref_sce)$diagnostic_hvgs <- diagnostic_hvgs
    
    return(ref_sce)
}
