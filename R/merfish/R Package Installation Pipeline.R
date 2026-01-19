# ------------------------------------------------
# MERFISH - R Package Installation Pipeline
# ------------------------------------------------

# This script checks for and installs all necessary R packages 
# for the MERFISH analysis pipeline on an HPC cluster like O2. 

# ____________________
# Set CRAN Repository
# ____________________

options(repos = c(CRAN = "https://cran.rstudio.com/"))

# _________________________
# Define Required Packages
# _________________________

# Packages to be installed from CRAN
cran_packages <- c(
    "dplyr",        # Data manipulation
    "tidyr",        # Data tidying
    "tibble",       # Modern data frames
    "Matrix",       # Handling sparse matrices
    "reticulate",   # R interface to Python
    "remotes",      # Needed to install packages from GitHub
    "here",         # For robust file paths
    "ggplot2",      # Core plotting (main + supp figures)
    "cowplot",      # Combining plots into grids (main + supp figures)
    "patchwork",    # Combining plots (main figures)
    "viridis",      # Color palettes (ECM scores, supplementary figures)
    "GGally",       # ggpairs for PCA comparisons (supplementary Fig S4)
    "ggridges",     # Ridge plots for distributions (supplementary Fig S4)
    "pheatmap",     # Simple heatmaps (supplementary Fig S3)
    "circlize"      # Color mapping functions (supplementary Fig S6)
)

# Packages to be installed from Bioconductor
bioc_packages <- c(
    "zellkonverter",          # Reading/Writing h5ad files
    "SingleCellExperiment",   # Core object for scRNA-seq data
    "HDF5Array",              # On-disk data representation
    "DelayedArray",           # Handling large, on-disk arrays
    "scran",                  # Single-cell analysis tools (supplementary Fig S6)
    "scater",                 # Single-cell analysis tools (PCA, UMAP) (supplementary Fig S5)
    "BiocParallel",           # Parallel processing
    "BiocSingular",           # Singular Value Decomposition (supplementary Fig S5)
    "SingleR",                # Reference-based annotation
    "biomaRt",                # Gene ID conversion
    "SpatialExperiment",      # Core object for Spatial data
    "MerfishData",            # MERFISH specific data structures
    "ComplexHeatmap"          # Advanced heatmap visualization (supplementary Fig S6)
)

# Packages to be installed from GitHub
github_packages <- list(
    "Seurat" = "satijalab/seurat@v5.3.1.0",     # Specific version required for Azimuth
    "Azimuth" = "satijalab/azimuth",            # Reference mapping
    
    # --- Anomaly Detection & Diagnostics ---
    "scDiagnostics" = "ccb-hms/scDiagnostics"   # Anomaly detection (main + supplementary)
)

# _______________________
# Installation from CRAN
# _______________________

message("--- Checking and installing CRAN packages ---")
for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing '", pkg, "' from CRAN...")
        install.packages(pkg)
    } else {
        message("Package '", pkg, "' is already installed.")
    }
}

# _______________________________
# Installation from Bioconductor
# _______________________________

message("\n--- Checking and installing Bioconductor packages ---")

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing 'BiocManager'...")
    install.packages("BiocManager")
}

for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing '", pkg, "' from Bioconductor...")
        BiocManager::install(pkg)
    } else {
        message("Package '", pkg, "' is already installed.")
    }
}

# _________________________
# Installation from GitHub
# _________________________

message("\n--- Checking and installing GitHub packages ---")

for (pkg_name in names(github_packages)) {
    github_repo <- github_packages[[pkg_name]]
    
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
        message("Installing '", pkg_name, "' from GitHub repository: ", github_repo)
        tryCatch({
            remotes::install_github(github_repo, force = TRUE)
            message("✓ Successfully installed ", pkg_name)
        }, error = function(e) {
            message("✗ Failed to install ", pkg_name, ": ", e$message)
        })
    } else {
        # Specific check for Seurat version pinning
        if (pkg_name == "Seurat") {
            current_version <- packageVersion("Seurat")
            target_version <- "5.3.1.0"
            
            if (as.character(current_version) != target_version) {
                message("Updating Seurat from version ", current_version, " to ", target_version)
                remotes::install_github(github_repo, force = TRUE)
            } else {
                message("Seurat version ", target_version, " is already installed.")
            }
        } else {
            message("Package '", pkg_name, "' is already installed.")
        }
    }
}

# ____________
# Final Check
# ____________

message("\n-------------------------------------------------")
message("All package checks are complete.")

# Verify specific critical versions and key packages
message("\nVersion verification:")
tryCatch({
    seurat_version <- packageVersion("Seurat")
    message("✓ Seurat version: ", seurat_version)
    
    if (requireNamespace("Azimuth", quietly = TRUE)) {
        azimuth_version <- packageVersion("Azimuth")
        message("✓ Azimuth version: ", azimuth_version)
    }
    
    if (requireNamespace("scDiagnostics", quietly = TRUE)) {
        diag_version <- packageVersion("scDiagnostics")
        message("✓ scDiagnostics version: ", diag_version)
    }
    
    # Verify visualization packages for supplementary figures
    viz_packages <- c("ggplot2", "cowplot", "GGally", "ggridges", "pheatmap", "ComplexHeatmap", "circlize")
    for (pkg in viz_packages) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            message("✓ ", pkg, " is installed")
        }
    }
    
}, error = function(e) {
    message("Could not verify some package versions: ", e$message)
})

message("\nThe R environment should now be ready for the MERFISH pipeline.")
message("Packages installed for:")
message("  • Main analysis pipeline (Seurat, SingleR, Azimuth)")
message("  • Spatial analysis (SpatialExperiment, MerfishData)")
message("  • Quality diagnostics (scDiagnostics)")
message("  • Supplementary figures visualization:")
message("    - PCA comparisons (GGally, ggridges)")
message("    - Heatmaps (pheatmap, ComplexHeatmap, circlize)")
message("    - Plot composition (cowplot, patchwork)")
message("    - Color schemes (viridis)")
message("-------------------------------------------------")