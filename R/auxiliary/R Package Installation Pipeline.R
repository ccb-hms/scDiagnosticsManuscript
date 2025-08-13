# ---------------------------------
# R Package Installation Pipeline
# ---------------------------------

# This script checks for and installs all necessary R packages 
# for the analysis pipeline on an HPC cluster like O2. It avoids 
# re-installing packages that are already present.

# __________________________________
# Set CRAN Repository
# __________________________________

# Using a reliable CRAN mirror is good practice on HPC systems.
options(repos = c(CRAN = "https://cran.rstudio.com/"))


# __________________________________
# Define Required Packages
# __________________________________

# Packages to be installed from CRAN
cran_packages <- c(
    "dplyr",        # Data manipulation
    "tidyr",        # Data tidying
    "tibble",       # Modern data frames
    "Matrix",       # Handling sparse matrices (often base, but good to check)
    "reticulate",   # R interface to Python
    "remotes",      # Needed to install packages from GitHub
    "here"          # For robust file paths
    # NOTE: Seurat moved to GitHub section for version 5.3.1.0
)

# Packages to be installed from Bioconductor
bioc_packages <- c(
    "zellkonverter",      # Reading h5ad files
    "SingleCellExperiment", # Core object for scRNA-seq data
    "HDF5Array",          # On-disk data representation
    "DelayedArray",       # Handling large, on-disk arrays
    "scran",              # Single-cell analysis tools (HVGs, etc.)
    "scater",             # Single-cell analysis tools (PCA, UMAP, etc.)
    "BiocParallel",       # For parallel processing in Bioconductor
    "SingleR",            # A reference-based annotation tool
    "biomaRt"             # For converting gene IDs
)

# Packages to be installed from GitHub
github_packages <- list(
    "Seurat" = "satijalab/seurat@v5.3.1.0",  # Specific version from GitHub
    "Azimuth" = "satijalab/azimuth"           # Latest from GitHub
)


# __________________________________
# Installation from CRAN
# __________________________________

message("--- Checking and installing CRAN packages ---")
for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Installing '", pkg, "' from CRAN...")
        install.packages(pkg)
    } else {
        message("Package '", pkg, "' is already installed.")
    }
}


# __________________________________
# Installation from Bioconductor
# __________________________________

message("\n--- Checking and installing Bioconductor packages ---")

# First, ensure BiocManager itself is installed
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


# __________________________________
# Installation from GitHub
# __________________________________

message("\n--- Checking and installing GitHub packages ---")

# Install packages from GitHub with specific versions
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
        # Check if we need to update to the specific version
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


# __________________
# Final Check
# __________________

message("\n-------------------------------------------------")
message("All package checks are complete.")

# Verify specific versions
message("\nVersion verification:")
tryCatch({
    seurat_version <- packageVersion("Seurat")
    message("✓ Seurat version: ", seurat_version)
    
    if (requireNamespace("Azimuth", quietly = TRUE)) {
        azimuth_version <- packageVersion("Azimuth")
        message("✓ Azimuth version: ", azimuth_version)
    }
}, error = function(e) {
    message("Could not verify package versions: ", e$message)
})

message("The R environment should now be ready for the pipeline.")
message("-------------------------------------------------")