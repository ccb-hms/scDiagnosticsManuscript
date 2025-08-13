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
    "Seurat"        # Required by Azimuth
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
github_repo <- "satijalab/azimuth"
github_pkg_name <- "Azimuth"


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

# Azimuth requires installation from GitHub via the 'remotes' package
if (!requireNamespace(github_pkg_name, quietly = TRUE)) {
    message("Installing '", github_pkg_name, "' from GitHub repository: ", github_repo)
    # The 'remotes' package should have been installed in the CRAN step
    remotes::install_github(github_repo)
} else {
    message("Package '", github_pkg_name, "' is already installed.")
}


# __________________
# Final Check
# __________________

message("\n-------------------------------------------------")
message("All package checks are complete.")
message("The R environment should now be ready for the pipeline.")
message("-------------------------------------------------")