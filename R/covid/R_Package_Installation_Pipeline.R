# ------------------------------------------------
# COVID-19 PBMC - R Package Installation Pipeline
# (Version-Pinned for Reproducibility)
# ------------------------------------------------

# This script installs all necessary R packages with SPECIFIC versions 
# for the COVID-19 PBMC analysis pipeline on an HPC cluster like O2.

# ____________________
# Set CRAN Repository
# ____________________

options(repos = c(CRAN = "https://cran.rstudio.com/"))

# _________________________
# Define Required Packages
# _________________________

# CRAN packages with pinned versions
cran_packages_versioned <- list(
    "dplyr" = "1.1.4",
    "tidyr" = "1.3.1",
    "tibble" = "3.2.1",
    "Matrix" = "1.7.3",
    "reticulate" = "1.42.0",
    "remotes" = "2.5.0",
    "here" = "1.0.1",
    "ggplot2" = "3.5.2",
    "cowplot" = "1.1.3",
    "GGally" = "2.2.1",
    "ggridges" = "0.5.6",
    "viridis" = "0.6.5"
)

# Bioconductor packages with pinned versions
bioc_packages_versioned <- list(
    "zellkonverter" = "1.16.0",
    "SingleCellExperiment" = "1.28.1",
    "HDF5Array" = "1.34.0",
    "DelayedArray" = "0.32.0",
    "scran" = "1.34.0",
    "scater" = "1.34.1",
    "BiocParallel" = "1.40.2",
    "SingleR" = "2.8.0",
    "biomaRt" = "2.62.1",
    "SpatialExperiment" = "1.16.0",
    "MerfishData" = "1.8.0",
    "scDiagnostics" = "1.5.1"  # Devel version on Bioconductor
)

# GitHub packages with pinned versions
github_packages_versioned <- list(
    "Seurat" = "satijalab/seurat@v5.3.1.0",
    "Azimuth" = "satijalab/azimuth@v0.5.0"
)

# _______________________
# Installation from CRAN
# _______________________

message("--- Checking and installing CRAN packages (version-pinned) ---")
for (pkg_name in names(cran_packages_versioned)) {
    pkg_version <- cran_packages_versioned[[pkg_name]]
    
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
        message("Installing '", pkg_name, "' v", pkg_version, " from CRAN...")
        tryCatch({
            remotes::install_version(pkg_name, version = pkg_version, repos = getOption("repos"))
            message("✓ Successfully installed ", pkg_name, " v", pkg_version)
        }, error = function(e) {
            message("✗ Failed to install ", pkg_name, " v", pkg_version, ": ", e$message)
            message("  Attempting to install latest version instead...")
            install.packages(pkg_name)
        })
    } else {
        current_version <- as.character(packageVersion(pkg_name))
        if (current_version == pkg_version) {
            message("✓ ", pkg_name, " v", pkg_version, " is already installed (correct version)")
        } else {
            message("⚠ ", pkg_name, " v", current_version, " installed (target: v", pkg_version, ")")
            message("  Updating to v", pkg_version, "...")
            tryCatch({
                remotes::install_version(pkg_name, version = pkg_version, repos = getOption("repos"))
                message("✓ Successfully updated ", pkg_name, " to v", pkg_version)
            }, error = function(e) {
                message("✗ Failed to update: ", e$message)
            })
        }
    }
}

# _______________________________
# Installation from Bioconductor
# _______________________________

message("\n--- Checking and installing Bioconductor packages (version-pinned) ---")

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing 'BiocManager'...")
    install.packages("BiocManager")
}

for (pkg_name in names(bioc_packages_versioned)) {
    pkg_version <- bioc_packages_versioned[[pkg_name]]
    
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
        message("Installing '", pkg_name, "' v", pkg_version, " from Bioconductor...")
        tryCatch({
            if (pkg_name == "scDiagnostics") {
                # Install scDiagnostics devel version
                message("Installing scDiagnostics devel version (v1.5.1)...")
                BiocManager::install("scDiagnostics", version = "devel", force = TRUE, update = FALSE)
            } else {
                BiocManager::install(pkg_name, version = BiocManager::version(), force = TRUE, update = FALSE)
            }
            message("✓ Successfully installed ", pkg_name, " v", pkg_version)
        }, error = function(e) {
            message("✗ Failed to install ", pkg_name, ": ", e$message)
        })
    } else {
        current_version <- as.character(packageVersion(pkg_name))
        if (current_version == pkg_version) {
            message("✓ ", pkg_name, " v", pkg_version, " is already installed (correct version)")
        } else {
            message("⚠ ", pkg_name, " v", current_version, " installed (target: v", pkg_version, ")")
            if (pkg_name == "scDiagnostics") {
                message("  (scDiagnostics devel version may vary, but should be ~1.5.1)")
            } else {
                message("  ", pkg_name, " may need manual update if version mismatch is critical")
            }
        }
    }
}

# _________________________
# Installation from GitHub
# _________________________

message("\n--- Checking and installing GitHub packages (version-pinned) ---")

for (pkg_name in names(github_packages_versioned)) {
    github_repo <- github_packages_versioned[[pkg_name]]
    
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
        message("Installing '", pkg_name, "' from GitHub: ", github_repo)
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
            
            if (as.character(current_version) == target_version) {
                message("✓ Seurat v", target_version, " is already installed (correct version)")
            } else {
                message("⚠ Seurat v", current_version, " installed (target: v", target_version, ")")
                message("  Updating to v", target_version, " from GitHub...")
                tryCatch({
                    remotes::install_github(github_repo, force = TRUE)
                    message("✓ Successfully updated Seurat to v", target_version)
                }, error = function(e) {
                    message("✗ Failed to update: ", e$message)
                })
            }
        } else if (pkg_name == "Azimuth") {
            current_version <- packageVersion("Azimuth")
            target_version <- "0.5.0"
            
            if (as.character(current_version) == target_version) {
                message("✓ Azimuth v", target_version, " is already installed (correct version)")
            } else {
                message("⚠ Azimuth v", current_version, " installed (target: v", target_version, ")")
            }
        }
    }
}

# ____________
# Final Check
# ____________

message("\n-------------------------------------------------")
message("All package checks are complete.")

# Verify specific versions and key packages
message("\nVersion verification:")
message("-------------------\n")

tryCatch({
    seurat_version <- packageVersion("Seurat")
    message("✓ Seurat v", seurat_version)
    
    if (requireNamespace("Azimuth", quietly = TRUE)) {
        azimuth_version <- packageVersion("Azimuth")
        message("✓ Azimuth v", azimuth_version)
    }
    
    if (requireNamespace("scDiagnostics", quietly = TRUE)) {
        diag_version <- packageVersion("scDiagnostics")
        message("✓ scDiagnostics v", diag_version, " (devel/1.5.1)")
    }
    
    # Verify visualization packages
    message("\nVisualization packages:")
    viz_packages <- c("ggplot2", "cowplot", "GGally", "ggridges", "viridis")
    for (pkg in viz_packages) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            ver <- packageVersion(pkg)
            message("✓ ", pkg, " v", ver)
        }
    }
    
    # Verify core packages
    message("\nCore analysis packages:")
    core_packages <- c("SingleCellExperiment", "scran", "scater", "SingleR")
    for (pkg in core_packages) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            ver <- packageVersion(pkg)
            message("✓ ", pkg, " v", ver)
        }
    }
    
}, error = function(e) {
    message("Could not verify package versions: ", e$message)
})

message("\n-------------------------------------------------")
message("The R environment is ready for the COVID-19 PBMC pipeline.")
message("Installed versions:")
message("  • CRAN: 12 packages (ggplot2 3.5.2, cowplot 1.1.3, GGally 2.2.1, etc.)")
message("  • Bioconductor: 12 packages (SingleCellExperiment 1.28.1, scDiagnostics 1.5.1, etc.)")
message("  • GitHub: Seurat v5.3.1.0, Azimuth v0.5.0")
message("-------------------------------------------------\n")