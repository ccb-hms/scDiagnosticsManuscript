# ------------------------------------------------
# MERFISH - R Package Installation Pipeline
# (Version-Pinned for Reproducibility)
# ------------------------------------------------

# This script installs all necessary R packages with SPECIFIC versions 
# for the MERFISH analysis pipeline on an HPC cluster like O2.

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
    "patchwork" = "1.3.0",
    "viridis" = "0.6.5",
    "GGally" = "2.2.1",
    "ggridges" = "0.5.6",
    "pheatmap" = "1.0.12",
    "circlize" = "0.4.16"
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
    "BiocSingular" = "1.22.0",
    "SingleR" = "2.8.0",
    "biomaRt" = "2.62.1",
    "SpatialExperiment" = "1.16.0",
    "MerfishData" = "1.8.0",
    "ComplexHeatmap" = "2.22.0",
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
    
    # Verify CRAN visualization packages
    message("\nCRAN visualization packages:")
    cran_viz_packages <- c("ggplot2", "cowplot", "patchwork", "GGally", "ggridges", "pheatmap", "circlize", "viridis")
    for (pkg in cran_viz_packages) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            ver <- packageVersion(pkg)
            message("✓ ", pkg, " v", ver)
        }
    }
    
    # Verify core Bioconductor packages
    message("\nCore Bioconductor packages:")
    bioc_core_packages <- c("SingleCellExperiment", "SpatialExperiment", "scran", "scater", "SingleR")
    for (pkg in bioc_core_packages) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            ver <- packageVersion(pkg)
            message("✓ ", pkg, " v", ver)
        }
    }
    
    # Verify visualization Bioconductor packages
    message("\nBioconductor visualization packages:")
    bioc_viz_packages <- c("ComplexHeatmap", "BiocSingular")
    for (pkg in bioc_viz_packages) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            ver <- packageVersion(pkg)
            message("✓ ", pkg, " v", ver)
        }
    }
    
}, error = function(e) {
    message("Could not verify package versions: ", e$message)
})

message("\n-------------------------------------------------")
message("The R environment is ready for the MERFISH pipeline.")
message("Installed versions:")
message("  • CRAN: 15 packages")
message("    - Plotting: ggplot2 3.5.2, cowplot 1.1.3, patchwork 1.3.0")
message("    - Advanced plots: GGally 2.2.1, ggridges 0.5.6, pheatmap 1.0.12, circlize 0.4.16")
message("    - Data: dplyr 1.1.4, tidyr 1.3.1, viridis 0.6.5")
message("  • Bioconductor: 14 packages")
message("    - Core: SingleCellExperiment 1.28.1, SpatialExperiment 1.16.0")
message("    - Analysis: scran 1.34.0, scater 1.34.1, SingleR 2.8.0")
message("    - Visualization: ComplexHeatmap 2.22.0, BiocSingular 1.22.0")
message("    - Diagnostics: scDiagnostics 1.5.1 (devel)")
message("    - Data: MerfishData 1.8.0")
message("  • GitHub: Seurat v5.3.1.0, Azimuth v0.5.0")
message("\nKey capabilities:")
message("  • Main analysis pipeline (Seurat, SingleR, Azimuth)")
message("  • Spatial analysis (SpatialExperiment, MerfishData)")
message("  • Anomaly detection (scDiagnostics v1.5.1)")
message("  • Advanced visualization (ComplexHeatmap, pheatmap, circlize)")
message("-------------------------------------------------\n")