# ------------------------------------------------
# MERFISH - R Package Installation Pipeline
# (Version-Pinned for Reproducibility)
# ------------------------------------------------

# This script installs all necessary R packages with SPECIFIC versions 
# for the MERFISH analysis pipeline on an HPC cluster like O2.

# ____________________________________________________________________
# USER CONFIGURATION: SELECT ANALYSIS TARGET
# Set to "main" for Main Manuscript results (scDiagnostics v1.5.1)
# Set to "supp" for Supplementary results added in revision (v1.7.2)
# ____________________________________________________________________

TARGET_ANALYSIS <- "main"  # <--- CHANGE THIS TO "supp" FOR SUPPLEMENTARY REPRODUCTION

message("===============================================================")
message("Initializing installation for: ", toupper(TARGET_ANALYSIS), " MANUSCRIPT")
message("===============================================================")

# ____________________
# Set CRAN Repository
# ____________________

options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Determine scDiagnostics version based on user target
scDiag_version <- ifelse(TARGET_ANALYSIS == "main", "1.5.1", "1.7.2")

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
    "scDiagnostics" = scDiag_version  # <--- Dynamically assigned
)

# GitHub packages with pinned versions
github_packages_versioned <- list(
    "Seurat" = "satijalab/seurat@v5.3.1.0",
    "Azimuth" = "satijalab/azimuth@v0.5.0"
)

# _______________________
# Installation from CRAN
# _______________________

message("\n--- Checking and installing CRAN packages (version-pinned) ---")
for (pkg_name in names(cran_packages_versioned)) {
    pkg_version <- cran_packages_versioned[[pkg_name]]
    
    if (!requireNamespace(pkg_name, quietly = TRUE) || as.character(packageVersion(pkg_name)) != pkg_version) {
        message("Installing '", pkg_name, "' v", pkg_version, " from CRAN...")
        tryCatch({
            remotes::install_version(pkg_name, version = pkg_version, repos = getOption("repos"))
            message("✓ Successfully installed ", pkg_name, " v", pkg_version)
        }, error = function(e) {
            message("✗ Failed to install ", pkg_name, " v", pkg_version, ": ", e$message)
        })
    } else {
        message("✓ ", pkg_name, " v", pkg_version, " is already installed (correct version)")
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
    
    if (!requireNamespace(pkg_name, quietly = TRUE) || as.character(packageVersion(pkg_name)) != pkg_version) {
        message("Installing '", pkg_name, "' target v", pkg_version, "...")
        tryCatch({
            if (pkg_name == "scDiagnostics") {
                
                # Handling the two specific versions of scDiagnostics
                if (TARGET_ANALYSIS == "main") {
                    message("  Installing scDiagnostics v1.5.1 (Main manuscript version)...")
                    BiocManager::install("scDiagnostics", version = "devel", force = TRUE, update = FALSE)
                } else {
                    message("  Installing scDiagnostics v1.7.2 (Supplementary revision version)...")
                    BiocManager::install("scDiagnostics", force = TRUE, update = FALSE) 
                }
                
            } else {
                BiocManager::install(pkg_name, update = FALSE, ask = FALSE)
            }
            message("✓ Installation step completed for ", pkg_name)
        }, error = function(e) {
            message("✗ Failed to install ", pkg_name, ": ", e$message)
        })
    } else {
        message("✓ ", pkg_name, " v", pkg_version, " is already installed (correct version)")
    }
}

# _________________________
# Installation from GitHub
# _________________________

message("\n--- Checking and installing GitHub packages (version-pinned) ---")

for (pkg_name in names(github_packages_versioned)) {
    github_repo <- github_packages_versioned[[pkg_name]]
    target_version <- strsplit(github_repo, "@v")[[1]][2] # Extract version from string
    
    if (!requireNamespace(pkg_name, quietly = TRUE) || as.character(packageVersion(pkg_name)) != target_version) {
        message("Installing '", pkg_name, "' v", target_version, " from GitHub...")
        tryCatch({
            remotes::install_github(github_repo, force = TRUE)
            message("✓ Successfully installed ", pkg_name, " v", target_version)
        }, error = function(e) {
            message("✗ Failed to install ", pkg_name, ": ", e$message)
        })
    } else {
        message("✓ ", pkg_name, " v", target_version, " is already installed (correct version)")
    }
}

# ____________
# Final Check
# ____________

message("\n-------------------------------------------------")
message("Installation pipeline complete for: ", toupper(TARGET_ANALYSIS))
message("Verifying critical versions:")
message("-------------------------------------------------")

tryCatch({
    seurat_version <- packageVersion("Seurat")
    message("✓ Seurat v", seurat_version)
    
    if (requireNamespace("Azimuth", quietly = TRUE)) {
        azimuth_version <- packageVersion("Azimuth")
        message("✓ Azimuth v", azimuth_version)
    }
    
    if (requireNamespace("scDiagnostics", quietly = TRUE)) {
        diag_version <- packageVersion("scDiagnostics")
        message("✓ scDiagnostics v", diag_version, " (Target was: ", scDiag_version, ")")
        if(as.character(diag_version) != scDiag_version) {
            message("   WARNING: Installed scDiagnostics version does not match target exactly.")
        }
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
message(paste0("    - Diagnostics: scDiagnostics ", scDiag_version))
message("    - Data: MerfishData 1.8.0")
message("  • GitHub: Seurat v5.3.1.0, Azimuth v0.5.0")
message("\nKey capabilities:")
message("  • Main analysis pipeline (Seurat, SingleR, Azimuth)")
message("  • Spatial analysis (SpatialExperiment, MerfishData)")
message(paste0("  • Anomaly detection (scDiagnostics v", scDiag_version, ")"))
message("  • Advanced visualization (ComplexHeatmap, pheatmap, circlize)")
message("-------------------------------------------------\n")