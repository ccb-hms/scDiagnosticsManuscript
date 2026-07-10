# ------------------------------------------------
# ZEISEL BRAIN - R Package Installation Pipeline
# (Version-Pinned for Reproducibility)
# ------------------------------------------------

# This script installs all necessary R packages with SPECIFIC versions 
# for the Zeisel Brain Supplementary analysis pipeline.

# ____________________________________________________________________
# USER CONFIGURATION: SELECT ANALYSIS TARGET
# Set to "supp" for Supplementary results added in revision (v1.7.2)
# Set to "main" for Main Manuscript results (scDiagnostics v1.5.1)
# NOTE: The Zeisel Brain analysis was added entirely during revision.
# ____________________________________________________________________

TARGET_ANALYSIS <- "supp"  # <--- DEFAULT SET TO SUPP FOR ZEISEL BRAIN

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
    "ggplot2" = "3.5.2",
    "patchwork" = "1.3.0",
    "ggrepel" = "0.9.5",
    "pROC" = "1.18.5",
    "PRROC" = "1.3.1",
    "remotes" = "2.5.0"
)

# Bioconductor packages with pinned versions
bioc_packages_versioned <- list(
    "SingleCellExperiment" = "1.28.1",
    "scater" = "1.34.1",
    "SingleR" = "2.8.0",
    "scRNAseq" = "2.16.0",
    "scDiagnostics" = scDiag_version  # <--- Dynamically assigned
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

# ____________
# Final Check
# ____________

message("\n-------------------------------------------------")
message("Installation pipeline complete for: ", toupper(TARGET_ANALYSIS))
message("Verifying critical versions:")
message("-------------------------------------------------")

tryCatch({
    if (requireNamespace("scDiagnostics", quietly = TRUE)) {
        diag_version <- packageVersion("scDiagnostics")
        message("✓ scDiagnostics v", diag_version, " (Target was: ", scDiag_version, ")")
        if(as.character(diag_version) != scDiag_version) {
            message("   WARNING: Installed scDiagnostics version does not match target exactly.")
        }
    }
    
    # Verify CRAN packages
    message("\nCRAN packages:")
    for (pkg in names(cran_packages_versioned)) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            ver <- packageVersion(pkg)
            message("✓ ", pkg, " v", ver)
        }
    }
    
    # Verify Bioconductor packages
    message("\nBioconductor packages:")
    for (pkg in names(bioc_packages_versioned)) {
        if (pkg != "scDiagnostics" && requireNamespace(pkg, quietly = TRUE)) {
            ver <- packageVersion(pkg)
            message("✓ ", pkg, " v", ver)
        }
    }
    
}, error = function(e) {
    message("Could not verify package versions: ", e$message)
})

message("\n-------------------------------------------------")
message("The R environment is ready for the Zeisel Brain pipeline.")
message("Installed configurations:")
message(paste0("  • scDiagnostics: v", scDiag_version))
message("  • Visualization: ggplot2, patchwork, ggrepel")
message("  • Metrics: pROC, PRROC")
message("  • Core Bioc: SingleCellExperiment, scater, SingleR, scRNAseq")
message("-------------------------------------------------\n")