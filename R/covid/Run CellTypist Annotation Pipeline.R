# ----------------------------------------------
# COVID-19 - CellTypist Annotation Pipeline
# ----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)

# Source files 
source("R/auxiliary/environmentSetupCellTypist.R")

# ----------------------------------------------

cat("=== COVID-19 CellTypist Annotation Pipeline ===\n\n")

# Create data directory if it doesn't exist
if (!dir.exists("data/covid")) {
  dir.create("data/covid", recursive = TRUE)
  cat("Created data/covid directory\n")
}

# Step 0: Setup Python environment
cat("Step 0: Setting up Python environment...\n")
environmentSetupCellTypist()

# Step 1: Data preparation
cat("\nStep 1: Preparing data for Python annotation...\n")
source("R/covid/Python Annotation Data Preparation.R")

# Step 2: Python annotation  
cat("\nStep 2: Running CellTypist annotation...\n")
python_script <- "python/covid/CellTypist_Annotation.py"
system(paste("conda run -n celltypist_env python", python_script))

# Step 3: Integration (adds annotations to original covid_data.rds)
cat("\nStep 3: Integrating annotations with original COVID data...\n") 
source("R/covid/CellTypist Results Integration.R")

# Step 4: Clean up ALL h5ad files
cat("\nStep 4: Cleaning up all intermediate h5ad files...\n")
h5ad_files <- c(
  "data/covid/reference_data.h5ad",
  "data/covid/query_data.h5ad", 
  "data/covid/annotated_query_data.h5ad"
)

for (file in h5ad_files) {
  if (file.exists(file)) {
    file.remove(file)
    cat("✓ Removed:", file, "\n")
  }
}

cat("\n=== Pipeline Complete! ===\n")
cat("✓ CellTypist annotations integrated into: data/covid/covid_data.rds\n")

# Final validation
cat("\nFinal validation:\n")
result <- readRDS("data/covid/covid_data.rds")
cat("Final dimensions:", dim(result), "\n")
cat("Total colData columns:", ncol(colData(result)), "\n")

# Show CellTypist annotation columns
annotation_cols <- grep("predicted|majority|conf", 
                        colnames(colData(result)), value = TRUE)
if(length(annotation_cols) > 0) {
  cat("CellTypist annotations added:", paste(annotation_cols, collapse = ", "), "\n")
} else {
  cat("⚠️  Warning: No CellTypist annotations found!\n")
}

# Show what's left in the directory
cat("\nFinal files in data/covid/:\n")
covid_files <- list.files("data/covid/", full.names = TRUE)
for(file in covid_files) {
  if(file.info(file)$isdir == FALSE) {
    size_mb <- round(file.size(file) / 1024^2, 1)
    cat("- ", basename(file), " (", size_mb, " MB)\n", sep = "")
  }
}

cat("\n🎉 Your COVID data now has CellTypist annotations integrated!\n")