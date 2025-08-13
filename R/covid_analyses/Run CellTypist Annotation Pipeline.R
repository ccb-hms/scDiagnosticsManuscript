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
if (!dir.exists("data")) {
  dir.create("data")
  cat("Created data directory\n")
}

# Step 0: Setup Python environment (one-time setup)
cat("Step 0: Setting up Python environment...\n")
environmentSetupCellTypist()

# Step 1: Data preparation (minimal colData approach)
cat("\nStep 1: Preparing data for Python annotation...\n")
source("COVID-19 - Python Annotation Data Preparation.R")

# Step 2: Python annotation  
cat("\nStep 2: Running CellTypist annotation...\n")
system("COVID-19 - CellTypist Cell Type Annotation.py")

# Step 3: Merge results with original data
cat("\nStep 3: Integrating annotations with original COVID data...\n") 
source("COVID-19 - CellTypist Results Integration.R")

cat("\n=== Pipeline Complete! ===\n")
cat("Your annotated data is saved as: data/covid_data_annotated.rds\n")

# Quick validation
cat("\nQuick validation:\n")
result <- readRDS("data/covid_data_annotated.rds")
cat("Final dimensions:", dim(result), "\n")
cat("Original colData preserved:", 
    ncol(colData(result)) > 3, "\n")  # Assuming original had >3 columns
cat("CellTypist annotations added:", 
    any(grepl("predicted|majority|conf", colnames(colData(result)))), "\n")

# Show annotation columns
annotation_cols <- grep("predicted|majority|conf", 
                        colnames(colData(result)), value = TRUE)
if(length(annotation_cols) > 0) {
  cat("Available annotations:", paste(annotation_cols, collapse = ", "), "\n")
}