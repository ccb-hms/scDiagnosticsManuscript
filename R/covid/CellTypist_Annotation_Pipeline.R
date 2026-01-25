# ---------------------------------------------------
# COVID-19 PBMC - CellTypist Annotation Pipeline
# ---------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)

# Source files 
source("R/auxiliary/environmentSetupCellTypist.R")

# ---------------------------------------------------

# _________________________
# Pipeline Initialization
# _________________________

cat("=== COVID-19 CellTypist Annotation Pipeline ===\n")

# Create data output directory if it doesn't exist
if (!dir.exists("data/covid")) {
  dir.create("data/covid", recursive = TRUE)
  cat("✓ Created output directory: data/covid\n")
}

# Setting the seed
set.seed(0)

# __________________________________
# Step 1: Python Environment Setup
# __________________________________

cat("\n--- Step 1: Setting up Python environment ---\n")
environmentSetupCellTypist()

# _______________________________________
# Step 2: Data Preparation (R -> Python)
# _______________________________________

cat("\n--- Step 2: Preparing data for Python annotation ---\n")
source("R/covid/CellTypist_Data_Preparation.R")

# ___________________________________________
# Step 3: Run CellTypist Annotation (Python)
# ___________________________________________

cat("\n--- Step 3: Running CellTypist annotation in Python ---\n")
python_script <- "python/covid/CellTypist_Annotation.py"
system(paste("conda run -n celltypist_env python", python_script))
cat("✓ Python script execution complete.\n")

# ____________________________________
# Step 4: Integrate Results (Python -> R)
# ____________________________________

cat("\n--- Step 4: Integrating annotations with original R object ---\n") 
source("R/covid/CellTypist_Results_Integration.R")

# _____________________________________
# Step 5: Clean Up Intermediate Files
# _____________________________________

cat("\n--- Step 5: Cleaning up all intermediate .h5ad files ---\n")
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

# _________________________
# Step 6: Final Validation
# _________________________

cat("\n--- Step 6: Final Validation ---\n")
result <- readRDS("data/covid/covid_data_sce.rds")

cat("Final object dimensions: ", nrow(result), " genes x ", ncol(result), " cells\n")

# Show CellTypist annotation columns by searching for the prefix
annotation_cols <- grep("celltypist", colnames(colData(result)), value = TRUE)

if(length(annotation_cols) > 0) {
  cat("✓ Successfully found CellTypist annotations in final object:\n")
  for(col in annotation_cols) {
    cat("- ", col, "\n")
  }
} else {
  cat("⚠️  Warning: No CellTypist annotations found in the final object!\n")
}

# __________________
# Pipeline Complete
# __________________

cat("\n🎉🎉🎉 Pipeline Complete! 🎉🎉🎉\n")
cat("✓ CellTypist annotations are now integrated into: data/covid/covid_data_sce.rds\n")