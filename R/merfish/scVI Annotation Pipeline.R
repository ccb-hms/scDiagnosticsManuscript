# ----------------------------------------------------------
# MERFISH Mouse Colon IBD - Run Full scVI/scArches Pipeline
# ----------------------------------------------------------

# ___________________________________
# Step 1: Prepare Data (R -> Python)
# ___________________________________

message("\n--- Step 1: Preparing data for scVI (R -> Python) ---")
# Ensure you save the script below as "scVI_Data_Preparation.R" in R/merfish/
source("R/merfish/scVI Data Preparation.R")

# ________________________________
# Step 2: Run Annotation (Python)
# ________________________________

message("\n--- Step 2: Running scVI/scArches in Python ---")
# This command assumes a conda environment named 'scvi-env' has been set up
python_script <- "python/merfish/scVI_Annotation.py"
system(paste("conda run -n scvi-env python", python_script))
message("--- Python script finished ---")

# ________________________________________
# Step 3: Integrate Results (Python -> R)
# ________________________________________

message("\n--- Step 3: Integrating results back into R (Python -> R) ---")
# Ensure you save the integration script in R/merfish/
source("R/merfish/scVI Results Integration.R")

# ____________________________________
# Step 4: Clean Up Intermediate Files
# ____________________________________

message("\n--- Step 4: Cleaning up all intermediate files ---")

files_to_remove <- c(
  "data/merfish/scvi_reference_data.h5ad",
  "data/merfish/scvi_query_data.h5ad",
  "data/merfish/scvi_predictions.csv",    
  "data/merfish/scvi_reference_umap.csv" 
)

for (file in files_to_remove) {
  if (file.exists(file)) {
    file.remove(file)
    message("✓ Removed: ", file)
  }
}
message("✓ Temporary files removed.")

# ___________________
# Pipeline Complete
# ___________________

message("\n🎉 --- scVI/scArches Pipeline Finished Successfully! --- 🎉")