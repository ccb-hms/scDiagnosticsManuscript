# -------------------------------------------------
# COVID-19 PBMC - Run Full scVI/scArches Pipeline
# -------------------------------------------------

# __________________________________________________
# Step 1: Prepare Data (R -> Python)
# __________________________________________________

message("\n--- Step 1: Preparing data for scVI (R -> Python) ---")
source("R/covid/scVI_Data_Preparation.R")

# __________________________________________________
# Step 2: Run Annotation (Python)
# __________________________________________________

message("\n--- Step 2: Running scVI/scArches in Python ---")
# This command assumes a conda environment named 'scvi-env' has been set up
# by the sbatch script.
system("python python/covid/scVI_Annotation.py")
message("--- Python script finished ---")

# __________________________________________________
# Step 3: Integrate Results (Python -> R)
# __________________________________________________

message("\n--- Step 3: Integrating results back into R (Python -> R) ---")
source("R/covid/scVI_Results_Integration.R")

# __________________________________________________
# Step 4: Clean Up Intermediate Files
# __________________________________________________

message("\n--- Step 4: Cleaning up all intermediate files ---")

# --- UPDATED: Add the two CSV files to the list of files to remove ---
files_to_remove <- c(
  "data/covid/scvi_reference_data.h5ad",
  "data/covid/scvi_query_data.h5ad",
  "data/covid/scvi_predictions.csv",    
  "data/covid/scvi_reference_umap.csv" 
)

for (file in files_to_remove) {
  if (file.exists(file)) {
    file.remove(file)
    message("✓ Removed: ", file)
  }
}
message("✓ Temporary files removed.")

# _________________________
# Pipeline Complete
# _________________________

message("\n🎉 --- scVI/scArches Pipeline Finished Successfully! --- 🎉")
