# -------------------------------------------------
# COVID-19 PBMC - Run Full scVI/scArches Pipeline
# -------------------------------------------------

# ___________________________________
# Step 1: Prepare Data (R -> Python)
# ___________________________________

message("\n--- Step 1: Preparing data for scVI (R -> Python) ---")
source("R/covid/scVI_Data_Preparation.R")

# _________________________________
# Step 2: Run Annotation (Python)
# _________________________________

message("\n--- Step 2: Running scVI/scArches in Python ---")

# Get Python from activated conda environment
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (conda_prefix == "") {
  stop("CONDA_PREFIX not set. Ensure conda environment is activated in SBATCH script.")
}

python_path <- file.path(conda_prefix, "bin", "python")
message("Using Python: ", python_path)

# Run Python script
python_script <- "python/covid/scVI_Annotation.py"
exit_code <- system(paste(python_path, python_script))

if (exit_code != 0) {
  stop("Python script failed with exit code ", exit_code)
}

message("--- Python script finished ---")

# _________________________________________
# Step 3: Integrate Results (Python -> R)
# _________________________________________

message("\n--- Step 3: Integrating results back into R (Python -> R) ---")
source("R/covid/scVI_Results_Integration.R")

# ____________________________________
# Step 4: Clean Up Intermediate Files
# ____________________________________

message("\n--- Step 4: Cleaning up all intermediate files ---")

# Add the two CSV files to the list of files to remove
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

# __________________
# Pipeline Complete
# __________________

message("\n🎉 --- scVI/scArches Pipeline Finished Successfully! --- 🎉")
