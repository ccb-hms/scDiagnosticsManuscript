# ----------------------------------------------------------
# MERFISH Mouse Colon IBD - CellTypist Results Integration
# ----------------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scater)

# Source auxiliary functions
source("R/auxiliary/addMergedCellTypes.R")

# ----------------------------------------------------------

# ________________________
# Setup and Data Loading
# ________________________

message("--- Initializing CellTypist Integration (MERFISH) ---")

# Define file paths
ORIGINAL_QUERY_FILE <- "data/merfish/dss9_data.rds"
ANNOTATED_H5AD      <- "data/merfish/annotated_query_data.h5ad"
TEMP_CSV            <- "data/merfish/temp_annotations.csv"
TEMP_PY_SCRIPT      <- "python/merfish/export_to_csv.py"

# __________________________________________________
# Step A: Export Annotations to CSV using Python
# __________________________________________________

message("\n--- Step A: Converting H5AD to CSV (External Python) ---")

# 1. Create the helper Python script dynamically
py_code <- c(
  "import scanpy as sc",
  "import pandas as pd",
  "import sys",
  "import os",
  "",
  paste0("h5ad_path = '", ANNOTATED_H5AD, "'"),
  paste0("csv_path = '", TEMP_CSV, "'"),
  "",
  "print(f'Loading {h5ad_path}...')",
  "if not os.path.exists(h5ad_path):",
  "    print('Error: H5AD file not found')",
  "    sys.exit(1)",
  "",
  "try:",
  "    adata = sc.read_h5ad(h5ad_path)",
  "    # Extract only the CellTypist annotation columns",
  "    cols_to_keep = [c for c in adata.obs.columns if 'predicted_labels' in c or 'majority_voting' in c or 'conf_score' in c]",
  "    print(f'Found columns: {cols_to_keep}')",
  "    ",
  "    # Save to CSV",
  "    adata.obs[cols_to_keep].to_csv(csv_path)",
  "    print(f'Saved annotations to {csv_path}')",
  "except Exception as e:",
  "    print(f'Error: {e}')",
  "    sys.exit(1)"
)

# Write python script
if(!dir.exists("python/merfish")) dir.create("python/merfish", recursive = TRUE)
writeLines(py_code, TEMP_PY_SCRIPT)

# 2. Run the script using conda run
# We use the same environment that worked for the annotation step
cmd <- paste("conda run -n celltypist_env python", TEMP_PY_SCRIPT)
message("Running command: ", cmd)

status <- system(cmd, intern = FALSE)

if (status != 0) {
  stop("Python CSV conversion failed! The environment 'celltypist_env' might not be reachable or the H5AD file is missing.")
}

if (!file.exists(TEMP_CSV)) {
  stop("Python ran but the CSV file was not found.")
}
message("✓ CSV conversion successful.")

# __________________________________________
# Step B: Load CSV and Integrate into R
# __________________________________________

message("\n--- Step B: Loading CSV into R ---")

# Load Original Data
message("Loading original DSS9 SCE object...")
dss9_data_original <- readRDS(ORIGINAL_QUERY_FILE)

# Load New Annotations
message("Loading annotation CSV...")
annotations <- read.csv(TEMP_CSV, row.names = 1)
message("Loaded annotations for ", nrow(annotations), " cells.")

# _______________________________
# Data Validation and Alignment
# _______________________________

message("\n--- Validating Cell Alignment ---")

cell_names_original <- colnames(dss9_data_original)
cell_names_annotated <- rownames(annotations)

# Intersect cells
common_cells <- intersect(cell_names_original, cell_names_annotated)
message("Common cells found: ", length(common_cells))

if(length(common_cells) == 0) {
  stop("No common cells found! The barcodes in the H5AD do not match the RDS.")
}

# Align both objects
dss9_data_original <- dss9_data_original[, common_cells]
annotations <- annotations[common_cells, , drop = FALSE]

message("Alignment complete.")

# ____________________________________
# Annotation Renaming and Merging
# ____________________________________

message("\n--- Formatting Annotations ---")

# 1. Rename columns (Add 'celltypist_' prefix)
name_map <- c(
  "predicted_labels"  = "celltypist_predicted_labels",
  "majority_voting"   = "celltypist_majority_voting",
  "conf_score"        = "celltypist_conf_score"
)

existing_cols <- colnames(annotations)
new_names <- existing_cols
for(col in existing_cols) {
  if(col %in% names(name_map)) {
    new_names[which(existing_cols == col)] <- name_map[col]
  }
}
colnames(annotations) <- new_names
message("Renamed columns: ", paste(new_names, collapse=", "))

# 2. Add to SCE
colData(dss9_data_original) <- cbind(colData(dss9_data_original), annotations)
message("✓ Annotations added to object.")

# 3. Create Merged Categories
message("Creating merged cell types...")

if ("celltypist_predicted_labels" %in% colnames(colData(dss9_data_original))) {
  dss9_data_original <- addMergedCellTypes(
    sce_object = dss9_data_original,
    input_col_name = "celltypist_predicted_labels", 
    dataset = "MERFISH"
  )
}

# __________________________
# Final Summary and Saving
# __________________________

message("\n--- Final Summary ---")

if ("celltypist_predicted_labels_merged" %in% colnames(colData(dss9_data_original))) {
  cat("Top 10 predicted cell types:\n")
  print(head(sort(table(dss9_data_original$celltypist_predicted_labels_merged), decreasing = TRUE), 10))
}

# SAVE THE FILE
message("\nSaving updated object...")
saveRDS(dss9_data_original, ORIGINAL_QUERY_FILE)
message("✓ Saved updated object to: ", ORIGINAL_QUERY_FILE)

# Cleanup
if(file.exists(TEMP_CSV)) file.remove(TEMP_CSV)
if(file.exists(TEMP_PY_SCRIPT)) file.remove(TEMP_PY_SCRIPT)

message("\n✓✓✓ Integration Complete! ✓✓✓")