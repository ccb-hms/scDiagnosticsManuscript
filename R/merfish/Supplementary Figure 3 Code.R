# -------------------------------------------------------
# MERFISH Mouse Colon IBD - Supplementary Figure 3 Code
# -------------------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(viridis)
library(scDiagnostics)

# -------------------------------------------------------

# _____________
# Data Loading
# _____________

cat("\nLoading datasets for inflamed fibroblast analysis...\n")
normal_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data   <- readRDS("data/merfish/dss9_data.rds")

cat("Datasets loaded.\n")

# ______
# Setup
# ______

inflamed_fb_mask <- dss9_data$tier2_merged == "Inflamed Fibroblast"
n_inflamed_fb <- sum(inflamed_fb_mask)
cat(sprintf("\n✓ Total inflamed fibroblasts: %d\n", n_inflamed_fb))

all_cell_types <- unique(normal_data$tier2_merged)
cat(sprintf("✓ Cell types in reference: %s\n", paste(all_cell_types, collapse = ", ")))

inflamed_fb_cellnames <- colnames(dss9_data)[inflamed_fb_mask]

# _______________________________
# Function to Analyze One Method
# _______________________________

analyze_method <- function(method_name, anomaly_result, pred_col_name) {
    
    cat(sprintf("\n=== %s ===\n", method_name))
    
    # Get predictions for inflamed fibroblasts
    names(dss9_data[[pred_col_name]]) <- colnames(dss9_data)
    pred_inflamed <- dss9_data[[pred_col_name]][inflamed_fb_mask]
    pred_unique <- unique(pred_inflamed)
    
    cat(sprintf("Predictions: %s\n", paste(pred_unique, collapse = ", ")))
    cat(sprintf("\n%% anomalous by predicted type (among inflamed fibroblasts):\n"))
    
    method_row <- c()
    
    for (pred_type in pred_unique) {
        # Which inflamed fibroblasts were predicted as this type?
        mask_inflamed <- pred_inflamed == pred_type
        inflamed_cells_of_type <- inflamed_fb_cellnames[mask_inflamed]
        
        # Get anomaly status from matching cell type detection
        if (!pred_type %in% names(anomaly_result)) {
            cat(sprintf("  WARNING: %s not in anomaly results\n", pred_type))
            next
        }
        
        anomaly_for_type <- anomaly_result[[pred_type]][["query_anomaly"]]
        
        # Remove Query_ prefix
        names(anomaly_for_type) <- gsub("Query_", "", names(anomaly_for_type))
        
        # Extract anomaly status
        inflamed_cells_clean <- gsub("Query_", "", inflamed_cells_of_type)
        inflamed_anomaly <- anomaly_for_type[inflamed_cells_clean]
        
        n_cells <- sum(mask_inflamed)
        n_anom <- sum(inflamed_anomaly, na.rm = TRUE)
        pct_anom <- (n_anom / n_cells) * 100
        
        method_row[pred_type] <- pct_anom
        cat(sprintf("  %s: %d/%d (%.1f%%)\n", pred_type, n_anom, n_cells, pct_anom))
    }
    
    return(method_row)
}

# ______________________________________
# Run Anomaly Detection for All Methods
# ______________________________________

cat("\n" %+% paste(rep("=", 60), collapse = "") %+% "\n")
cat("ANOMALY DETECTION\n")
cat(paste(rep("=", 60), collapse = "") %+% "\n")

cat("\nDetecting anomalies for Azimuth...\n")
anomaly_azimuth <- detectAnomaly(
    reference_data = normal_data,
    query_data = dss9_data,
    ref_cell_type_col = "tier2_merged",
    query_cell_type_col = "azimuth_celltype_l1_merged",
    cell_types = all_cell_types,
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_query = NULL,
    max_cells_ref = NULL
)

cat("\nDetecting anomalies for SingleR...\n")
anomaly_singler <- detectAnomaly(
    reference_data = normal_data,
    query_data = dss9_data,
    ref_cell_type_col = "tier2_merged",
    query_cell_type_col = "singler_annotations_merged",
    cell_types = all_cell_types,
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_query = NULL,
    max_cells_ref = NULL
)

cat("\nDetecting anomalies for CellTypist...\n")
anomaly_celltypist <- detectAnomaly(
    reference_data = normal_data,
    query_data = dss9_data,
    ref_cell_type_col = "tier2_merged",
    query_cell_type_col = "celltypist_predicted_labels_merged",
    cell_types = all_cell_types,
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_query = NULL,
    max_cells_ref = NULL
)

cat("\nDetecting anomalies for scArches...\n")
anomaly_scarches <- detectAnomaly(
    reference_data = normal_data,
    query_data = dss9_data,
    ref_cell_type_col = "tier2_merged",
    query_cell_type_col = "scvi_prediction_merged",
    cell_types = all_cell_types[1:9],
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_query = NULL,
    max_cells_ref = NULL
)

cat("\n✓ All anomaly detections complete\n")

# _____________________
# Analyze Each Method
# _____________________

cat("\n" %+% paste(rep("=", 60), collapse = "") %+% "\n")
cat("ANALYSIS RESULTS\n")
cat(paste(rep("=", 60), collapse = "") %+% "\n")

azimuth_row <- analyze_method("Azimuth", anomaly_azimuth, "azimuth_celltype_l1_merged")
singler_row <- analyze_method("SingleR", anomaly_singler, "singler_annotations_merged")
celltypist_row <- analyze_method("CellTypist", anomaly_celltypist, "celltypist_predicted_labels_merged")
scarches_row <- analyze_method("scArches", anomaly_scarches, "scvi_prediction_merged")

# _____________________________
# Create Combined Heatmap Data
# _____________________________

cat("\n" %+% paste(rep("=", 60), collapse = "") %+% "\n")
cat("CREATING COMBINED HEATMAP\n")
cat(paste(rep("=", 60), collapse = "") %+% "\n\n")

# Get all unique predicted types across all methods
all_pred_types <- unique(c(names(azimuth_row), names(singler_row), 
                            names(celltypist_row), names(scarches_row)))
all_pred_types <- sort(all_pred_types)

# ===== Matrix 1: Total counts (for coloring by sample size) =====
heatmap_counts <- matrix(NA, nrow = 4, ncol = length(all_pred_types))
rownames(heatmap_counts) <- c("Azimuth", "SingleR", "CellTypist", "scArches")
colnames(heatmap_counts) <- all_pred_types

# ===== Matrix 2: Display text (x / y format) =====
heatmap_display <- matrix(NA, nrow = 4, ncol = length(all_pred_types))
rownames(heatmap_display) <- c("Azimuth", "SingleR", "CellTypist", "scArches")
colnames(heatmap_display) <- all_pred_types

# ===== Azimuth =====
cat("Processing Azimuth...\n")
names(dss9_data$azimuth_celltype_l1_merged) <- colnames(dss9_data)
pred_inflamed_azimuth <- dss9_data$azimuth_celltype_l1_merged[inflamed_fb_mask]
for (pred_type in all_pred_types) {
    mask_inflamed <- pred_inflamed_azimuth == pred_type
    if (sum(mask_inflamed) > 0) {
        inflamed_cells_of_type <- inflamed_fb_cellnames[mask_inflamed]
        anomaly_for_type <- anomaly_azimuth[[pred_type]][["query_anomaly"]]
        names(anomaly_for_type) <- gsub("Query_", "", names(anomaly_for_type))
        inflamed_cells_clean <- gsub("Query_", "", inflamed_cells_of_type)
        inflamed_anomaly <- anomaly_for_type[inflamed_cells_clean]
        n_anom <- sum(inflamed_anomaly, na.rm = TRUE)
        n_total <- sum(mask_inflamed)
        heatmap_display["Azimuth", pred_type] <- sprintf("%d/%d", n_anom, n_total)
        heatmap_counts["Azimuth", pred_type] <- n_total
    }
}

# ===== SingleR =====
cat("Processing SingleR...\n")
names(dss9_data$singler_annotations_merged) <- colnames(dss9_data)
pred_inflamed_singler <- dss9_data$singler_annotations_merged[inflamed_fb_mask]
for (pred_type in all_pred_types) {
    mask_inflamed <- pred_inflamed_singler == pred_type
    if (sum(mask_inflamed) > 0) {
        inflamed_cells_of_type <- inflamed_fb_cellnames[mask_inflamed]
        anomaly_for_type <- anomaly_singler[[pred_type]][["query_anomaly"]]
        names(anomaly_for_type) <- gsub("Query_", "", names(anomaly_for_type))
        inflamed_cells_clean <- gsub("Query_", "", inflamed_cells_of_type)
        inflamed_anomaly <- anomaly_for_type[inflamed_cells_clean]
        n_anom <- sum(inflamed_anomaly, na.rm = TRUE)
        n_total <- sum(mask_inflamed)
        heatmap_display["SingleR", pred_type] <- sprintf("%d/%d", n_anom, n_total)
        heatmap_counts["SingleR", pred_type] <- n_total
    }
}

# ===== CellTypist =====
cat("Processing CellTypist...\n")
names(dss9_data$celltypist_predicted_labels_merged) <- colnames(dss9_data)
pred_inflamed_celltypist <- dss9_data$celltypist_predicted_labels_merged[inflamed_fb_mask]
for (pred_type in all_pred_types) {
    mask_inflamed <- pred_inflamed_celltypist == pred_type
    if (sum(mask_inflamed) > 0) {
        inflamed_cells_of_type <- inflamed_fb_cellnames[mask_inflamed]
        anomaly_for_type <- anomaly_celltypist[[pred_type]][["query_anomaly"]]
        names(anomaly_for_type) <- gsub("Query_", "", names(anomaly_for_type))
        inflamed_cells_clean <- gsub("Query_", "", inflamed_cells_of_type)
        inflamed_anomaly <- anomaly_for_type[inflamed_cells_clean]
        n_anom <- sum(inflamed_anomaly, na.rm = TRUE)
        n_total <- sum(mask_inflamed)
        heatmap_display["CellTypist", pred_type] <- sprintf("%d/%d", n_anom, n_total)
        heatmap_counts["CellTypist", pred_type] <- n_total
    }
}

# ===== scArches =====
cat("Processing scArches...\n")
names(dss9_data$scvi_prediction_merged) <- colnames(dss9_data)
pred_inflamed_scarches <- dss9_data$scvi_prediction_merged[inflamed_fb_mask]
for (pred_type in all_pred_types) {
    mask_inflamed <- pred_inflamed_scarches == pred_type
    if (sum(mask_inflamed) > 0) {
        inflamed_cells_of_type <- inflamed_fb_cellnames[mask_inflamed]
        anomaly_for_type <- anomaly_scarches[[pred_type]][["query_anomaly"]]
        names(anomaly_for_type) <- gsub("Query_", "", names(anomaly_for_type))
        inflamed_cells_clean <- gsub("Query_", "", inflamed_cells_of_type)
        inflamed_anomaly <- anomaly_for_type[inflamed_cells_clean]
        n_anom <- sum(inflamed_anomaly, na.rm = TRUE)
        n_total <- sum(mask_inflamed)
        heatmap_display["scArches", pred_type] <- sprintf("%d/%d", n_anom, n_total)
        heatmap_counts["scArches", pred_type] <- n_total
    }
}

cat("\nCombined heatmap data (colored by sample size, text shows anomalous/total):\n")
print(heatmap_display)

# ===== Fill NA values with 0/0 =====
cat("Filling NA cells with 0/0...\n")

# Replace NA in display matrix with "0/0"
heatmap_display[is.na(heatmap_display)] <- "0/0"

# Replace NA in count matrix with 0 (will get colored as lowest value)
heatmap_counts[is.na(heatmap_counts)] <- 0

cat("\nCombined heatmap data (colored by sample size, text shows anomalous/total):\n")
print(heatmap_display)

# ________________
# Create Heatmap
# ________________

cat("\nCreating combined heatmap...\n")

dir.create("figures/supp/merfish", showWarnings = FALSE, recursive = TRUE)

png("figures/supp/merfish/Fig_S3_annotation_heatmap.png", 
    width = 2100, height = 700, res = 150)

pheatmap(heatmap_counts,
         main = "Inflamed Fibroblasts: Assigned Cell Type & Anomaly Detection\n(Color Intensity = # Inflamed Fibroblasts, Text = Anomalous / Total)",
         display_numbers = heatmap_display,
         breaks = seq(0, max(heatmap_counts, na.rm = TRUE), length.out = 100),
         color = viridis(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize = 12,
         fontsize_row = 13,
         fontsize_col = 11,
         fontsize_number = 10,
         angle_col = 45,
         cellwidth = 60,
         cellheight = 50,
         margins = c(12, 8))

dev.off()

cat("\n✓ Combined heatmap saved to figures/supp/merfish/Fig_S3_annotation_heatmap.png\n")