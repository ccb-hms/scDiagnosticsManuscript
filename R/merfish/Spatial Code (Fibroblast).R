# --------------------------------------------------
# MERFISH Mouse Colon IBD - Figure 2 Spatial Plots 
# --------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(scDiagnostics)
library(scater) 

# --------------------------------------------------

# Load the processed SCE objects
# Assuming healthy_data is your reference and dss9_data is the single slice you want to plot.
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

# Make sure annotation columns are characters
healthy_data$tier2 <- as.character(healthy_data$tier2)
dss9_data$tier2 <- as.character(dss9_data$tier2)
dss9_data$azimuth_celltype_l1 <- as.character(dss9_data$azimuth_celltype_l1)

# __________________________________________________________________
# 1. Run Anomaly Detection and Process Results
# __________________________________________________________________

# First, get the counts of each cell type in both datasets
ref_cell_counts <- table(healthy_data$tier2)
query_cell_counts <- table(dss9_data$azimuth_celltype_l1)

# Identify cell types with at least 20 cells in EACH dataset
abundant_ref_types <- names(ref_cell_counts[ref_cell_counts >= 10])
abundant_query_types <- names(query_cell_counts[query_cell_counts >= 10])

# The final list of cell types to run is the intersection of these two abundant sets
final_cell_types_to_run <- intersect(abundant_ref_types, abundant_query_types)

cat(
    "Found", length(final_cell_types_to_run), 
    "cell types with >= 20 cells in both reference and query datasets.\n"
)
cat("Running anomaly detection on these cell types...\n")

# Run anomaly detection using the filtered list of cell types
anomaly_output <- detectAnomaly(
    query_data = dss9_data,
    reference_data = healthy_data, 
    query_cell_type_col = "azimuth_celltype_l1", 
    ref_cell_type_col = "tier2", 
    cell_types = final_cell_types_to_run, # Use the robust, filtered list
    pc_subset = 1:3, 
    anomaly_threshold = 0.5, 
    max_cells_query = NULL,
    max_cells_ref = NULL
)

# Extract ALL barcodes that the package flagged as anomalous from the valid runs
all_anomalous_barcodes_with_prefix <- unlist(lapply(anomaly_output, function(res) {
    names(res$query_anomaly[res$query_anomaly == TRUE])
}))
all_anomalous_barcodes <- gsub("Query_", "", all_anomalous_barcodes_with_prefix)

cat(length(all_anomalous_barcodes), "total anomalous cells detected by the package.\n")

# __________________________________________________________________
# 2. Prepare a Single Data Frame for Plotting
# __________________________________________________________________
# Filter the detected anomalous cells to only include ground-truth IAFs

# Condition 1: Get barcodes of cells that are ground-truth IAFs
iaf_barcodes <- colnames(dss9_data)[grepl("IAF", dss9_data$tier2)]
cat(length(iaf_barcodes), "total ground-truth IAF cells in the data.\n")

# Condition 2: Find the intersection of all detected anomalous cells and the IAF cells
final_anomalous_barcodes <- intersect(all_anomalous_barcodes, iaf_barcodes)
cat(length(final_anomalous_barcodes), "cells are BOTH anomalous AND ground-truth IAFs. These will be plotted in red.\n")

# Add the refined anomaly status to the dss9_data object
dss9_data$anomaly_status <- "Typical"
dss9_data$anomaly_status[colnames(dss9_data) %in% final_anomalous_barcodes] <- "Anomalous"

# Create one clean data frame containing everything ggplot needs.
plot_df <- data.frame(
    x = spatialCoords(dss9_data)[, "x"],
    y = spatialCoords(dss9_data)[, "y"],
    ground_truth_type = dss9_data$tier2,
    anomaly_status = dss9_data$anomaly_status
)

cat("Created plotting data frame with", nrow(plot_df), "cells.\n")


# __________________________________________________________________
# 3. Create the Two-Panel Figure (No changes needed here)
# __________________________________________________________________

# --- Define consistent aesthetic choices ---
point_size <- 0.5
plot_theme <- theme_void() + 
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom"
    )

# --- Panel A: The Ground Truth Map ---
plot_df$ground_truth_highlight <- ifelse(
    grepl("IAF", plot_df$ground_truth_type), 
    "IAF (All)",
    "Other"
)

panel_a_colors <- c("IAF (All)" = "#E31A1C", "Other" = "grey85")

panel_a <- ggplot() +
    geom_point(data = subset(plot_df, ground_truth_highlight == "Other"), 
               aes(x = x, y = y), color = "grey85", size = point_size) +
    geom_point(data = subset(plot_df, ground_truth_highlight == "IAF (All)"), 
               aes(x = x, y = y, color = ground_truth_highlight), size = point_size) +
    scale_color_manual(values = panel_a_colors, name = "Ground Truth") +
    labs(title = "Ground Truth Pathological Cells") +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size=4)))

# --- Panel B: Your Package's Anomaly Map ---
panel_b_colors <- c("Typical" = "grey85", "Anomalous" = "#E31A1C")

panel_b <- ggplot() +
    geom_point(data = subset(plot_df, anomaly_status == "Typical"), 
               aes(x = x, y = y), color = "grey85", size = point_size) +
    geom_point(data = subset(plot_df, anomaly_status == "Anomalous"), 
               aes(x = x, y = y, color = anomaly_status), size = point_size) +
    scale_color_manual(values = panel_b_colors, name = "Package Output") +
    labs(title = "Anomalous IAFs Detected") +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size=4)))

# __________________________________________________________________
# 4. Combine and Save the Final Figure
# __________________________________________________________________

final_spatial_figure <- panel_a + panel_b

# Display the combined figure
print(final_spatial_figure)

# --- Save the plot ---
ggsave("MERFISH_Figure2_Spatial_2panel_IAF_filtered.pdf", 
       plot = final_spatial_figure, 
       width = 10, 
       height = 5.5, 
       units = "in")

ggsave("MERFISH_Figure2_Spatial_2panel_IAF_filtered_300dpi.png", 
       plot = final_spatial_figure, 
       width = 10, 
       height = 5.5, 
       units = "in", 
       dpi = 300)

cat("\nSpatial figure saved.\n")