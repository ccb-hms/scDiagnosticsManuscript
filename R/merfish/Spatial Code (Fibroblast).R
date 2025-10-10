# -----------------------------------------------
# MERFISH Mouse Colon IBD - Figure 2 Spatial Plots 
# -----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(scDiagnostics)
library(scater) 

# -----------------------------------------------

# Load the processed SCE objects
# Assuming healthy_data is your reference and dss9_data is the single slice you want to plot.
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

# Make sure annotation columns are characters
healthy_data$tier2 <- as.character(healthy_data$tier2)
dss9_data$tier2 <- as.character(dss9_data$tier2)
dss9_data$azimuth_celltype_l1 <- as.character(dss9_data$azimuth_celltype_l1)

# __________________________________________________________________
# 1. Run Anomaly Detection and Add Status to Data
# __________________________________________________________________
# This is the essential step to generate the data for Panel B.

cat("Running anomaly detection on the 'azimuth_celltype_l1' column...\n")
anomaly_output <- detectAnomaly(
    query_data = dss9_data,
    reference_data = healthy_data, 
    query_cell_type_col = "azimuth_celltype_l1", 
    ref_cell_type_col = "tier2", 
    cell_types = c("Fibro 2", "Fibro 6"), 
    pc_subset = 1:10, 
    anomaly_threshold = 0.4
)

# Add the anomaly results to the dss9_data object
anomalous_barcodes_with_prefix <- c(
    names(anomaly_output$`Fibro 2`$query_anomaly[anomaly_output$`Fibro 2`$query_anomaly == TRUE]),
    names(anomaly_output$`Fibro 6`$query_anomaly[anomaly_output$`Fibro 6`$query_anomaly == TRUE])
)
anomalous_barcodes <- gsub("Query_", "", anomalous_barcodes_with_prefix)
dss9_data$anomaly_status <- "Typical"
dss9_data$anomaly_status[colnames(dss9_data) %in% anomalous_barcodes] <- "Anomalous"

cat("Anomaly status added to dss9_data.\n")


# __________________________________________________________________
# 2. Prepare a Single Data Frame for Plotting
# __________________________________________________________________
# We create one clean data frame containing everything ggplot needs.

plot_df <- data.frame(
    x = spatialCoords(dss9_data)[, "x"],
    y = spatialCoords(dss9_data)[, "y"],
    ground_truth_type = dss9_data$tier2,
    anomaly_status = dss9_data$anomaly_status
)

cat("Created plotting data frame with", nrow(plot_df), "cells.\n")


# __________________________________________________________________
# 3. Create the Two-Panel Figure
# __________________________________________________________________

# --- Define consistent aesthetic choices ---
point_size <- 0.5
plot_theme <- theme_void() + 
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom"
    )

# --- Panel A: The Ground Truth Map ---
# We will highlight the key pathological cell types (IAFs) and make others gray.
plot_df$ground_truth_highlight <- ifelse(
    grepl("IAF", plot_df$ground_truth_type), 
    "IAF (All)", # Group all IAFs together for a clear visual
    "Other"
)

panel_a_colors <- c("IAF (All)" = "#E31A1C", "Other" = "grey85")

# Plotting in two layers brings the highlighted points to the front.
panel_a <- ggplot() +
    geom_point(data = subset(plot_df, ground_truth_highlight == "Other"), 
               aes(x = x, y = y), color = "grey85", size = point_size) +
    geom_point(data = subset(plot_df, ground_truth_highlight == "IAF (All)"), 
               aes(x = x, y = y, color = ground_truth_highlight), size = point_size) +
    scale_color_manual(values = panel_a_colors, name = "Ground Truth") +
    labs(title = "A) Ground Truth Pathological Cells") +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size=4))) # Make legend points bigger


# --- Panel B: Your Package's Anomaly Map ---
panel_b_colors <- c("Typical" = "grey85", "Anomalous" = "#E31A1C") # Use the same red

panel_b <- ggplot() +
    geom_point(data = subset(plot_df, anomaly_status == "Typical"), 
               aes(x = x, y = y), color = "grey85", size = point_size) +
    geom_point(data = subset(plot_df, anomaly_status == "Anomalous"), 
               aes(x = x, y = y, color = anomaly_status), size = point_size) +
    scale_color_manual(values = panel_b_colors, name = "Package Output") +
    labs(title = "B) Anomalous Cells Detected") +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size=4)))

# __________________________________________________________________
# 4. Combine and Save the Final Figure
# __________________________________________________________________

final_spatial_figure <- panel_a + panel_b

# Display the combined figure
print(final_spatial_figure)

# --- Save the plot ---
ggsave("MERFISH_Figure2_Spatial_2panel.pdf", 
       plot = final_spatial_figure, 
       width = 10, 
       height = 5.5, 
       units = "in")

ggsave("MERFISH_Figure2_Spatial_2panel_300dpi.png", 
       plot = final_spatial_figure, 
       width = 10, 
       height = 5.5, 
       units = "in", 
       dpi = 300)

cat("\nSpatial figure saved.\n")