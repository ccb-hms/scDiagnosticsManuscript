# -----------------------------------------------
# MERFISH Mouse Colon IBD - IAE Spatial Figure
# -----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(scater) 
library(scDiagnostics)

# -----------------------------------------------

# Load the processed SCE objects
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

# Make sure annotation columns are characters
healthy_data$tier2 <- as.character(healthy_data$tier2)
dss9_data$tier2 <- as.character(dss9_data$tier2)
dss9_data$azimuth_celltype_l2 <- as.character(dss9_data$azimuth_celltype_l2) 

# __________________________________________________________________
# 1. Run Anomaly Detection for EPITHELIAL Lineage
# __________________________________________________________________
# This is the key change. We now target the mislabeled epithelial cells.
# IAEs were mislabeled as Colonocytes, TA, and Stem cells.

cat("--- Running anomaly detection on naive epithelial predictions ---\n")
anomaly_output_iae <- detectAnomaly(
    query_data = dss9_data,
    reference_data = healthy_data, 
    query_cell_type_col = "azimuth_celltype_l2", 
    ref_cell_type_col = "tier2", 
    cell_types = c("Colonocytes", "TA", "Stem cells"), 
    pc_subset = 1:5, 
    anomaly_threshold = 0.5
)

# Add the anomaly results to the dss9_data object in a NEW column
anomalous_barcodes_iae <- gsub("Query_", "", c(
    names(anomaly_output_iae$Colonocytes$query_anomaly[anomaly_output_iae$Colonocytes$query_anomaly == TRUE]),
    names(anomaly_output_iae$TA$query_anomaly[anomaly_output_iae$TA$query_anomaly == TRUE]),
    names(anomaly_output_iae$`Stem cells`$query_anomaly[anomaly_output_iae$`Stem cells`$query_anomaly == TRUE])
))
dss9_data$anomaly_status_iae <- "Typical"
dss9_data$anomaly_status_iae[colnames(dss9_data) %in% anomalous_barcodes_iae] <- "Anomalous"

cat("Epithelial anomaly status added to dss9_data.\n")


# __________________________________________________________________
# 2. Prepare a Single Data Frame for Plotting
# __________________________________________________________________
plot_df_iae <- data.frame(
    x = spatialCoords(dss9_data)[, "x"],
    y = spatialCoords(dss9_data)[, "y"],
    ground_truth_type = dss9_data$tier2,
    anomaly_status = dss9_data$anomaly_status_iae
)

cat("Created plotting data frame with", nrow(plot_df_iae), "cells.\n")


# __________________________________________________________________
# 3. Create the Two-Panel Figure for IAEs
# __________________________________________________________________

# --- Define consistent aesthetic choices ---
point_size <- 0.5
plot_theme <- theme_void() + 
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom"
    )

# --- Panel A: The Ground Truth Map (Highlighting IAEs) ---
# We will highlight all IAE types (1, 2, and 3).
plot_df_iae$ground_truth_highlight <- ifelse(
    grepl("IAE", plot_df_iae$ground_truth_type), 
    "IAE (All)", 
    "Other"
)

panel_a_colors <- c("IAE (All)" = "#E31A1C", "Other" = "grey85")

panel_a_iae <- ggplot() +
    geom_point(data = subset(plot_df_iae, ground_truth_highlight == "Other"), 
               aes(x = x, y = y), color = "grey85", size = point_size) +
    geom_point(data = subset(plot_df_iae, ground_truth_highlight == "IAE (All)"), 
               aes(x = x, y = y, color = ground_truth_highlight), size = point_size) +
    scale_color_manual(values = panel_a_colors, name = "Ground Truth") +
    labs(title = "A) Ground Truth Pathological Epithelial Cells (IAEs)") +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size=4)))


# --- Panel B: Your Package's Anomaly Map ---
# This panel shows the results from our IAE-focused anomaly detection run.
panel_b_colors <- c("Typical" = "grey85", "Anomalous" = "#E31A1C")

panel_b_iae <- ggplot() +
    geom_point(data = subset(plot_df_iae, anomaly_status == "Typical"), 
               aes(x = x, y = y), color = "grey85", size = point_size) +
    geom_point(data = subset(plot_df_iae, anomaly_status == "Anomalous"), 
               aes(x = x, y = y, color = anomaly_status), size = point_size) +
    scale_color_manual(values = panel_b_colors, name = "Package Output") +
    labs(title = "B) Anomalous Cells Detected in Epithelium") +
    plot_theme +
    guides(color = guide_legend(override.aes = list(size=4)))

# __________________________________________________________________
# 4. Combine and Save the Final Figure
# __________________________________________________________________

final_iae_figure <- panel_a_iae + panel_b_iae

# Display the combined figure
print(final_iae_figure)

# --- Save the plot ---
ggsave("MERFISH_IAE_Spatial_Figure.pdf", 
       plot = final_iae_figure, 
       width = 10, 
       height = 5.5, 
       units = "in")

ggsave("MERFISH_IAE_Spatial_Figure_300dpi.png", 
       plot = final_iae_figure, 
       width = 10, 
       height = 5.5, 
       units = "in", 
       dpi = 300)

cat("\nIAE spatial figure saved.\n")