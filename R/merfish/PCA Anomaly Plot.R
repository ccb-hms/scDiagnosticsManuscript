# ----------------------------------------------------------
# MERFISH Mouse Colon IBD - Anomaly Detection & PCA Plot
# ----------------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scDiagnostics)
library(ggplot2)

# ----------------------------------------------------------

# ________________________
# Data Loading & Lumping
# ________________________

cat("\n--- Loading MERFISH Datasets ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

# Setting the seed
set.seed(0)

cat("\n--- Grouping Fibroblast Lineage ---\n")
# Logic: Grab anything that looks like a Fibroblast, IAF, FRC, or Pericyte
# This ensures we capture the transition states.
is_fibro_lineage <- function(x) {
    grepl("^Fibro|^IAF", x)
}

# Apply to Healthy Data
healthy_data$analysis_class <- ifelse(is_fibro_lineage(healthy_data$tier2), 
                                      "Fibroblast_Lineage", "Other")

# Apply to DSS9 Data (Includes IAFs)
dss9_data$analysis_class <- ifelse(is_fibro_lineage(dss9_data$tier2), 
                                   "Fibroblast_Lineage", "Other")

cat(sprintf("Healthy Fibroblasts in Reference: %d\n", sum(healthy_data$analysis_class == "Fibroblast_Lineage")))
cat(sprintf("Total Fibroblasts in Query:       %d\n", sum(dss9_data$analysis_class == "Fibroblast_Lineage")))


# ________________________________________
# Anomaly Detection on Fibroblast Lineage
# ________________________________________

cat("\n--- Running Anomaly Detection on Fibroblast Lineage ---\n")

# We compare the lumped "Fibroblast_Lineage" from Query (Sick) 
# against the "Fibroblast_Lineage" from Reference (Healthy).
anomaly_output <- detectAnomaly(
    query_data = dss9_data,
    reference_data = healthy_data, 
    query_cell_type_col = "analysis_class", 
    ref_cell_type_col = "analysis_class", 
    cell_types = "Fibroblast_Lineage", # The lumped category
    pc_subset = 1:4,
    n_tree = 500,
    anomaly_threshold = 0.5, 
    max_cells_query = NULL, 
    max_cells_ref = NULL
)


# ______________________
# Generate Anomaly Plot
# ______________________

cat("\n--- Generating PCA Anomaly Plot ---\n")

# Use the plot() method from scDiagnostics.
# This visualizes the Anomaly Score (Red = High Anomaly, Blue/Grey = Low).
anomaly_plot <- plot(
    anomaly_output,
    cell_type = "Fibroblast_Lineage", # Plot the lumped category
    pc_subset = 1:4,
    data_type = "query",  # Show the query cells colored by anomaly score
    n_tree = 500,
    diagonal_facet = "density",
    upper_facet = "blank",
    max_cells_query = 2000 # Downsample slightly for clearer dots
)

# Add styling and titles
anomaly_plot <- anomaly_plot + 
    labs(
        title = "Anomaly Detection: Fibroblast Lineage",
        subtitle = "Red indicates deviation from the Healthy Fibroblast Reference"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, face = "italic")
    )


# ______________________
# Display and Save Plot
# ______________________

cat("\n--- Displaying and Saving Plot ---\n")

print(anomaly_plot)

dir.create("figures/merfish", showWarnings = FALSE, recursive = TRUE)

ggsave("figures/merfish/merfish_pca_anomaly_plot_lumped.png", 
       plot = anomaly_plot, 
       width = 12, 
       height = 8, 
       dpi = 600)

cat("✓ Anomaly PCA plot saved to figures/merfish/ directory.\n")