# --------------------------------------
# COVID-19 PBMC - Preliminary Analyses
# --------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scDiagnostics)
library(ggplot2)
library(scran)
library(scater)
library(ComplexHeatmap)
library(circlize)

# --------------------------------------

# __________
# Load Data
# __________

# Load the processed SCE objects
normal_data <- readRDS("data/covid/normal_data_sce.rds")
covid_data <- readRDS("data/covid/covid_data_sce.rds")

# Setting the seed
set.seed(0)

# __________________
# Anomaly Detection
# __________________

# Anomaly Detection
anomaly_output <- detectAnomaly(query_data = covid_data,
                                reference_data = normal_data, 
                                query_cell_type_col = "azimuth_celltype_l1", 
                                ref_cell_type_col = "author_cell_type", 
                                cell_types = c("CD14_mono"), 
                                pc_subset = 1:3,
                                n_tree = 500,
                                anomaly_threshold = 0.5, 
                                max_cells_query = NULL, 
                                max_cells_ref = NULL)

anomaly_plot <- plot(anomaly_output,
                     cell_type = "CD14_mono",
                     pc_subset = 1:3,
                     data_type = "query",
                     n_tree = 500,
                     diagonal_facet = "ridge",
                     upper_facet = "blank", 
                     max_cells_query = 1000)

# Print the final plot
anomaly_plot

# Save the final plot
ggsave("figures/covid/anomaly_plot.png", width = 12, height = 8, dpi = 600)

