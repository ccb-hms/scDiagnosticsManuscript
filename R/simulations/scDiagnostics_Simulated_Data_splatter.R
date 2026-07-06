# ------------------------------------------------
# scDiagnostics - Simulated Data (splatter)
# ------------------------------------------------

# Load libraries
library(splatter)
library(scuttle)
library(scater)
library(SingleR)
library(scDiagnostics)
library(ggplot2)
library(ggridges)
library(dplyr)
library(patchwork)

# Source files
source("R/simulations/plotCellTypePCA_custom.R")
source("R/simulations/plotIsoForest_custom.R")

# ------------------------------------------------

# ________________
# Data Simulation
# ________________

# Setting the seed
set.seed(100)

# Load a mock reference dataset
sce_ref <- mockSCE()

# Estimate parameters from the dataset
params <- splatEstimate(sce_ref)

# Customize SplatParams for batches, groups, differential expression, and outliers
params <- setParams(params,
                    batchCells = c(500, 500),              # Two batches, each with 500 cells
                    batch.facLoc = 0,                      # Moderate batch effect location
                    batch.facScale = 0,                    # Low variability in batch effects
                    group.prob = c(1/3, 1/3, 1/3),         # Cell types
                    de.prob = c(0.1, 0.2, 0.2),            # Group-specific DE probabilities
                    de.facLoc = c(0.250, 0.375, 0.375),    # Different DE factor locations
                    de.facScale = c(0.2, 0.3, 0.4),        # Different DE factor scales
                    out.prob = 0,                          # 0% probability of a gene being an outlier
                    out.facLoc = 4,                        # Outlier factor location
                    out.facScale = 0.5                     # Outlier factor scale
)

# Simulate reference and query data 
simulated_data <- splatSimulate(params, 
                                method = "groups",
                                verbose = FALSE)

# ________________________________
# Identical Cell Types - PCA Plot
# ________________________________

# Normalize counts and run PCA for visualization
reference_data <- simulated_data[, simulated_data$Batch == "Batch1"]
reference_data <- logNormCounts(reference_data)
reference_data <- runPCA(reference_data, ncomponents = 10)

query_data <- simulated_data[, simulated_data$Batch == "Batch2"]
query_data <- logNormCounts(query_data)
query_data <- runPCA(query_data, ncomponents = 10)

# Modify SCE objects column data
colData(reference_data)$Cell_Type <- colData(reference_data)$Group
colData(reference_data)$Group <- NULL
levels(colData(reference_data)$Cell_Type) <- c("Cell Type C", "Cell Type A", "Cell Type B")
colData(reference_data)$Cell_Type <- as.character(colData(reference_data)$Cell_Type)

colData(query_data)$Cell_Type <- colData(query_data)$Group
colData(query_data)$Group <- NULL
levels(colData(query_data)$Cell_Type) <- c("Cell Type C", "Cell Type A", "Cell Type B")
colData(query_data)$Cell_Type <- as.character(colData(query_data)$Cell_Type)

# Add SingleR annotations
SingleR_annotation <- SingleR(query_data, reference_data, 
                              labels = reference_data$Cell_Type)
query_data$SingleR_annotation <- SingleR_annotation$labels

# Generate color of cell types
cell_type_colors <- scDiagnostics:::generateColors(1:20, 
                                                   paired = TRUE)
cell_type_colors <- cell_type_colors[c(1:2, 5:6, 9:10)]

# PCA projection diagnostics (same cell types)
plotCellTypePCA_custom(
    query_data = query_data,
    reference_data = reference_data, 
    query_cell_type_col = "SingleR_annotation",
    ref_cell_type_col = "Cell_Type", 
    pc_subset = 1:2,
    assay_name = "logcounts", 
    diagonal_facet = "ridge", 
    upper_facet = "blank", 
    cell_type_colors = cell_type_colors)
ggsave("figures/simulation/cell_type_pca_original.png", 
    dpi = 600, 
    width = 6, height = 4)

# _____________________________
# Missing Cell Type - PCA Plot
# _____________________________

# Normalize counts and run PCA for visualization
reference_data_missing <- reference_data[, reference_data$Cell_Type != "Cell Type C"]
reference_data_missing <- runPCA(reference_data_missing, ncomponents = 10)

# Add SingleR annotations
SingleR_annotation <- SingleR(query_data, reference_data_missing, 
                              labels = reference_data_missing$Cell_Type)
query_data_missing <- query_data
query_data_missing$SingleR_annotation <- SingleR_annotation$labels
query_data_missing$SingleR_annotation[
    query_data_missing$SingleR_annotation == "Cell Type A" &
        query_data$SingleR_annotation == "Cell Type C"] <- 
    "Cell Type C \n Misannotated as Cell Type A"
query_data_missing$SingleR_annotation[
    query_data_missing$SingleR_annotation == "Cell Type B" &
        query_data$SingleR_annotation == "Cell Type C"] <- 
    "Cell Type C \n Misannotated as Cell Type B"

# Generate color of cell types
cell_type_colors <- scDiagnostics:::generateColors(1:20, paired = TRUE)
cell_type_colors <- c(cell_type_colors[1:2], "#3D7FFF", 
                      cell_type_colors[5:6], "#4DE280")

# PCA projection diagnostics (same cell types)
plotCellTypePCA_custom(
    query_data = query_data_missing,
    reference_data = reference_data_missing, 
    query_cell_type_col = "SingleR_annotation",
    ref_cell_type_col = "Cell_Type", 
    pc_subset = 1:2,
    assay_name = "logcounts", 
    diagonal_facet = "ridge", 
    upper_facet = "blank", 
    cell_type_colors = cell_type_colors)
ggsave("figures/simulation/cell_type_pca.png", 
    dpi = 600, 
    width = 6, height = 4)

# _____________________________________________
# Identical Cell Types - Isolation Forest Plot
# _____________________________________________

# Isolation Forest (same cell types)
anomaly_output <- detectAnomaly(reference_data = reference_data,
                                query_data = query_data,
                                ref_cell_type_col = "Cell_Type",
                                query_cell_type_col = "SingleR_annotation",
                                pc_subset = 1:2,
                                n_tree = 1000,
                                anomaly_threshold = 0.625)

# Plot the output for cell type A
gradient_colors <- c(
    "#7B9FC9", "#7B9FC9", "#7B9FC9",   
    "#9EB7D9", "#BAD0EB", "#CEE4F5",   
    "#E2F1FB", "white", "white")
stops <- c(
    0.00, 0.1, 0.2,
    0.3, 0.375, 0.425, 0.525,
    0.625, 0.675, 1.00
)
plotIsoForest_custom(
    anomaly_output,
    cell_type = "Cell Type A",
    pc_subset = 1:2,
    data_type = "query",
    n_tree = 1000,
    diagonal_facet = "blank",
    upper_facet = "blank",
    tile_colors = gradient_colors, 
    tile_stops = stops)
ggsave("figures/simulation/cell_type_A_isoForest_original.png", 
    dpi = 600, 
    width = 6, height = 4)

# Plot the output for cell type B
gradient_colors <- colorRampPalette(c(
    "#A1D99B", "#A1D99B", "#A1D99B",  
    "#A1D99B", "#B4E0AC", "#BDE5B3", 
    "#C7E9C0",  "white", "white"))(9)
stops <- c(
    0.00, 0.1, 0.2,
    0.3, 0.375, 0.425, 0.525,
    0.6, 0.625, 1.00
)
plotIsoForest_custom(
    anomaly_output,
    cell_type = "Cell Type B",
    pc_subset = 1:2,
    data_type = "query",
    n_tree = 1000,
    diagonal_facet = "blank",
    upper_facet = "blank",
    tile_colors = gradient_colors, 
    tile_stops = stops)
ggsave("figures/simulation/cell_type_B_isoForest_original.png", 
    dpi = 600, 
    width = 6, height = 4)


# Add density plot for Cell Type A anomaly scores
anomaly_scores_A <- anomaly_output$`Cell Type A`$query_anomaly_scores
df_scores_A <- data.frame(
    anomaly_score = anomaly_scores_A,
    cell_type = "Cell Type A"
)

# Calculate density
density_calc <- density(df_scores_A$anomaly_score)
density_df <- data.frame(
    x = density_calc$x,
    y = density_calc$y,
    fill_color = ifelse(density_calc$x < 0.6, "#7B9FC9", "#FF6B6B"),
    line_color = ifelse(density_calc$x < 0.6, "#4A6FA5", "#D63031")
)

density_A <- ggplot(density_df, aes(x = x, y = y)) +
    geom_area(aes(fill = fill_color), alpha = 0.7) +
    geom_line(aes(color = line_color), linewidth = 1.2) +
    scale_fill_identity() +
    scale_color_identity() +
    geom_vline(xintercept = 0.6, linetype = "dashed", 
               color = "black", linewidth = 1) +
    labs(x = "Anomaly Score", y = "Density") +
    theme_minimal() +
    theme(
        strip.background = element_rect(fill = "white", 
                                        color = "black", linewidth = 0.5),
        strip.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = 1)
    )

ggsave("figures/simulation/cell_type_A_anomaly_density_original.png", 
       density_A, 
       dpi = 600, 
       width = 6, height = 4)

# Add density plot for Cell Type B anomaly scores
anomaly_scores_B <- anomaly_output$`Cell Type B`$query_anomaly_scores
df_scores_B <- data.frame(
    anomaly_score = anomaly_scores_B,
    cell_type = "Cell Type B"
)

# Calculate density
density_calc <- density(df_scores_B$anomaly_score)
density_df <- data.frame(
    x = density_calc$x,
    y = density_calc$y,
    fill_color = ifelse(density_calc$x < 0.6, "#A1D99B", "#FF6B6B"),
    line_color = ifelse(density_calc$x < 0.6, "#2B8C3F", "#D63031")
)

density_B <- ggplot(density_df, aes(x = x, y = y)) +
    geom_area(aes(fill = fill_color), alpha = 0.7) +
    geom_line(aes(color = line_color), linewidth = 1.2) +
    scale_fill_identity() +
    scale_color_identity() +
    geom_vline(xintercept = 0.6, linetype = "dashed", 
               color = "black", linewidth = 1) +
    labs(x = "Anomaly Score", y = "Density") +
    theme_minimal() +
    theme(
        strip.background = element_rect(fill = "white", 
                                        color = "black", linewidth = 0.5),
        strip.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = 1)
    )

ggsave("figures/simulation/cell_type_B_anomaly_density_original.png", 
       density_B, 
       dpi = 600, 
       width = 6, height = 4)

# __________________________________________
# Missing Cell Type - Isolation Forest Plot
# __________________________________________

# Query data with original misannotations
query_data_missing_if <- query_data_missing
query_data_missing_if$SingleR_annotation[
    grepl("Misannotated as Cell Type A", query_data_missing_if$SingleR_annotation) 
] <- "Cell Type A"
query_data_missing_if$SingleR_annotation[
    grepl("Misannotated as Cell Type B", query_data_missing_if$SingleR_annotation) 
] <- "Cell Type B"

# Isolation Forest (same cell types)
anomaly_output <- detectAnomaly(reference_data = reference_data_missing,
                                query_data = query_data_missing_if,
                                ref_cell_type_col = "Cell_Type",
                                query_cell_type_col = "SingleR_annotation",
                                pc_subset = 1:2,
                                n_tree = 1000,
                                anomaly_threshold = 0.6)

# Plot the output for cell type A
gradient_colors <- c(
    "#7B9FC9", "#7B9FC9", "#7B9FC9",   
    "#9EB7D9", "#BAD0EB", "#CEE4F5",   
    "#E2F1FB", "white", "white")
stops <- c(
    0.00, 0.1, 0.2,
    0.3, 0.375, 0.425, 0.525,
    0.6, 0.65, 1.00
)
plotIsoForest_custom(
    anomaly_output,
    cell_type = "Cell Type A",
    pc_subset = 1:2,
    data_type = "query",
    n_tree = 1000,
    diagonal_facet = "blank",
    upper_facet = "blank",
    tile_colors = gradient_colors, 
    tile_stops = stops)
ggsave("figures/simulation/cell_type_A_isoForest.png", 
    dpi = 600, 
    width = 6, height = 4)

# Plot the output for cell type B
gradient_colors <- colorRampPalette(c(
    "#A1D99B", "#A1D99B", "#A1D99B",  
    "#A1D99B", "#B4E0AC", "#BDE5B3", 
    "#C7E9C0",  "white", "white"))(9)
stops <- c(
    0.00, 0.1, 0.2,
    0.3, 0.375, 0.425, 0.525,
    0.6, 0.65, 1.00
)
plotIsoForest_custom(
    anomaly_output,
    cell_type = "Cell Type B",
    pc_subset = 1:2,
    data_type = "query",
    n_tree = 1000,
    diagonal_facet = "blank",
    upper_facet = "blank",
    tile_colors = gradient_colors, 
    tile_stops = stops)
ggsave("figures/simulation/cell_type_B_isoForest.png", 
       dpi = 600, 
       width = 6, height = 4)


# Add density plot for Cell Type A anomaly scores (missing cell type scenario)
anomaly_scores_A_missing <- anomaly_output$`Cell Type A`$query_anomaly_scores
df_scores_A_missing <- data.frame(
    anomaly_score = anomaly_scores_A_missing,
    cell_type = "Cell Type A"
)

# Calculate density
density_calc <- density(df_scores_A_missing$anomaly_score)
density_df <- data.frame(
    x = density_calc$x,
    y = density_calc$y,
    fill_color = ifelse(density_calc$x < 0.6, "#7B9FC9", "#FF6B6B"),
    line_color = ifelse(density_calc$x < 0.6, "#4A6FA5", "#D63031")
)

density_A_missing <- ggplot(density_df, aes(x = x, y = y)) +
    geom_area(aes(fill = fill_color), alpha = 0.7) +
    geom_line(aes(color = line_color), linewidth = 1.2) +
    scale_fill_identity() +
    scale_color_identity() +
    geom_vline(xintercept = 0.6, linetype = "dashed", 
               color = "black", linewidth = 1) +
    labs(x = "Anomaly Score", y = "Density") +
    theme_minimal() +
    theme(
        strip.background = element_rect(fill = "white", 
                                        color = "black", linewidth = 0.5),
        strip.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = 1)
    )

ggsave("figures/simulation/cell_type_A_anomaly_density_missing.png", 
       density_A_missing, 
       dpi = 600,
       width = 6,  height = 4)

# Add density plot for Cell Type B anomaly scores (missing cell type scenario)
anomaly_scores_B_missing <- anomaly_output$`Cell Type B`$query_anomaly_scores
df_scores_B_missing <- data.frame(
    anomaly_score = anomaly_scores_B_missing,
    cell_type = "Cell Type B"
)

# Calculate density
density_calc <- density(df_scores_B_missing$anomaly_score)
density_df <- data.frame(
    x = density_calc$x,
    y = density_calc$y,
    fill_color = ifelse(density_calc$x < 0.6, "#A1D99B", "#FF6B6B")
)

density_B_missing <- ggplot(density_df, aes(x = x, y = y)) +
    geom_area(aes(fill = fill_color), alpha = 0.7) +
    scale_fill_identity() +
    geom_line(color = "black", linewidth = 1.2) +
    geom_vline(xintercept = 0.6, linetype = "dashed", 
               color = "black", linewidth = 1) +
    labs(x = "Anomaly Score", y = "Density") +
    theme_minimal() +
    theme(
        strip.background = element_rect(fill = "white", 
                                        color = "black", linewidth = 0.5),
        strip.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = 1)
    )


ggsave("figures/simulation/cell_type_B_anomaly_density_missing.png", 
       density_B_missing,
       dpi = 600, 
       width = 6, height = 4)
