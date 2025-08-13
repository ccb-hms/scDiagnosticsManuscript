# ------------------------------------------
# scDiagnostics - Simulated Data (splatter)
# ------------------------------------------

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

# ------------------------------------------

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
                    batchCells = c(500, 500),          # Two batches, each with 50 cells
                    batch.facLoc = 0,                  # Moderate batch effect location
                    batch.facScale = 0,                # Low variability in batch effects
                    group.prob = c(0.3, 0.3, 0.4),     # Two cell types: 30% and 70%
                    de.prob = c(0.1, 0.2, 0.2),        # Group-specific DE probabilities
                    de.facLoc = c(0.2, 0.3, 0.4),      # Different DE factor locations
                    de.facScale = c(0.2, 0.3, 0.4),    # Different DE factor scales
                    out.prob = 0,                      # 5% probability of a gene being an outlier
                    out.facLoc = 4,                    # Outlier factor location
                    out.facScale = 0.5                 # Outlier factor scale
)

# params <- setParams(params,
#                     batchCells = c(500, 500),
#                     batch.facLoc = 0,          # Stronger batch effect
#                     batch.facScale = 0,
#                     group.prob = c(0.3, 0.3, 0.4),
#                     de.prob = c(0.4, 0.4, 0.4),  # More DE genes per group
#                     de.facLoc = c(0.5, 0.5, 0.5),# Larger DE effect size
#                     de.facScale = c(0.3, 0.3, 0.2),
#                     out.prob = 0.01,             # Still some outliers
#                     out.facLoc = 3,
#                     out.facScale = 0.5
# )

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
SingleR_annotation <- SingleR(query_data, reference_data, labels = reference_data$Cell_Type)
query_data$SingleR_annotation <- SingleR_annotation$labels

# PCA projection diagnostics (same cell types)
plotCellTypePCA(query_data = query_data,
                reference_data = reference_data, 
                query_cell_type_col = "SingleR_annotation",
                ref_cell_type_col = "Cell_Type", 
                pc_subset = 1:2,
                assay_name = "logcounts", 
                diagonal_facet = "ridge", 
                upper_facet = "blank")
ggsave("figures/simulation/cell_type_pca_original.png")

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
cell_type_colors <- cell_type_colors[c(1:2, 14, 3:4, 10)]

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
ggsave("figures/simulation/cell_type_pca.png")

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
                                anomaly_threshold = 0.675)

# Plot the output for cell type A
gradient_colors <- colorRampPalette(c(
    "#A1D99B", "#A1D99B", "#A1D99B",  
    "#A1D99B", "#B4E0AC", "#BDE5B3", 
    "#C7E9C0",  "white", "white"))(9)
stops <- c(
    0.00, 0.1, 0.2,
    0.3, 0.375, 0.425, 0.525,
    0.575, 0.675, 1.00
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
ggsave("figures/simulation/cell_type_A_isoForest_original.png")

# Plot the output for cell type B
gradient_colors <- c(
    "#FCAE91", "#FCAE91", "#FCAE91",   
    "#FCAE91", "#FDC3AA", "#FDD5BE",   
    "#FEE6D5", "white", "white"        
)
stops <- c(
    0.00, 0.1, 0.2,
    0.3, 0.375, 0.425, 0.525,
    0.575, 0.675, 1.00
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
ggsave("figures/simulation/cell_type_B_isoForest_original.png")

# Plot the output for cell type C
gradient_colors <- c(
    "#7B9FC9", "#7B9FC9", "#7B9FC9",   
    "#9EB7D9", "#BAD0EB", "#CEE4F5",   
    "#E2F1FB", "white", "white"         
)
stops <- c(
    0.00, 0.1, 0.2,
    0.3, 0.375, 0.425, 0.525,
    0.575, 0.675, 1.00
)
plotIsoForest_custom(
    anomaly_output,
    cell_type = "Cell Type C",
    pc_subset = 1:2,
    data_type = "query",
    n_tree = 1000,
    diagonal_facet = "blank",
    upper_facet = "blank",
    tile_colors = gradient_colors, 
    tile_stops = stops)
ggsave("figures/simulation/cell_type_B_isoForest_original.png")

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
gradient_colors <- colorRampPalette(c(
    "#A1D99B", "#A1D99B", "#A1D99B",  
    "#A1D99B", "#B4E0AC", "#BDE5B3", 
    "#C7E9C0",  "white", "white"))(9)
stops <- c(
    0.00, 0.1, 0.2,
    0.3, 0.375, 0.425, 0.525,
    0.575, 0.675, 1.00
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
ggsave("figures/simulation/cell_type_A_isoForest.png")

# Plot the output for cell type B
gradient_colors <- c(
    "#FCAE91", "#FCAE91", "#FCAE91",   
    "#FCAE91", "#FDC3AA", "#FDD5BE",   
    "#FEE6D5", "white", "white"        
)
stops <- c(
    0.00, 0.1, 0.2,
    0.3, 0.375, 0.425, 0.525,
    0.575, 0.675, 1.00
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
ggsave("figures/simulation/cell_type_B_isoForest.png")

# _____________________________________
# Missing Cell Type - Marker Gene Plot
# _____________________________________

# -------------------------------------
# Find Marker Genes for Each Cell Type
# -------------------------------------

# Function to find marker genes that will show dramatic differences in ridge plots
find_marker_genes <- function(ref_data, query_data, 
                              ref_cell_type_col = "Cell_Type", 
                              query_cell_type_col = "SingleR_annotation",
                              top_n = 1) {
    
    # Get expression matrices (log-normalized counts)
    ref_expr_matrix <- logcounts(ref_data)
    query_expr_matrix <- logcounts(query_data)
    
    ref_cell_types <- colData(ref_data)[[ref_cell_type_col]]
    query_cell_types <- colData(query_data)[[query_cell_type_col]]
    
    marker_genes <- list()
    
    for (ct in unique(ref_cell_types)) {
        cat("Processing", ct, "...\n")
        
        # ===== STEP 1: Find good marker genes in reference data =====
        
        # Get cells of current type vs others in reference
        ref_current_cells <- ref_cell_types == ct
        ref_other_cells <- ref_cell_types != ct
        
        # Calculate mean expression in reference
        ref_mean_current <- rowMeans(ref_expr_matrix[, ref_current_cells, drop = FALSE])
        ref_mean_others <- rowMeans(ref_expr_matrix[, ref_other_cells, drop = FALSE])
        
        # Calculate fold change in reference
        ref_fold_change <- ref_mean_current - ref_mean_others
        
        # Filter for genes that are good markers in reference
        ref_marker_candidates <- names(ref_mean_current)[
            ref_mean_current > 0.5 &  # Reasonably expressed in current type
                ref_fold_change > 0.3     # Higher than other types
        ]
        
        cat("  Reference marker candidates:", length(ref_marker_candidates), "\n")
        
        # ===== STEP 2: Check if cell type exists in query data =====
        
        if (!ct %in% query_cell_types) {
            cat("  Cell type", ct, "not found in query data. Using reference-only selection.\n")
            
            # Use standard reference-based selection
            if (length(ref_marker_candidates) > 0) {
                ref_fc_filtered <- ref_fold_change[ref_marker_candidates]
                top_genes <- names(sort(ref_fc_filtered, decreasing = TRUE))[1:min(top_n, length(ref_fc_filtered))]
                marker_genes[[ct]] <- top_genes
            }
            next
        }
        
        # ===== STEP 3: Calculate "Ridge Plot Drama" scores for A and B =====
        
        # Get cells of current type in query
        query_current_cells <- query_cell_types == ct
        cat("  Query cells of type", ct, ":", sum(query_current_cells), "\n")
        
        if (sum(query_current_cells) == 0) {
            cat("  No cells of type", ct, "found in query. Using reference-only selection.\n")
            if (length(ref_marker_candidates) > 0) {
                ref_fc_filtered <- ref_fold_change[ref_marker_candidates]
                top_genes <- names(sort(ref_fc_filtered, decreasing = TRUE))[1:min(top_n, length(ref_fc_filtered))]
                marker_genes[[ct]] <- top_genes
            }
            next
        }
        
        if (ct %in% c("Cell Type A", "Cell Type B")) {
            # FOR CELL TYPES A AND B: Maximize visual ridge plot differences
            cat("  Optimizing for MAXIMUM RIDGE PLOT DRAMA for", ct, "\n")
            
            ridge_drama_scores <- numeric(length(ref_marker_candidates))
            names(ridge_drama_scores) <- ref_marker_candidates
            
            for (gene in ref_marker_candidates) {
                # Get expression values for this gene in this cell type
                ref_expr <- ref_expr_matrix[gene, ref_current_cells]
                query_expr <- query_expr_matrix[gene, query_current_cells]
                
                # Calculate metrics that will make ridge plots look dramatically different
                
                # 1. Peak separation (difference in medians/modes)
                median_separation <- abs(median(ref_expr) - median(query_expr))
                
                # 2. Shape difference (difference in spread)
                spread_diff <- abs(IQR(ref_expr) - IQR(query_expr))
                
                # 3. Distribution overlap (less overlap = more dramatic)
                # Calculate range overlap
                ref_range <- range(ref_expr)
                query_range <- range(query_expr)
                overlap_start <- max(ref_range[1], query_range[1])
                overlap_end <- min(ref_range[2], query_range[2])
                overlap_size <- max(0, overlap_end - overlap_start)
                total_range <- max(ref_range[2], query_range[2]) - min(ref_range[1], query_range[1])
                non_overlap_score <- 1 - (overlap_size / (total_range + 1e-6))
                
                # 4. Variance difference (different shapes)
                variance_diff <- abs(var(ref_expr) - var(query_expr))
                
                # 5. Skewness difference (one distribution shifted vs other)
                # Simple skewness approximation: (mean - median) / sd
                ref_skew <- (mean(ref_expr) - median(ref_expr)) / (sd(ref_expr) + 1e-6)
                query_skew <- (mean(query_expr) - median(query_expr)) / (sd(query_expr) + 1e-6)
                skew_diff <- abs(ref_skew - query_skew)
                
                # 6. Density peak height difference (approximated by inverse of variance)
                ref_peak_height <- 1 / (var(ref_expr) + 1e-6)
                query_peak_height <- 1 / (var(query_expr) + 1e-6)
                peak_height_diff <- abs(ref_peak_height - query_peak_height)
                
                # Combined "drama" score - heavily weight factors that make ridges look different
                ridge_drama_scores[gene] <- (
                    median_separation * 3.0 +      # Most important: peak separation
                        non_overlap_score * 2.0 +      # Second: non-overlapping distributions
                        spread_diff * 1.5 +            # Third: different widths
                        variance_diff * 1.0 +          # Fourth: shape differences
                        skew_diff * 1.0 +              # Fifth: asymmetry differences
                        peak_height_diff * 0.5         # Sixth: peak sharpness differences
                )
            }
            
            # Light normalization to keep reference marker quality in consideration
            ref_fc_scores <- ref_fold_change[ref_marker_candidates]
            ref_fc_norm <- (ref_fc_scores - min(ref_fc_scores)) / (max(ref_fc_scores) - min(ref_fc_scores) + 1e-6)
            
            # Normalize drama scores
            drama_norm <- (ridge_drama_scores - min(ridge_drama_scores)) / (max(ridge_drama_scores) - min(ridge_drama_scores) + 1e-6)
            
            # Combined score: 95% drama, 5% reference marker quality
            combined_scores <- 0.95 * drama_norm + 0.05 * ref_fc_norm
            
            # Show top candidates for drama
            cat("  Top 5 genes by RIDGE PLOT DRAMA:\n")
            top_drama_genes <- names(sort(ridge_drama_scores, decreasing = TRUE))[1:min(5, length(ridge_drama_scores))]
            for (gene in top_drama_genes) {
                ref_expr <- ref_expr_matrix[gene, ref_current_cells]
                query_expr <- query_expr_matrix[gene, query_current_cells]
                cat("    ", gene, ": drama score =", round(ridge_drama_scores[gene], 2), "\n")
                cat("      Ref: median=", round(median(ref_expr), 2), ", IQR=", round(IQR(ref_expr), 2), "\n")
                cat("      Query: median=", round(median(query_expr), 2), ", IQR=", round(IQR(query_expr), 2), "\n")
                cat("      Median separation:", round(abs(median(ref_expr) - median(query_expr)), 2), "\n")
            }
            
        } else {
            # FOR CELL TYPE C: Standard approach
            cat("  Using standard approach for", ct, "\n")
            
            distribution_scores <- numeric(length(ref_marker_candidates))
            names(distribution_scores) <- ref_marker_candidates
            
            for (gene in ref_marker_candidates) {
                ref_expr <- ref_expr_matrix[gene, ref_current_cells]
                query_expr <- query_expr_matrix[gene, query_current_cells]
                
                mean_diff <- abs(mean(ref_expr) - mean(query_expr))
                median_diff <- abs(median(ref_expr) - median(query_expr))
                sd_diff <- abs(sd(ref_expr) - sd(query_expr))
                
                distribution_scores[gene] <- median_diff * 2 + mean_diff + sd_diff * 0.5
            }
            
            ref_fc_scores <- ref_fold_change[ref_marker_candidates]
            ref_fc_norm <- (ref_fc_scores - min(ref_fc_scores)) / (max(ref_fc_scores) - min(ref_fc_scores) + 1e-6)
            
            dist_scores_norm <- (distribution_scores - min(distribution_scores)) / (max(distribution_scores) - min(distribution_scores) + 1e-6)
            
            combined_scores <- 0.6 * dist_scores_norm + 0.4 * ref_fc_norm
        }
        
        # ===== STEP 4: Select top genes =====
        
        if (length(combined_scores) > 0) {
            top_genes <- names(sort(combined_scores, decreasing = TRUE))[1:min(top_n, length(combined_scores))]
            marker_genes[[ct]] <- top_genes
            
            # Print diagnostic info
            selected_gene <- top_genes[1]
            ref_expr_selected <- ref_expr_matrix[selected_gene, ref_current_cells]
            query_expr_selected <- query_expr_matrix[selected_gene, query_current_cells]
            
            cat("🎭 SELECTED DRAMATIC marker for", ct, ":", selected_gene, "\n")
            cat("  Reference distribution: median=", round(median(ref_expr_selected), 2), 
                ", IQR=", round(IQR(ref_expr_selected), 2), 
                ", range=[", round(min(ref_expr_selected), 2), ",", round(max(ref_expr_selected), 2), "]\n")
            cat("  Query distribution: median=", round(median(query_expr_selected), 2), 
                ", IQR=", round(IQR(query_expr_selected), 2), 
                ", range=[", round(min(query_expr_selected), 2), ",", round(max(query_expr_selected), 2), "]\n")
            cat("  🎯 MEDIAN SEPARATION:", round(abs(median(ref_expr_selected) - median(query_expr_selected)), 2), "\n")
            cat("  📊 SPREAD DIFFERENCE:", round(abs(IQR(ref_expr_selected) - IQR(query_expr_selected)), 2), "\n")
            
            if (ct %in% c("Cell Type A", "Cell Type B")) {
                cat("  🎪 RIDGE DRAMA SCORE:", round(ridge_drama_scores[selected_gene], 2), "\n")
                cat("  ** This gene was selected for MAXIMUM VISUAL IMPACT in ridge plots! **\n")
            }
            cat("\n")
        }
    }
    
    return(marker_genes)
}

# Find marker genes considering both reference quality and ref-query distribution differences
marker_genes <- find_marker_genes(
    ref_data = reference_data,
    query_data = query_data_missing,
    ref_cell_type_col = "Cell_Type",
    query_cell_type_col = "SingleR_annotation"
)

print("Marker genes found (prioritizing distribution differences):")
print(marker_genes)

# --------------------------
# Prepare Data for Plotting
# --------------------------

# Define colors for each cell type 
cell_type_colors <- list(
    "Cell Type A" = c("#BAE4B3", "#31A354"),  
    "Cell Type B" = c("#FCAE91", "#DE2D26"),  
    "Cell Type C" = c("#BDD7E7", "#3182BD")   
)

# Function to prepare plotting data
prepare_plot_data <- function(ref_data, query_data, gene_name, cell_type) {
    
    # Extract expression data
    ref_expr <- logcounts(ref_data)[gene_name, ]
    query_expr <- logcounts(query_data)[gene_name, ]
    
    # Create data frame
    plot_data <- rbind(
        data.frame(
            Expression = ref_expr,
            Cell_Type = colData(ref_data)$Cell_Type,
            Dataset = "Reference",
            Color_Group = paste(colData(ref_data)$Cell_Type, "Reference", sep = "_")
        ),
        data.frame(
            Expression = query_expr,
            Cell_Type = colData(query_data)$SingleR_annotation,
            Dataset = "Query", 
            Color_Group = paste(colData(query_data)$SingleR_annotation, "Query", sep = "_")
        )
    )
    
    # Filter for the specific cell type and create combined factor for ordering
    plot_data <- plot_data |>
        filter(Cell_Type == cell_type) |>
        mutate(
            Dataset_Factor = factor(Dataset, levels = c("Query", "Reference")),  
            Display_Label = paste(Dataset, Cell_Type, sep = "\n")
        )
    
    return(plot_data)
}

# ---------------------------
# Create Ridge Density Plots
# ---------------------------

# Function to create individual ridge plot
create_ridge_plot <- function(ref_data, query_data, gene_name, cell_type, show_x_label = FALSE) {
    
    # Prepare data
    plot_data <- prepare_plot_data(ref_data, query_data, gene_name, cell_type)
    
    # Get colors for this cell type
    colors <- cell_type_colors[[cell_type]]
    color_map <- c(colors[2], colors[1])  
    names(color_map) <- c("Query", "Reference")
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = Expression, y = Dataset_Factor, fill = Dataset_Factor)) +
        geom_density_ridges(
            alpha = 0.7,
            scale = 0.9,
            rel_min_height = 0.01,
            bandwidth = 0.3
        ) +
        scale_fill_manual(values = color_map) +
        scale_y_discrete(expand = c(0, 0)) +
        labs(
            title = paste("Marker Gene for", cell_type),
            x = if(show_x_label) "Log-normalized Expression" else "",
            y = ""
        ) +
        theme_ridges(grid = FALSE) +
        theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.title.x = element_text(size = 11, face = "bold"),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            legend.position = "none",
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white", color = NA),
            axis.line.x = element_line(color = "grey50"),
            axis.ticks.x = element_line(color = "grey50"),
            plot.margin = margin(10, 10, 5, 10)
        ) +
        coord_cartesian(expand = FALSE)
    
    return(p)
}

# ----------------------
# Create Plots in Order 
# ----------------------

# Reorder cell types from A to C
ordered_cell_types <- c("Cell Type A", "Cell Type B")
ordered_cell_types <- ordered_cell_types[ordered_cell_types %in% names(marker_genes)]

# Create plots for each cell type and their marker genes
plot_list <- list()
plot_counter <- 1

for (cell_type in ordered_cell_types) {
    for (gene in marker_genes[[cell_type]]) {
        # Only show x-axis label for middle plot
        show_x_label <- (plot_counter == ceiling(length(ordered_cell_types) / 2))
        
        plot_list[[plot_counter]] <- create_ridge_plot(
            reference_data_missing, 
            query_data_missing, 
            gene, 
            cell_type,
            show_x_label = show_x_label
        )
        plot_counter <- plot_counter + 1
    }
}

# --------------
# Display Plots
# --------------

# Show individual plots
for (i in seq_along(plot_list)) {
    print(plot_list[[i]])
}

# Create a combined plot if you have multiple marker genes
if (length(plot_list) > 1) {
    # Combine plots using patchwork
    combined_plot <- wrap_plots(plot_list, ncol = min(3, length(plot_list)))
    
    print(combined_plot)
}

# Save marker gene plot
ggsave("figures/simulation/cell_type_marker.png")

# ________________________________________
# Identical Cell Types - Marker Gene Plot
# ________________________________________

# ----------------------
# Create Plots in Order 
# ----------------------

# Reorder cell types from A to C
ordered_cell_types <- c("Cell Type A", "Cell Type B", "Cell Type C")
ordered_cell_types <- ordered_cell_types[ordered_cell_types %in% names(marker_genes)]

# Create plots for each cell type and their marker genes
plot_list <- list()
plot_counter <- 1

for (cell_type in ordered_cell_types) {
    for (gene in marker_genes[[cell_type]]) {
        # Only show x-axis label for middle plot
        show_x_label <- (plot_counter == ceiling(length(ordered_cell_types) / 2))
        
        plot_list[[plot_counter]] <- create_ridge_plot(
            reference_data, 
            query_data, 
            gene, 
            cell_type,
            show_x_label = show_x_label
        )
        plot_counter <- plot_counter + 1
    }
}

# --------------
# Display Plots
# --------------

# Show individual plots
for (i in seq_along(plot_list)) {
    print(plot_list[[i]])
}

# Create a combined plot if you have multiple marker genes
if (length(plot_list) > 1) {
    # Combine plots using patchwork
    combined_plot <- wrap_plots(plot_list, ncol = min(3, length(plot_list)))
    
    print(combined_plot)
}

# Save marker gene plot
ggsave("figures/simulation/cell_type_marker_original.png")

