# ---------------------------------------------------------------
# MERFISH Mouse Colon IBD - PCA Plot for Fibroblast ECM Score
# ---------------------------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scDiagnostics)
library(ggplot2)
library(GGally)
library(ggridges)
library(viridis)

# ---------------------------------------------------------------

# ________________________
# Data Loading & Lumping
# ________________________

cat("\n--- Loading MERFISH Datasets ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

healthy_data$tier2 <- as.character(healthy_data$tier2)
dss9_data$tier2 <- as.character(dss9_data$tier2)
set.seed(0)

cat("\n--- Grouping Fibroblast Lineage ---\n")
# Regex to catch all fibroblast types (Healthy, IAF)
is_fibro_lineage <- function(x) {
    grepl("^Fibro|^IAF", x)
}

# Apply lumping label
healthy_data$analysis_class <- ifelse(is_fibro_lineage(healthy_data$tier2), "Fibroblast_Lineage", "Other")
dss9_data$analysis_class <- ifelse(is_fibro_lineage(dss9_data$tier2), "Fibroblast_Lineage", "Other")

cat(sprintf("Healthy Fibroblasts: %d\n", sum(healthy_data$analysis_class == "Fibroblast_Lineage")))
cat(sprintf("DSS9 Fibroblasts:    %d\n", sum(dss9_data$analysis_class == "Fibroblast_Lineage")))

# ______________________
# Parameters & Analysis
# ______________________

ecm_homeostasis_signature <- c("Col1a2", "Timp2", "Col6a1", "Sparc", "Dpt")

# Parameters
CELL_TYPE_TO_PLOT <- "Fibroblast_Lineage" 
PC_SUBSET <- 1:4
REF_CELL_TYPE_COL <- "analysis_class"
QUERY_CELL_TYPE_COL <- "analysis_class"
ASSAY_NAME <- "logcounts"
MAX_CELLS <- 1000

cat("--- Projecting data and calculating signature scores ---\n")

# Project PCA
pca_output <- scDiagnostics::projectPCA(
    query_data = dss9_data,
    reference_data = healthy_data,
    query_cell_type_col = QUERY_CELL_TYPE_COL,
    ref_cell_type_col = REF_CELL_TYPE_COL,
    cell_types = CELL_TYPE_TO_PLOT,
    pc_subset = PC_SUBSET,
    assay_name = ASSAY_NAME,
    max_cells_query = MAX_CELLS,
    max_cells_ref = MAX_CELLS
)

# Clean IDs and Calculate Score
pca_output$original_cell_id <- gsub("Reference_|Query_", "", rownames(pca_output))
common_genes <- intersect(rownames(dss9_data), rownames(healthy_data))
signature_genes_avail <- intersect(ecm_homeostasis_signature, common_genes)

# Extract expression for score
ref_cells_toplot <- pca_output$original_cell_id[pca_output$dataset == "Reference"]
query_cells_toplot <- pca_output$original_cell_id[pca_output$dataset == "Query"]

full_expr_matrix <- cbind(
    assay(healthy_data, ASSAY_NAME)[signature_genes_avail, ref_cells_toplot, drop = FALSE],
    assay(dss9_data, ASSAY_NAME)[signature_genes_avail, query_cells_toplot, drop = FALSE]
)

pca_output$ecm_score <- colMeans(full_expr_matrix, na.rm = TRUE)[match(pca_output$original_cell_id, colnames(full_expr_matrix))]
pca_output$dataset <- factor(pca_output$dataset, levels = c("Query", "Reference"))

# __________________________
# Custom Plotting Functions 
# __________________________

# BOTH Reference and Query are colored by ECM Score
scatter_fn <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
        # 1. Plot Reference (Open Circle) - Colored by Score
        geom_point(data = ~subset(., dataset == "Reference"), 
                   aes(color = ecm_score),
                   shape = 1,          
                   alpha = 0.5, 
                   size = 0.75) +
        
        # 2. Plot Query (Filled Point) - Colored by Score
        geom_point(data = ~subset(., dataset == "Query"), 
                   aes(color = ecm_score), 
                   shape = 16,         
                   alpha = 0.8, 
                   size = 0.75) +
        
        # Shared Gradient (Plasma)
        viridis::scale_color_viridis(option = "plasma", direction = -1) + 
        theme_minimal() +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
}

# This keeps the datasets distinct in the histograms
ridge_fn <- function(data, mapping, ...) {
    ggplot(data = data, mapping = aes(x = !!mapping$x, y = dataset, fill = dataset, color = dataset)) +
        ggridges::geom_density_ridges(alpha = 0.6, scale = 1.5, rel_min_height = 0.01) +
        scale_fill_manual(values = c("Reference" = "#5A9BD8", "Query" = "#B565D8"), guide = "none") +
        scale_color_manual(values = c("Reference" = "#5A9BD8", "Query" = "#B565D8"), guide = "none") +
        theme_minimal() +
        theme(
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
}

blank_fn <- function(data, mapping, ...) {
    ggplot() + theme_void() +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
}

# Shows the Colorbar for ECM Score
legend_plot <- ggplot(pca_output, aes(x = PC1, y = PC2, color = ecm_score)) +
    geom_point() +
    viridis::scale_color_viridis(option = "plasma", direction = -1, name = "ECM Homeostasis\nScore") +
    theme(
        legend.position = "right", 
        legend.box = "vertical",
        legend.key = element_rect(fill = NA, color = NA) 
    )
plot_legend <- GGally::grab_legend(legend_plot)

# ________________
# Generate & Save
# ________________

cat("\n--- Generating final pairs plot ---\n")

pc_plot_names <- paste0(
    "PC", PC_SUBSET, " (",
    sprintf("%.1f%%", attributes(reducedDim(healthy_data, "PCA"))[["percentVar"]][PC_SUBSET]), ")"
)

pca_ecm_plot <- GGally::ggpairs(
    pca_output,
    columns = paste0("PC", PC_SUBSET),
    columnLabels = pc_plot_names,
    mapping = aes(color = ecm_score, shape = dataset), 
    lower = list(continuous = scatter_fn),
    diag = list(continuous = ridge_fn),
    upper = list(continuous = blank_fn),
    progress = FALSE,
    legend = plot_legend
) + 
theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(color = "black", face = "bold")
)

print(pca_ecm_plot)

dir.create("figures/merfish", showWarnings = FALSE, recursive = TRUE)
ggsave("figures/merfish/pca_fibroblast_lineage_ecm_corrected.png", pca_ecm_plot, width = 12, height = 8, dpi = 600)
cat("✓ Plot saved.\n")