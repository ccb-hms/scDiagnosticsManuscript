# -------------------------------
# COVID-19 PBMC - PCA IFN Plot
# -------------------------------

# Load libraries
library(SingleCellExperiment)
library(scDiagnostics)
library(ggplot2)
library(scran)
library(scater)
library(GGally)
library(ggridges)
library(viridis)

# -------------------------------

# Load the processed SCE objects
normal_data <- readRDS("data/covid/normal_data_sce.rds")
covid_data <- readRDS("data/covid/covid_data_sce.rds")

# ______________________
# Setup and Parameters
# ______________________

# Define the Yoshida IFN gene signature
yoshida_ifn_signature <- c(
  "BST2", "CMPK2", "EIF2AK2", "EPSTI1", "HERC5", "IFI35", "IFI44L", 
  "IFI6", "IFIT3", "ISG15", "LY6E", "MX1", "MX2", "OAS1", "OAS2", 
  "PARP9", "PLSCR1", "SAMD9", "SAMD9L", "SP110", "STAT1", "TRIM22", 
  "UBE2L6", "XAF1", "IRF7"
)

# Define plot parameters
CELL_TYPES_TO_PLOT <- c("CD14_mono")
PC_SUBSET <- 1:3
REF_CELL_TYPE_COL <- "author_cell_type"
QUERY_CELL_TYPE_COL <- "azimuth_celltype_l1"
ASSAY_NAME <- "logcounts"
MAX_CELLS <- 2000 # Max cells per dataset (ref/query)

# _____________________________________
# Project Data and Calculate IFN Score
# _____________________________________

# Project data using scDiagnostics::projectPCA
pca_output <- scDiagnostics::projectPCA(
    query_data = covid_data,
    reference_data = normal_data,
    query_cell_type_col = QUERY_CELL_TYPE_COL,
    ref_cell_type_col = REF_CELL_TYPE_COL,
    cell_types = CELL_TYPES_TO_PLOT,
    pc_subset = PC_SUBSET,
    assay_name = ASSAY_NAME,
    max_cells_query = MAX_CELLS,
    max_cells_ref = MAX_CELLS
)

# Extract original cell names and add to the data frame
pca_output$original_cell_id <- gsub("Reference_|Query_", "", rownames(pca_output))

# Find available signature genes
common_genes <- intersect(rownames(covid_data), rownames(normal_data))
signature_genes_avail <- intersect(yoshida_ifn_signature, common_genes)
if (length(signature_genes_avail) < length(yoshida_ifn_signature)) {
    warning(paste(length(yoshida_ifn_signature) - length(signature_genes_avail), "IFN genes not found."))
}

# Create a combined expression matrix for the plotted cells
ref_cells_toplot <- pca_output$original_cell_id[pca_output$dataset == "Reference"]
query_cells_toplot <- pca_output$original_cell_id[pca_output$dataset == "Query"]
full_expr_matrix <- cbind(
    assay(normal_data, ASSAY_NAME)[signature_genes_avail, ref_cells_toplot, drop = FALSE],
    assay(covid_data, ASSAY_NAME)[signature_genes_avail, query_cells_toplot, drop = FALSE]
)

# Calculate and add IFN score
pca_output$ifn_score <- colMeans(full_expr_matrix, na.rm = TRUE)[pca_output$original_cell_id]

# Set dataset factor levels to control ridge plot order (first level is at the bottom)
pca_output$dataset <- factor(pca_output$dataset, levels = c("Query", "Reference"))

# _______________________________________
# Define Final Custom Plotting Functions
# _______________________________________

# Lower panels: Scatter plot with SOLID (Reference) vs. HOLLOW (Query) shapes
scatter_fn <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
        geom_point(data = ~subset(., dataset == "Query"), alpha = 0.5, size = 1.5, show.legend = FALSE) +
        geom_point(data = ~subset(., dataset == "Reference"), alpha = 0.7, size = 1.5, show.legend = FALSE) +
        scale_shape_manual(values = c("Reference" = 16, "Query" = 1)) +
        viridis::scale_color_viridis(option = "B") +
        theme_minimal() +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
}

# Diagonal panels: Ridge plot with specified colors and no Y-axis labels
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

# Upper panels: Blank
blank_fn <- function(data, mapping, ...) {
    ggplot() + theme_void() +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
}

# ________________________
# Create a Legend Object
# ________________________

# Create a simple "bogus" plot with the desired aesthetics to generate the legend
legend_plot <- ggplot(pca_output, aes(x = PC1, y = PC2, color = ifn_score, shape = dataset)) +
    geom_point() +
    scale_shape_manual(name = "Dataset", values = c("Reference" = 16, "Query" = 1)) +
    viridis::scale_color_viridis(option = "B", name = "IFN Signature\nScore") +
    theme(
        legend.position = "right", 
        legend.box = "vertical",
        legend.key = element_rect(fill = "white", color = NA)
    )

# Grab the legend from this plot
plot_legend <- GGally::grab_legend(legend_plot)

# __________________________________
# Create and Display the Final Plot
# __________________________________

# Set up column names with variance explained
pc_plot_names <- paste0(
    "PC", PC_SUBSET, " (",
    sprintf("%.1f%%", attributes(reducedDim(normal_data, "PCA"))[["percentVar"]][PC_SUBSET]), ")"
)

# Create the pairs plot, now explicitly passing the legend object
pca_ifn_plot <- GGally::ggpairs(
    pca_output,
    columns = paste0("PC", PC_SUBSET),
    columnLabels = pc_plot_names,
    # Mapping is still needed for the functions to work, but legend is now separate
    mapping = aes(color = ifn_score, shape = dataset),
    lower = list(continuous = scatter_fn),
    diag = list(continuous = ridge_fn),
    upper = list(continuous = blank_fn),
    progress = FALSE,
    legend = plot_legend
)

# Apply final theme for publication quality
pca_ifn_plot <- pca_ifn_plot + theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(color = "black")
)

# Print the final plot
print(pca_ifn_plot)

# Save the final plot
ggsave("figures/covid/pca_ifn.png", width = 12, height = 8, dpi = 600)
