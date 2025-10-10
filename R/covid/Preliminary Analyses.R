# ----------------------------------
# COVID-19 - Preliminary Analyses
# ----------------------------------

# Load libraries
library(SingleCellExperiment)
library(scDiagnostics)
library(ggplot2)
library(scran)
library(scater)
library(ComplexHeatmap)
library(circlize)

# ----------------------------------

# Load the processed SCE objects
normal_data <- readRDS("data/covid/normal_data_sce.rds")
covid_data <- readRDS("data/covid/covid_data_sce.rds")

# _________________________
# PCA Projections Boxplots
# _________________________

# PCA boxplot
boxplot_pca <- boxplotPCA(query_data = covid_data,
                          reference_data = normal_data, 
                          query_cell_type_col = "azimuth_celltype_l1", 
                          ref_cell_type_col = "author_cell_type", 
                          cell_types = c("CD14_mono", "CD83_CD14_mono"), 
                          pc_subset = 1:9,
                          max_cells_ref = 5000,
                          max_cells_query = 5000)

boxplot_pca

# _________________
# PCA Projections
# _________________

# PCA scatterplot
scatter_pca <- plotCellTypePCA(query_data  = covid_data,
                               reference_data = normal_data, 
                               query_cell_type_col = "azimuth_celltype_l1", 
                               ref_cell_type_col =  "author_cell_type", 
                               cell_types = c("CD83_CD14_mono"), 
                               pc_subset = 1:3, 
                               upper_facet = "blank",
                               diagonal_facet = "ridge", 
                               max_cells_query = 500, 
                               max_cells_ref = 500)

scatter_pca
ggsave("cell_type_pca.png", dpi = 600)

# __________________
# Anomaly Detection
# __________________

# Anomaly Detection
anomaly_output <- detectAnomaly(query_data = covid_data,
                                reference_data = normal_data, 
                                query_cell_type_col = "azimuth_celltype_l1", 
                                ref_cell_type_col = "author_cell_type", 
                                cell_types = c("CD14_mono", "CD83_CD14_mono"), 
                                pc_subset = 1:10,
                                n_tree = 500,
                                anomaly_threshold = 0.5, 
                                max_cells_query = NULL, 
                                max_cells_ref = NULL)

anomaly_plot <- plot(anomaly_output,
                     cell_type = "CD83_CD14_mono",
                     pc_subset = 1:3,
                     data_type = "query",
                     n_tree = 500,
                     diagonal_facet = "ridge",
                     upper_facet = "blank", 
                     max_cells_query = 1000)

anomaly_plot

ggsave("figures/covid/anomaly_plot.png")

# _______________________________
# Top Genes Distributional Shift
# _______________________________

# Compute distributional shifts for genes with top loadings
gene_shifts <- calculateTopLoadingGeneShifts(
     query_data = covid_data,
     reference_data = normal_data,
     query_cell_type_col = "azimuth_celltype_l1",
     ref_cell_type_col = "author_cell_type",
     cell_types = c("CD14_mono", "CD83_CD14_mono"), 
     pc_subset = 1:5,
     n_top_loadings = 30,
     assay_name = "logcounts",
     p_value_threshold = 0.05,
     adjust_method = "fdr", 
     detect_anomalies = TRUE,
     anomaly_comparison = TRUE,
     anomaly_threshold = 0.5, 
     max_cells_query = 5000, 
     max_cells_ref = 5000)

png("heatmap.png", width = 1200, height = 1200, res = 300)
plot(gene_shifts,
     cell_type = "CD83_CD14_mono",
     pc_subset = 1:5,
     plot_type = c("heatmap", "barplot", "boxplot")[2],
     plot_by = c("p_adjusted", "top_loading")[1],
     n_genes = 5,
     significance_threshold = 0.05,
     show_anomalies = TRUE,
     pseudo_bulk = TRUE,
     cluster_cols = TRUE,
     max_cells_ref = 5000,
     max_cells_query = 5000)
dev.off()

ggsave("figures/covid/gene_shifts_heatmap.png")


# __________________
# Marker Expression
# __________________

# Plot the marker expression
plotMarkerExpression(reference_data = normal_data,
                     query_data = covid_data,
                     query_cell_type_col = "singler_annotations",
                     ref_cell_type_col = "cell_type",
                     gene_name = "IFITM2",
                     cell_type = "CD14-positive monocyte",
                     normalization = "min_max")
ggsave("figures/covid/marker_expression.png")

# Plot the marker expression (mild and asymptomatic)
plotMarkerExpression(reference_data = normal_data,
                     query_data = covid_data[, covid_data$Status_on_day_collection_summary %in%
                       c("Mild", "Asymptomatic")],
                     query_cell_type_col = "singler_annotations",
                     ref_cell_type_col = "cell_type",
                     gene_name = "IFITM2",
                     cell_type = "CD14-positive monocyte",
                     normalization = "min_max")
ggsave("figures/covid/marker_expression_mild_asymptomatic.png")

# Regress the PC data (with reference data)
regress_res2 <- regressPC(query_data = covid_data,
                          reference_data = normal_data, 
                          query_cell_type_col = "singler_annotations", 
                          ref_cell_type_col = "cell_type", 
                          query_batch_col = "Site",
                          cell_types = c("CD14-positive monocyte",
                            "naive B cell"), 
                          pc_subset = 1:10, 
                                max_cells_query = NULL, 
                                max_cells_ref = NULL)
plot(regress_res2, plot_type = "r_squared")
plot(regress_res2, plot_type = "variance_contribution")
plot(regress_res2, plot_type = "coefficient_heatmap")

# ______________________________
# Graph Integration Diagnostics
# ______________________________

# Graph integration diagnostics
graph_diagnostics <- calculateGraphIntegration(
    query_data = covid_data,
    reference_data = normal_data,
    query_cell_type_col = "singler_annotations",
    ref_cell_type_col = "cell_type",
    cell_types = c("CD14-positive monocyte"), 
    pc_subset = 1:10,
    k_neighbors = 30,
    resolution = 0.15,
    high_query_prop_threshold = 0.9,
    cross_type_threshold = 0.1,
    local_consistency_threshold = 0.6,
    local_confidence_threshold = 0.2, 
    max_cells_query = NULL, 
    max_cells_ref = NULL
)

# Network graph showing all issue types (color by cell type)
plot(graph_diagnostics, plot_type = "community_network",
     color_by = "cell_type")

# ______________________
# Wasserstein Distances
# ______________________

# Compute Wasserstein distances 
wasserstein_data <- calculateWassersteinDistance(
    query_data = covid_data,
    reference_data = normal_data,
    query_cell_type_col = "singler_annotations",
    ref_cell_type_col = "cell_type",
    cell_types = c("CD14-positive monocyte"), 
    pc_subset = 1:5,
    n_resamples = 500)
plot(wasserstein_data)

# ________________
# SIR Projections
# ________________

# Compute SIR projections
sir_output <- calculateSIRSpace(reference_data = normal_data,
                                query_data = covid_data,
                                query_cell_type_col = "singler_annotations",
                                ref_cell_type_col = "cell_type",
                                cell_types = c("CD14-positive monocyte"), 
                                multiple_cond_means = TRUE,
                                cumulative_variance_threshold = 0.9,
                                n_neighbor = 1)

plot(sir_output,
     sir_subset = 1:5,
     cell_types = c("CD14-positive monocyte"),
     lower_facet = "scatter",
     diagonal_facet = "ridge",
     upper_facet = "blank")

# _____________
# Unified UMAP
# _____________

# Create dataset column
colData(normal_data)$dataset <- "Reference"
colData(covid_data)$dataset <- "Query"

# Create unified cell_group column
colData(normal_data)$cell_group <- colData(normal_data)$cell_type
colData(covid_data)$cell_group <- colData(covid_data)$singler_annotations

# Create combined group_label column
colData(normal_data)$group_label <- paste("Reference", normal_data$cell_group)
colData(covid_data)$group_label <- paste("Query", covid_data$cell_group)

# Ensure column consistency
common_cols <- intersect(colnames(colData(normal_data)), colnames(colData(covid_data)))
colData(normal_data) <- colData(normal_data)[, common_cols]
colData(covid_data) <- colData(covid_data)[, common_cols]

# Ensure row consistency
common_genes <- intersect(rownames(normal_data), rownames(covid_data))
normal_data <- normal_data[common_genes, ]
covid_data <- covid_data[common_genes, ]

# Remove counts assay from covid_data
assay(covid_data, "counts") <- NULL

# Combine SCE objects
combined_sce <- cbind(normal_data, covid_data)

# Feature selection and batch correction with PCA
gene_variances <- modelGeneVar(combined_sce, block = combined_sce$dataset)
hvg <- getTopHVGs(gene_variances, n = 2000)
combined_sce <- runPCA(combined_sce, subset_row = hvg)

# Run UMAP on corrected PCs
combined_sce <- runUMAP(combined_sce, dimred = "PCA")


# Create UMAP plot
umap_plot <- plotUMAP(
  combined_sce[, combined_sce$group_label %in%
    c("Reference CD14-positive monocyte",
      "Query CD14-positive monocyte",
      "Reference naive B cell", 
      "Query naive B cell", 
      "Reference effector CD8-positive, alpha-beta T cell",
      "Query effector CD8-positive, alpha-beta T cell",
      "Reference effector memory CD8-positive, alpha-beta T cell",
      "Query effector memory CD8-positive, alpha-beta T cell")],
  colour_by = "group_label",
  point_size = 0.5,
  text_size = 3
) +
  labs(
    title = "UMAP of Combined Reference and Query PBMCs",
    colour = "Dataset: Cell Type"
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3), ncol = 1))

umap_plot

ggsave("figures/covid/umap_visualization.png")

# _____________________________________
# Heatmap IFN Gene Signature (Yoshida)
# _____________________________________

# Define gene lists
yoshida_ifn_signature <- c(
  "BST2", "CMPK2", "EIF2AK2", "EPSTI1", "HERC5", "IFI35", "IFI44L", 
  "IFI6", "IFIT3", "ISG15", "LY6E", "MX1", "MX2", "OAS1", "OAS2", 
  "PARP9", "PLSCR1", "SAMD9", "SAMD9L", "SP110", "STAT1", "TRIM22", 
  "UBE2L6", "XAF1", "IRF7"
)

ifitm_genes <- c("IFITM1", "IFITM2", "IFITM3")
other_genes <- c("NKG7")

# Combine gene lists
all_genes <- c(ifitm_genes, other_genes, yoshida_ifn_signature)

# Function to create heatmap for a given SCE object
create_ifn_heatmap <- function(sce_obj, dataset_name) {
  
  # Check which genes are available
  genes_in_data <- all_genes[all_genes %in% rownames(sce_obj)]
  missing_genes <- setdiff(all_genes, genes_in_data)
  
  if (length(missing_genes) > 0) {
    cat("Missing genes in", dataset_name, ":", paste(missing_genes, collapse = ", "), "\n")
  }
  
  # Extract expression matrix (using logcounts)
  expr_matrix <- as.matrix(assay(sce_obj, "logcounts")[genes_in_data, ])
  
  # Filter out genes with zero variance to avoid clustering issues
  gene_vars <- apply(expr_matrix, 1, var)
  genes_to_keep <- gene_vars > 0
  expr_matrix <- expr_matrix[genes_to_keep, ]
  
  if (sum(!genes_to_keep) > 0) {
    cat("Removed", sum(!genes_to_keep), "genes with zero variance in", dataset_name, "\n")
  }
  
  # Sample cells if too many (for visualization purposes)
  max_cells <- 2000
  if (ncol(expr_matrix) > max_cells) {
    sampled_cells <- sample(ncol(expr_matrix), max_cells)
    expr_matrix <- expr_matrix[, sampled_cells]
    cell_annotations <- colData(sce_obj)[sampled_cells, ]
  } else {
    cell_annotations <- colData(sce_obj)
  }
  
  # Create cell type annotation
  cell_type_col <- if ("cell_type" %in% colnames(cell_annotations)) "cell_type" else "singler_annotations"
  
  # Get unique cell types and create colors
  unique_cell_types <- unique(cell_annotations[[cell_type_col]])
  n_types <- length(unique_cell_types)
  
  # Create color palette for cell types
  type_colors <- rainbow(n_types)
  names(type_colors) <- unique_cell_types
  
  # Create column annotation
  col_ha <- HeatmapAnnotation(
    CellType = cell_annotations[[cell_type_col]],
    col = list(CellType = type_colors),
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  # Create row annotation to distinguish gene groups using nested ifelse
  row_groups <- ifelse(rownames(expr_matrix) %in% ifitm_genes, "IFITM",
                ifelse(rownames(expr_matrix) %in% other_genes, "Other", "IFN_Signature"))
  
  row_ha <- rowAnnotation(
    GeneGroup = row_groups,
    col = list(GeneGroup = c("IFITM" = "red", "Other" = "green", "IFN_Signature" = "blue")),
    width = unit(5, "mm")
  )
  
  # Create color function for expression values
  col_fun <- colorRamp2(
    c(0, quantile(expr_matrix, 0.5), quantile(expr_matrix, 0.95)), 
    c("blue", "white", "red")
  )
  
  # Create heatmap with safer clustering options
  ht <- Heatmap(
    expr_matrix,
    name = "Log Expression",
    col = col_fun,
    top_annotation = col_ha,
    left_annotation = row_ha,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_title = paste("IFITM, NKG7 & IFN Signature Genes -", dataset_name),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete",
    clustering_method_columns = "complete",
    row_split = row_groups,
    column_split = cell_annotations[[cell_type_col]],
    cluster_column_slices = FALSE,
    cluster_row_slices = FALSE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8)
    )
  )
  
  return(ht)
}

# Create heatmaps for both datasets
ht_normal <- create_ifn_heatmap(normal_data, "Reference (Normal)")
ht_covid <- create_ifn_heatmap(covid_data, "Query (COVID)")

# Draw heatmaps
draw(ht_normal, heatmap_legend_side = "right", annotation_legend_side = "right")
draw(ht_covid, heatmap_legend_side = "right", annotation_legend_side = "right")

ggsave("figures/covid/ifn_signature_heatmap.png")

# Put both heatmaps side by side for direct comparison
ht_list <- ht_normal + ht_covid
draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "right")
