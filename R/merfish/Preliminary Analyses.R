# -----------------------------------------------
# MERFISH Mouse Colon IBD - Preliminary Analyses
# -----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scDiagnostics)
 
# -----------------------------------------------

# Load the processed SCE objects
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

# PCA boxplot
boxplot_pca <- boxplotPCA(query_data = dss9_data,
                          reference_data = healthy_data, 
                          query_cell_type_col = "predicted_labels", 
                          ref_cell_type_col = "tier2", 
                          cell_types = c("IASMC 1", "IAE 2", "IAF 3"), 
                          pc_subset = 1:5, 
                          max_cells = 10000)
boxplot_pca

# PCA scatterplot
scatter_pca <- plotCellTypePCA(query_data  = dss9_data,
                               reference_data = healthy_data, 
                               query_cell_type_col = "predicted_labels", 
                               ref_cell_type_col = "tier2", 
                               cell_types = c("IASMC 1", "IAE 2", "IAF 3"), 
                               pc_subset = 1:5, 
                               diagonal_facet = "ridge",
                               max_cells = 10000)
scatter_pca

# Anomaly Detection
anomaly_output <- detectAnomaly(query_data = dss9_data,
                                reference_data = healthy_data, 
                                query_cell_type_col = "singler_annotations", 
                                ref_cell_type_col = "cell_type", 
                                cell_types = c("CD14-positive monocyte", "myeloid dendritic cell", 
                                               "erythrocyte", "CD4-positive, alpha-beta T cell"), 
                                pc_subset = 1:5,
                                n_tree = 500,
                                anomaly_threshold = 0.6,
                                max_cells = 5000)
plot(anomaly_output,
     cell_type = "CD14-positive monocyte",
     pc_subset = 1:5,
     data_type = "query",
     n_tree = 500,
     diagonal_facet = "ridge",
     upper_facet = "blank")

# Compute distributional shifts for genes with top loadings
gene_shifts <- calculateTopLoadingGeneShifts(reference_data = healthy_data,
                                             query_data = dss9_data,
                                             query_cell_type_col = "singler_annotations",
                                             ref_cell_type_col = "cell_type",
                                             cell_types = c("CD14-positive monocyte", "myeloid dendritic cell", 
                                                            "erythrocyte", "CD4-positive, alpha-beta T cell"), 
                                             pc_subset = 1:5,
                                             n_top_loadings = 50,
                                             assay_name = "logcounts",
                                             p_value_threshold = 0.05,
                                             adjust_method = "fdr")
plot(gene_shifts,
     cell_type = "CD14-positive monocyte",
     pc_subset = 1:4,
     plot_type = c("heatmap", "boxplot")[1],
     plot_by = c("p_adjusted", "top_loading")[1],
     n_genes = 10,
     significance_threshold = 0.05)

# Plot the marker expression
plotMarkerExpression(reference_data = healthy_data,
                     query_data = dss9_data,
                     query_cell_type_col = "singler_annotations",
                     ref_cell_type_col = "cell_type",
                     gene_name = "S100A8",
                     cell_type = "CD14-positive monocyte",
                     normalization = "none")

# Regress the PC data (with reference data)
regress_res2 <- regressPC(query_data = dss9_data,
                          reference_data = healthy_data, 
                          query_cell_type_col = "singler_annotations", 
                          ref_cell_type_col = "cell_type", 
                          query_batch_col = "Site",
                          cell_types = c("CD14-positive monocyte", "myeloid dendritic cell", 
                                         "erythrocyte", "CD4-positive, alpha-beta T cell"), 
                          pc_subset = 1:10, 
                          max_cells = 5000)
plot(regress_res2, plot_type = "r_squared")
plot(regress_res2, plot_type = "variance_contribution")
plot(regress_res2, plot_type = "coefficient_heatmap")

