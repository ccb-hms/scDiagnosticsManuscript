# -----------------------------------------------
# MERFISH Mouse Colon IBD - Preliminary Analyses
# -----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scDiagnostics)
library(ggplot2)
 
# -----------------------------------------------

# Load the processed SCE objects
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

# Cast the annotation columns as characters
healthy_data$tier2 <- as.character(healthy_data$tier2)
dss9_data$predicted_labels <- as.character(dss9_data$azimuth_celltype_l1)

# _________________________
# PCA Projections Boxplots
# _________________________

# Check distribution of cell types 
cbind(table(healthy_data$tier2),
      table(dss9_data$azimuth_celltype_l1))

# Check IAE and IAF predicted classes by Azimuth
table(dss9_data$azimuth_celltype_l1[(grepl("IAE", dss9_data$tier2))])
table(dss9_data$azimuth_celltype_l1[(grepl("IAE 1", dss9_data$tier2))])
table(dss9_data$azimuth_celltype_l1[(grepl("IAE 2", dss9_data$tier2))])
table(dss9_data$azimuth_celltype_l1[(grepl("IAE 3", dss9_data$tier2))])

table(dss9_data$azimuth_celltype_l1[(grepl("IAF", dss9_data$tier2))])
table(dss9_data$azimuth_celltype_l1[(grepl("IAF 1", dss9_data$tier2))])
table(dss9_data$azimuth_celltype_l1[(grepl("IAF 2", dss9_data$tier2))])
table(dss9_data$azimuth_celltype_l1[(grepl("IAF 3", dss9_data$tier2))])
table(dss9_data$azimuth_celltype_l1[(grepl("IAF 4", dss9_data$tier2))])
table(dss9_data$azimuth_celltype_l1[(grepl("Fibro 4", dss9_data$tier2))])

# PCA boxplot
boxplot_pca <- boxplotPCA(query_data = dss9_data,
                          reference_data = healthy_data, 
                          query_cell_type_col = "azimuth_celltype_l1", 
                          ref_cell_type_col = "tier2", 
                          cell_types = c("Fibro 2", "Fibro 6", "Fibro 7", "Fibro 13",
                                         "SMC 1", "Pericyte 1",
                                         "Colonocytes", "Stem cells", "TA"), 
                          pc_subset = 1:5,
                          max_cells_ref = NULL,
                          max_cells_query = NULL)
boxplot_pca

ggsave("figures/merfish/boxplot_pca.png")

# _________________
# PCA Projections
# _________________

# PCA scatterplot
scatter_pca <- plotCellTypePCA(query_data = dss9_data,
                               reference_data = healthy_data, 
                               query_cell_type_col = "azimuth_celltype_l1", 
                               ref_cell_type_col = "tier2", 
                               cell_types = c("Colonocytes", "Stem cells", "TA"), 
                               pc_subset = 1:4, 
                               diagonal_facet = "ridge",
                               max_cells_ref = 1000,
                               max_cells_query = 1000)
scatter_pca

ggsave("figures/merfish/cell_type_pca.png")

# __________________
# Anomaly Detection
# __________________

# Anomaly Detection
anomaly_output <- detectAnomaly(query_data = dss9_data,
                                reference_data = healthy_data, 
                                query_cell_type_col = "azimuth_celltype_l1", 
                                ref_cell_type_col = "tier2", 
                                cell_types = c("Fibro 2", "Fibro 6", "Fibro 7", "Fibro 13",
                                              "SMC 1", "Pericyte 1",
                                              "Colonocytes"), 
                                pc_subset = 1:5,
                                n_tree = 500,
                                anomaly_threshold = 0.55,
                                max_cells_ref = NULL,
                                max_cells_query = NULL)
plot(anomaly_output,
     cell_type = "Colonocytes",
     pc_subset = 1:5,
     data_type = "query",
     n_tree = 500,
     diagonal_facet = "ridge",
     upper_facet = "blank", 
     max_cells_query = 2000)

ggsave("figures/merfish/anomaly_plot.png")

# _______________________________
# Top Genes Distributional Shift
# _______________________________

# Compute distributional shifts for genes with top loadings
gene_shifts <- calculateTopLoadingGeneShifts(
     query_data = dss9_data,
     reference_data = healthy_data, 
     query_cell_type_col = "azimuth_celltype_l1", 
     ref_cell_type_col = "tier2", 
     cell_types = c("Fibro 2", "Fibro 6", "Fibro 7", "Fibro 13",
                    "SMC 1", "Pericyte 1",
                    "Colonocytes", "Stem cells", "TA"), 
     pc_subset = 1:5,
     n_top_loadings = 50,
     assay_name = "logcounts",
     p_value_threshold = 0.05,
     adjust_method = "fdr", 
     detect_anomalies = TRUE, 
     anomaly_comparison = TRUE,
     anomaly_threshold = 0.5,
     max_cells_query = NULL,
     max_cells_ref = NULL)

plot(gene_shifts,
     cell_type = "Fibro 6",
     pc_subset = 1:5,
     plot_type = c("heatmap", "boxplot")[1],
     plot_by = c("p_adjusted", "top_loading")[1],
     n_genes = 5,
     significance_threshold = 0.05, 
     show_anomalies = TRUE,
     pseudo_bulk = TRUE,
     max_cells_query = NULL,
     max_cells_ref = NULL)

ggsave("figures/merfish/gene_shifts_heatmap.png")

# __________________
# Marker Expression
# __________________

# Plot the marker expression
plotMarkerExpression(reference_data = healthy_data,
                     query_data = dss9_data,
                     query_cell_type_col = "azimuth_celltype_l1",
                     ref_cell_type_col = "tier2",
                     gene_name = "Lgr4",
                     cell_type = "Colonocytes",
                     normalization = "min_max")

# Plot the marker expression
plotMarkerExpression(reference_data = healthy_data,
                     query_data = dss9_data,
                     query_cell_type_col = "azimuth_celltype_l1",
                     ref_cell_type_col = "tier2",
                     gene_name = "Vil1",
                     cell_type = "Colonocytes",
                     normalization = "min_max")

ggsave("figures/merfish/marker_expression.png")

# ______________
# PC Regression
# ______________

# Regress the PC data (with reference data)
regress_res2 <- regressPC(query_data = dss9_data,
                          reference_data = healthy_data, 
                          query_cell_type_col = "azimuth_celltype_l1", 
                          ref_cell_type_col = "tier2", 
                          cell_types = c("Fibro 2", "Fibro 6", "Fibro 7", "Fibro 13",
                                         "SMC 1", "Pericyte 1",
                                         "Colonocytes"), 
                          pc_subset = 1:10,
                          max_cells_query = NULL,
                          max_cells_ref = NULL)
plot(regress_res2, plot_type = "r_squared")
plot(regress_res2, plot_type = "variance_contribution")
plot(regress_res2, plot_type = "coefficient_heatmap")

