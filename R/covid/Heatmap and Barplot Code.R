# ----------------------------------
# COVID-19 PBMC - Heatmap Code
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

# __________
# Load Data
# __________

# Load the processed SCE objects
normal_data <- readRDS("data/covid/normal_data_sce.rds")
covid_data <- readRDS("data/covid/covid_data_sce.rds")

# Setting the seed
set.seed(0)

# _______________________________
# Top Genes Distributional Shift
# _______________________________

#  Define gene lists
yoshida_ifn_signature <- c(
  "BST2", "CMPK2", "EIF2AK2", "EPSTI1", "HERC5", "IFI35", "IFI44L", 
  "IFI6", "IFIT3", "ISG15", "LY6E", "MX1", "MX2", "OAS1", "OAS2", 
  "PARP9", "PLSCR1", "SAMD9", "SAMD9L", "SP110", "STAT1", "TRIM22", 
  "UBE2L6", "XAF1", "IRF7"
)

# Compute distributional shifts for genes with top loadings
gene_shifts <- calculateGeneShifts(
     query_data = covid_data[yoshida_ifn_signature,],
     reference_data = normal_data[yoshida_ifn_signature,],
     query_cell_type_col = "azimuth_celltype_l1",
     ref_cell_type_col = "author_cell_type",
     cell_types = c("CD14_mono"), 
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

# Heatmap of IFN genes
heatmap_object <- plot(gene_shifts,
     cell_type = "CD14_mono",
     pc_subset = 1:5,
     plot_type = c("heatmap", "barplot", "boxplot")[1],
     plot_by = c("p_adjusted", "top_loading")[1],
     n_genes = 25,
     significance_threshold = 0.05,
     show_anomalies = TRUE,
     pseudo_bulk = TRUE,
     cluster_cols = TRUE,
     max_cells_ref = 5000,
     max_cells_query = 5000,
    draw_plot = FALSE)

# Save the heatmap plot
png("figures/covid/covid_ifn_heatmap.png")
draw(heatmap_object)
dev.off()

# Barplot of IFN genes
plot(gene_shifts,
     cell_type = "CD14_mono",
     pc_subset = 1:5,
     plot_type = c("heatmap", "barplot", "boxplot")[2],
     plot_by = c("p_adjusted", "top_loading")[1],
     n_genes = 25,
     significance_threshold = 0.05,
     show_anomalies = TRUE,
     pseudo_bulk = TRUE,
     cluster_cols = TRUE,
     max_cells_ref = 5000,
     max_cells_query = 5000)

# Save the barplot
ggsave("figures/covid/gene_shifts_barplot.png")

# __________________________________
# Calculate Pseudo-Bulk Fold Change
# __________________________________

# Key genes
key_ifn_genes <- c("IFI6", "ISG15", "LY6E")

# Run anomaly detection on CD14 monocytes
anomaly_output <- detectAnomaly(
    reference_data = normal_data,
    query_data = covid_data,
    ref_cell_type_col = "author_cell_type",
    query_cell_type_col = "azimuth_celltype_l1",
    cell_types = "CD14_mono",
    pc_subset = 1:5,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts",
    max_cells_ref = 5000,
    max_cells_query = 5000
)

# Extract anomalous cell indices
anomaly_logical <- anomaly_output[["CD14_mono"]][["query_anomaly"]]

# Subset to CD14 monocytes
normal_cd14 <- normal_data[, normal_data$author_cell_type == "CD14_mono"]
covid_cd14 <- covid_data[, covid_data$azimuth_celltype_l1 == "CD14_mono"]

# Extract logcounts and calculate pseudo-bulk means
normal_expr <- logcounts(normal_cd14)[key_ifn_genes, ]
covid_anom_expr <- logcounts(covid_cd14)[key_ifn_genes, anomaly_logical]

mean_normal <- rowMeans(normal_expr)
mean_covid_anom <- rowMeans(covid_anom_expr)
log2fc <- mean_covid_anom - mean_normal

# Print results
cat("===== log2-fold change (pseudo-bulk) =====\n")
for (i in seq_along(key_ifn_genes)) {
  cat(sprintf("%s: %.2f\n", key_ifn_genes[i], log2fc[i]))
}
