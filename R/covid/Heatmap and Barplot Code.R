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
     query_data = covid_data,
     reference_data = normal_data,
     query_cell_type_col = "azimuth_celltype_l1",
     ref_cell_type_col = "author_cell_type",
     cell_types = c("CD14_mono"), 
     pc_subset = 1:3,
     n_top_loadings = 50,
     genes_to_analyze = yoshida_ifn_signature,
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
     pc_subset = 1:3,
     plot_type = c("heatmap", "barplot", "boxplot")[1],
     plot_by = c("p_adjusted", "top_loading")[1],
     n_genes = 50,
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
     pc_subset = 1:3,
     plot_type = c("heatmap", "barplot", "boxplot")[2],
     plot_by = c("p_adjusted", "top_loading")[1],
     n_genes = 50,
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

# Key genes for fold change calculation
key_ifn_genes <- yoshida_ifn_signature

# Run anomaly detection on FULL DATA (all genes, genome-wide)
# to get the same anomalies as Panel B
anomaly_output_full <- detectAnomaly(
    reference_data = normal_data,
    query_data = covid_data,
    ref_cell_type_col = "author_cell_type",
    query_cell_type_col = "azimuth_celltype_l1",
    cell_types = "CD14_mono",
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.45,
    assay_name = "logcounts",
    max_cells_ref = NULL,
    max_cells_query = NULL
)

# Extract anomalous cell indices (from genome-wide detection)
ref_anomaly_logical <- anomaly_output_full[["CD14_mono"]][["reference_anomaly"]]
query_anomaly_logical <- anomaly_output_full[["CD14_mono"]][["query_anomaly"]]

# Subset to CD14 monocytes
normal_cd14 <- normal_data[, normal_data$author_cell_type == "CD14_mono"]
covid_cd14 <- covid_data[, covid_data$azimuth_celltype_l1 == "CD14_mono"]

# Apply anomaly filtering (non-anomalous reference vs anomalous query, matching anomaly_comparison=TRUE)
normal_cd14_non_anom <- normal_cd14[, !ref_anomaly_logical]
covid_cd14_anom <- covid_cd14[, query_anomaly_logical]
covid_cd14_non_anom <- covid_cd14[, !query_anomaly_logical]
covid_cd14_all <- covid_cd14  # All query cells

# Extract logcounts for the three key genes
normal_expr <- logcounts(normal_cd14_non_anom)[key_ifn_genes, ]
covid_anom_expr <- logcounts(covid_cd14_anom)[key_ifn_genes, ]
covid_non_anom_expr <- logcounts(covid_cd14_non_anom)[key_ifn_genes, ]
covid_all_expr <- logcounts(covid_cd14_all)[key_ifn_genes, ]

# Calculate pseudo-bulk means
mean_normal <- rowMeans(normal_expr)
mean_covid_anom <- rowMeans(covid_anom_expr)
mean_covid_non_anom <- rowMeans(covid_non_anom_expr)
mean_covid_all <- rowMeans(covid_all_expr)

# Calculate log2-fold changes (pseudo-bulk)
log2fc_anom <- mean_covid_anom - mean_normal
log2fc_non_anom <- mean_covid_non_anom - mean_normal
log2fc_all <- mean_covid_all - mean_normal

# Print results
cat("===== log2-fold change (pseudo-bulk, all-gene anomaly detection) =====\n")
cat("\nAnomaly vs Reference (non-anom):\n")
for (i in seq_along(key_ifn_genes)) {
  cat(sprintf("  %s: %.2f\n", key_ifn_genes[i], log2fc_anom[key_ifn_genes[i]]))
}

cat("\nNon-anomalous vs Reference (non-anom):\n")
for (i in seq_along(key_ifn_genes)) {
  cat(sprintf("  %s: %.2f\n", key_ifn_genes[i], log2fc_non_anom[key_ifn_genes[i]]))
}

cat("\nAll CD14+ monocytes vs Reference (non-anom):\n")
for (i in seq_along(key_ifn_genes)) {
  cat(sprintf("  %s: %.2f\n", key_ifn_genes[i], log2fc_all[key_ifn_genes[i]]))
}