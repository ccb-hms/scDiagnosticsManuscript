# ----------------------------------------
# COVID-19 PBMC - Heatmap & Barplot Code
# ----------------------------------------

# Load libraries
library(SingleCellExperiment)
library(scDiagnostics)
library(ggplot2)
library(scran)
library(scater)
library(ComplexHeatmap)
library(circlize)

# ----------------------------------------

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

# Key genes - SAME as used in calculateGeneShifts
key_ifn_genes <- c("IFI6", "ISG15", "LY6E")
yoshida_ifn_signature <- c(
  "BST2", "CMPK2", "EIF2AK2", "EPSTI1", "HERC5", "IFI35", "IFI44L", 
  "IFI6", "IFIT3", "ISG15", "LY6E", "MX1", "MX2", "OAS1", "OAS2", 
  "PARP9", "PLSCR1", "SAMD9", "SAMD9L", "SP110", "STAT1", "TRIM22", 
  "UBE2L6", "XAF1", "IRF7"
)

# Run anomaly detection (same as in calculateGeneShifts)
anomaly_output <- detectAnomaly(
    reference_data = normal_data[yoshida_ifn_signature,],
    query_data = covid_data[yoshida_ifn_signature,],
    ref_cell_type_col = "author_cell_type",
    query_cell_type_col = "azimuth_celltype_l1",
    cell_types = "CD14_mono",
    pc_subset = 1:5,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts"
)

# Extract anomalous cell indices
ref_anomaly_logical <- anomaly_output[["CD14_mono"]][["reference_anomaly"]]
query_anomaly_logical <- anomaly_output[["CD14_mono"]][["query_anomaly"]]

# Subset to CD14 monocytes
normal_cd14 <- normal_data[yoshida_ifn_signature, normal_data$author_cell_type == "CD14_mono"]
covid_cd14 <- covid_data[yoshida_ifn_signature, covid_data$azimuth_celltype_l1 == "CD14_mono"]

# Apply anomaly filtering
normal_cd14_non_anom <- normal_cd14[, !ref_anomaly_logical]
covid_cd14_anom <- covid_cd14[, query_anomaly_logical]
covid_cd14_non_anom <- covid_cd14[, !query_anomaly_logical]
covid_cd14_all <- covid_cd14  # All query cells

# Calculate pseudo-bulk means using the 25 genes
mean_normal <- rowMeans(logcounts(normal_cd14_non_anom))
mean_covid_anom <- rowMeans(logcounts(covid_cd14_anom))
mean_covid_non_anom <- rowMeans(logcounts(covid_cd14_non_anom))
mean_covid_all <- rowMeans(logcounts(covid_cd14_all))

# Calculate log2-fold changes
log2fc_anom <- mean_covid_anom - mean_normal
log2fc_non_anom <- mean_covid_non_anom - mean_normal
log2fc_all <- mean_covid_all - mean_normal

# Print for your three genes
cat("===== log2-fold change (anomaly_comparison=TRUE) =====\n")
cat("\nAnomaly vs Reference (non-anom):\n")
for (i in seq_along(key_ifn_genes)) {
  if (key_ifn_genes[i] %in% yoshida_ifn_signature) {
    cat(sprintf("  %s: %.2f\n", key_ifn_genes[i], log2fc_anom[key_ifn_genes[i]]))
  }
}

cat("\nNon-anomaly vs Reference (non-anom):\n")
for (i in seq_along(key_ifn_genes)) {
  if (key_ifn_genes[i] %in% yoshida_ifn_signature) {
    cat(sprintf("  %s: %.2f\n", key_ifn_genes[i], log2fc_non_anom[key_ifn_genes[i]]))
  }
}

cat("\nAll CD14+ vs Reference (non-anom):\n")
for (i in seq_along(key_ifn_genes)) {
  if (key_ifn_genes[i] %in% yoshida_ifn_signature) {
    cat(sprintf("  %s: %.2f\n", key_ifn_genes[i], log2fc_all[key_ifn_genes[i]]))
  }
}

# _____________________________
# Pairwise Comparison Analysis 
# _____________________________

# Calculate pairwise log2-fold changes for all 25 genes
# Comparison 1: Anomalous vs Reference (non-anom)
lfc_anom_vs_ref <- log2fc_anom

# Comparison 2: Non-anomalous vs Reference (non-anom)
lfc_nonanom_vs_ref <- log2fc_non_anom

# Comparison 3: Anomalous vs Non-anomalous (within query)
mean_covid_nonanom_direct <- rowMeans(logcounts(covid_cd14_non_anom))
mean_covid_anom_direct <- rowMeans(logcounts(covid_cd14_anom))
lfc_anom_vs_nonanom <- mean_covid_anom_direct - mean_covid_nonanom_direct

# Create data frame for all comparisons
pairwise_comparisons <- data.frame(
    gene = names(lfc_anom_vs_ref),
    anom_vs_ref = lfc_anom_vs_ref,
    nonanom_vs_ref = lfc_nonanom_vs_ref,
    anom_vs_nonanom = lfc_anom_vs_nonanom,
    stringsAsFactors = FALSE
)

cat("===== Pairwise log2-fold change comparisons (all 25 genes) =====\n\n")
print(pairwise_comparisons)

# Statistical testing: Are the three comparisons significantly different?
cat("\n===== Statistical Tests =====\n\n")

# 1. Paired t-test: Anomalous vs Reference vs Non-anomalous vs Reference
cat("1. Paired t-test: Anomalous vs Ref vs Non-anomalous vs Ref\n")
t_test_anom_nonanom <- t.test(lfc_anom_vs_ref, lfc_nonanom_vs_ref, paired = TRUE)
cat(sprintf("   t-statistic: %.4f, p-value: %.4e\n", t_test_anom_nonanom$statistic, t_test_anom_nonanom$p.value))
cat(sprintf("   Mean difference: %.4f\n", mean(lfc_anom_vs_ref - lfc_nonanom_vs_ref)))

# 2. Wilcoxon signed-rank test (non-parametric alternative)
cat("\n2. Wilcoxon signed-rank test: Anomalous vs Ref vs Non-anomalous vs Ref\n")
wilcox_test_anom_nonanom <- wilcox.test(lfc_anom_vs_ref, lfc_nonanom_vs_ref, paired = TRUE)
cat(sprintf("   V-statistic: %.4f, p-value: %.4e\n", wilcox_test_anom_nonanom$statistic, wilcox_test_anom_nonanom$p.value))

# 3. One-way ANOVA on all three comparisons (treating them as repeated measures)
cat("\n3. Summary statistics for each comparison\n")
cat(sprintf("   Anomalous vs Reference:\n      Mean: %.4f, SD: %.4f, Range: [%.4f, %.4f]\n",
            mean(lfc_anom_vs_ref), sd(lfc_anom_vs_ref), 
            min(lfc_anom_vs_ref), max(lfc_anom_vs_ref)))
cat(sprintf("   Non-anomalous vs Reference:\n      Mean: %.4f, SD: %.4f, Range: [%.4f, %.4f]\n",
            mean(lfc_nonanom_vs_ref), sd(lfc_nonanom_vs_ref), 
            min(lfc_nonanom_vs_ref), max(lfc_nonanom_vs_ref)))
cat(sprintf("   Anomalous vs Non-anomalous:\n      Mean: %.4f, SD: %.4f, Range: [%.4f, %.4f]\n",
            mean(lfc_anom_vs_nonanom), sd(lfc_anom_vs_nonanom), 
            min(lfc_anom_vs_nonanom), max(lfc_anom_vs_nonanom)))

# 4. Effect size (Cohen's d for paired samples)
cat("\n4. Effect size (Cohen's d): Anomalous vs Non-anomalous (paired)\n")
mean_diff <- mean(lfc_anom_vs_ref - lfc_nonanom_vs_ref)
sd_diff <- sd(lfc_anom_vs_ref - lfc_nonanom_vs_ref)
cohens_d <- mean_diff / sd_diff
cat(sprintf("   Cohen's d: %.4f\n", cohens_d))
if (abs(cohens_d) < 0.2) {
    effect_interpretation <- "negligible"
} else if (abs(cohens_d) < 0.5) {
    effect_interpretation <- "small"
} else if (abs(cohens_d) < 0.8) {
    effect_interpretation <- "medium"
} else {
    effect_interpretation <- "large"
}
cat(sprintf("   Interpretation: %s effect\n", effect_interpretation))

# 5. Gene-by-gene comparison: identify genes driving the differences
cat("\n5. Gene-by-gene ranking (by anomalous vs non-anomalous difference)\n")
ranked_genes <- pairwise_comparisons[order(pairwise_comparisons$anom_vs_nonanom, decreasing = TRUE), ]
cat("   Top 5 genes with largest difference (Anom vs Non-anom):\n")
print(head(ranked_genes[, c("gene", "anom_vs_nonanom", "anom_vs_ref", "nonanom_vs_ref")], 5))
cat("   Bottom 5 genes with smallest difference (Anom vs Non-anom):\n")
print(tail(ranked_genes[, c("gene", "anom_vs_nonanom", "anom_vs_ref", "nonanom_vs_ref")], 5))

# 6. Summary statistics for ALL fold changes (including all query)
cat("\n6. Summary statistics for ALL fold changes (all 25 genes)\n")
cat(sprintf("   Anomalous vs Reference:\n      Mean: %.4f, SD: %.4f\n",
            mean(lfc_anom_vs_ref), sd(lfc_anom_vs_ref)))
cat(sprintf("   All Query vs Reference:\n      Mean: %.4f, SD: %.4f\n",
            mean(log2fc_all), sd(log2fc_all)))

# 7. Paired t-test comparing anomalous vs all query fold changes
cat("\n7. Paired t-test: Anomalous vs Ref vs All Query vs Ref\n")
t_test_anom_all <- t.test(lfc_anom_vs_ref, log2fc_all, paired = TRUE)
cat(sprintf("   t-statistic: %.4f, p-value: %.4e\n", t_test_anom_all$statistic, t_test_anom_all$p.value))
cat(sprintf("   Mean difference: %.4f\n", mean(lfc_anom_vs_ref - log2fc_all)))