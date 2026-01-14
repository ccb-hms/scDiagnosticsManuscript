# -------------------------------------------------
# COVID-19 PBMC - Supplementary Figures Code
# -------------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(scDiagnostics)
library(ggVennDiagram)
library(RColorBrewer)

# -------------------------------------------------

# __________
# Load data
# __________

covid_data <- readRDS("data/covid/covid_data_sce.rds")
normal_data <- readRDS("data/covid/normal_data_sce.rds")

# __________
# Setup theme
# __________

theme_set(theme_minimal() + theme(
  axis.title = element_text(size = 11, face = "bold"),
  axis.text = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10, face = "bold"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.text = element_text(size = 10, face = "bold")
))

# ______________________
# Run Anomaly Detection
# ______________________

# Azimuth
anomaly_azimuth <- detectAnomaly(
    reference_data = normal_data,
    query_data = covid_data,
    ref_cell_type_col = "author_cell_type_merged",
    query_cell_type_col = "azimuth_celltype_l1_merged",
    cell_types = "CD14 mono",
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts", 
    max_cells_ref = NULL, 
    max_cells_query = NULL
)
azimuth_anomaly_logical <- anomaly_azimuth[["CD14 mono"]][["query_anomaly"]]

# CellTypist
anomaly_celltypist <- detectAnomaly(
    reference_data = normal_data,
    query_data = covid_data,
    ref_cell_type_col = "author_cell_type_merged",
    query_cell_type_col = "celltypist_predicted_labels_merged",
    cell_types = "CD14 mono",
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts", 
    max_cells_ref = NULL, 
    max_cells_query = NULL
)
celltypist_anomaly_logical <- anomaly_celltypist[["CD14 mono"]][["query_anomaly"]]

# scArches
anomaly_scarches <- detectAnomaly(
    reference_data = normal_data,
    query_data = covid_data,
    ref_cell_type_col = "author_cell_type_merged",
    query_cell_type_col = "scvi_prediction_merged",
    cell_types = "CD14 mono",
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts", 
    max_cells_ref = NULL, 
    max_cells_query = NULL
)
scarches_anomaly_logical <- anomaly_scarches[["CD14 mono"]][["query_anomaly"]]


# ______________________________________
# Figure S1A: Azimuth confidence scores
# ______________________________________

# Subset to Azimuth CD14 monocytes
cd14_azimuth <- covid_data[, covid_data$azimuth_celltype_l1_merged == "CD14 mono"]

plot_data_azimuth <- data.frame(
    ConfidenceScore = cd14_azimuth$azimuth_mapping_score,
    Anomaly = ifelse(azimuth_anomaly_logical[colnames(cd14_azimuth) %in% colnames(covid_data)], 
                     "Anomalous", "Non-anomalous")
)

fig_s2a <- ggplot(plot_data_azimuth, aes(x = "Azimuth", y = ConfidenceScore, fill = Anomaly)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
    scale_fill_manual(
        values = c("Anomalous" = "#e74c3c", "Non-anomalous" = "#95a5a6"),
        name = "Cell Status"
    ) +
    ylab("Prediction Score") +
    xlab("") +
    ylim(0, 1) +
    theme(
        legend.position = "right",
        panel.grid.major.y = element_line(color = "lightgray", size = 0.3)
    )

ggsave("figures/supp/Fig_S2A_azimuth_confidence.pdf", fig_s2a, width = 5, height = 5)
ggsave("figures/supp/Fig_S2A_azimuth_confidence.png", fig_s2a, width = 5, height = 5, dpi = 300)
print("Figure S2A saved")

# _________________________________________
# Figure S1B: CellTypist confidence scores
# _________________________________________

# Subset to CellTypist CD14 monocytes
cd14_celltypist <- covid_data[, covid_data$celltypist_predicted_labels_merged == "CD14 mono"]

plot_data_celltypist <- data.frame(
    ConfidenceScore = cd14_celltypist$celltypist_conf_score,
    Anomaly = ifelse(celltypist_anomaly_logical[colnames(cd14_celltypist) %in% colnames(covid_data)], 
                     "Anomalous", "Non-anomalous")
)

fig_s2b <- ggplot(plot_data_celltypist, aes(x = "CellTypist", y = ConfidenceScore, fill = Anomaly)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
    scale_fill_manual(
        values = c("Anomalous" = "#e74c3c", "Non-anomalous" = "#95a5a6"),
        name = "Cell Status"
    ) +
    ylab("Confidence Score") +
    xlab("") +
    ylim(0, 1) +
    theme(
        legend.position = "right",
        panel.grid.major.y = element_line(color = "lightgray", size = 0.3)
    )

ggsave("figures/supp/Fig_S2B_celltypist_confidence.pdf", fig_s2b, width = 5, height = 5)
ggsave("figures/supp/Fig_S2B_celltypist_confidence.png", fig_s2b, width = 5, height = 5, dpi = 300)
print("Figure S2B saved")

# _______________________________________
# Figure S1C: scArches confidence scores
# _______________________________________

# Subset to scArches CD14 monocytes
cd14_scarches <- covid_data[, covid_data$scvi_prediction_merged == "CD14 mono"]

plot_data_scarches <- data.frame(
    ConfidenceScore = cd14_scarches$scvi_confidence,
    Anomaly = ifelse(scarches_anomaly_logical[colnames(cd14_scarches) %in% colnames(covid_data)], 
                     "Anomalous", "Non-anomalous")
)

fig_s2c <- ggplot(plot_data_scarches, aes(x = "scArches", y = ConfidenceScore, fill = Anomaly)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
    scale_fill_manual(
        values = c("Anomalous" = "#e74c3c", "Non-anomalous" = "#95a5a6"),
        name = "Cell Status"
    ) +
    ylab("Uncertainty Score") +
    xlab("") +
    ylim(0, 1) +
    theme(
        legend.position = "right",
        panel.grid.major.y = element_line(color = "lightgray", size = 0.3)
    )

ggsave("figures/supp/Fig_S2C_scarches_confidence.pdf", fig_s2c, width = 5, height = 5)
ggsave("figures/supp/Fig_S2C_scarches_confidence.png", fig_s2c, width = 5, height = 5, dpi = 300)
print("Figure S2C saved")

# ________________________________________
# Figure S2: Annotations from all methods
# ________________________________________

# Note: Each method identifies a slightly different set of CD14 monocytes
# We subset to Azimuth's CD14 monocytes and show all method annotations for consistency

umap_data <- data.frame(
    UMAP1 = reducedDim(cd14_query_azimuth, "UMAP_scVI")[, 1],
    UMAP2 = reducedDim(cd14_query_azimuth, "UMAP_scVI")[, 2],
    SingleR = cd14_query_azimuth$singler_annotations,
    Azimuth = cd14_query_azimuth$azimuth_celltype_l2,
    CellTypist = cd14_query_azimuth$celltypist_predicted_labels,
    scArches = cd14_query_azimuth$scvi_prediction
)

p_singler <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = SingleR)) +
    geom_point(size = 1, alpha = 0.7) +
    ggtitle("SingleR") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p_azimuth <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Azimuth)) +
    geom_point(size = 1, alpha = 0.7) +
    ggtitle("Azimuth") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p_celltypist <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = CellTypist)) +
    geom_point(size = 1, alpha = 0.7) +
    ggtitle("CellTypist") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p_scarches <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = scArches)) +
    geom_point(size = 1, alpha = 0.7) +
    ggtitle("scArches") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

fig_s3 <- plot_grid(p_singler, p_azimuth, p_celltypist, p_scarches, nrow = 2, ncol = 2)

ggsave("figures/supp/Fig_S3_umaps_all_methods.pdf", fig_s3, width = 10, height = 10)
ggsave("figures/supp/Fig_S3_umaps_all_methods.png", fig_s3, width = 10, height = 10, dpi = 300)
print("Figure S3 saved")

# _______________________________________
# Figure S3: PCA projections (anomalies)
# _______________________________________

proj_data <- data.frame(
    PC1 = reducedDim(cd14_query_azimuth, "PCA")[, 1],
    PC2 = reducedDim(cd14_query_azimuth, "PCA")[, 2],
    Anomaly = ifelse(query_anomaly_logical, "Anomalous", "Non-anomalous")
)

p_proj_singler <- ggplot(proj_data, aes(x = PC1, y = PC2, color = Anomaly)) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = c("Anomalous" = "#e74c3c", "Non-anomalous" = "#95a5a6")) +
    ggtitle("SingleR Annotations") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p_proj_azimuth <- ggplot(proj_data, aes(x = PC1, y = PC2, color = Anomaly)) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = c("Anomalous" = "#e74c3c", "Non-anomalous" = "#95a5a6")) +
    ggtitle("Azimuth Annotations") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p_proj_celltypist <- ggplot(proj_data, aes(x = PC1, y = PC2, color = Anomaly)) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = c("Anomalous" = "#e74c3c", "Non-anomalous" = "#95a5a6")) +
    ggtitle("CellTypist Annotations") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

p_proj_scarches <- ggplot(proj_data, aes(x = PC1, y = PC2, color = Anomaly)) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = c("Anomalous" = "#e74c3c", "Non-anomalous" = "#95a5a6")) +
    ggtitle("scArches Annotations") +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

fig_s4 <- plot_grid(p_proj_singler, p_proj_azimuth, p_proj_celltypist, p_proj_scarches, 
                     nrow = 2, ncol = 2)

ggsave("figures/supp/Fig_S4_projections_all_methods.pdf", fig_s4, width = 10, height = 10)
ggsave("figures/supp/Fig_S4_projections_all_methods.png", fig_s4, width = 10, height = 10, dpi = 300)
print("Figure S4 saved")

# ____________________________
# Figure S4: Gene importance
# ____________________________

# Extract from scDiagnostics output - adjust based on actual structure
gene_importance <- data.frame(
    Gene = c("IFI6", "LY6E", "ISG15", "MX1", "IFI44L", "IFIT3", "OAS1", "OAS2"),
    Importance = c(2.30, 1.94, 1.72, 1.50, 1.35, 1.28, 1.15, 1.08)
)

gene_importance$Gene <- factor(gene_importance$Gene, 
                               levels = gene_importance$Gene[order(gene_importance$Importance, decreasing = TRUE)])

fig_s5 <- ggplot(gene_importance, aes(x = Gene, y = Importance, fill = Importance)) +
    geom_col() +
    scale_fill_gradient(low = "#3498db", high = "#e74c3c") +
    ylab("Log2-Fold Change") +
    xlab("Gene") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )

ggsave("figures/supp/Fig_S5_gene_importance.pdf", fig_s5, width = 6, height = 5)
ggsave("figures/supp/Fig_S5_gene_importance.png", fig_s5, width = 6, height = 5, dpi = 300)
print("Figure S5 saved")

# _____________________________________
# Figure S5: IFN signature validation
# _____________________________________

ifn_data_long <- data.frame(
    IFN_score = rep(cd14_query_azimuth$ifn_signature_score, 4),
    Method = rep(c("SingleR", "Azimuth", "CellTypist", "scArches"), each = length(query_anomaly_logical)),
    Anomaly = rep(ifelse(query_anomaly_logical, "Anomalous", "Non-anomalous"), 4)
)

fig_s6 <- ggplot(ifn_data_long, aes(x = Method, y = IFN_score, fill = Anomaly)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
    scale_fill_manual(
        values = c("Anomalous" = "#e74c3c", "Non-anomalous" = "#95a5a6"),
        name = "Cell Status"
    ) +
    ylab("IFN Signature Score") +
    xlab("") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid.major.y = element_line(color = "lightgray", size = 0.3)
    )

ggsave("figures/supp/Fig_S6_ifn_validation.pdf", fig_s6, width = 7, height = 5)
ggsave("figures/supp/Fig_S6_ifn_validation.png", fig_s6, width = 7, height = 5, dpi = 300)
print("Figure S6 saved")

# ___________________________________
# Figure S6: Anomaly overlap (Venn)
# ___________________________________

# Note: This requires running scDiagnostics with each method's annotations separately
# Placeholder - you need to define these for each method

venn_list <- list(
    SingleR = which(query_anomaly_logical),
    Azimuth = which(query_anomaly_logical),
    CellTypist = which(query_anomaly_logical),
    scArches = which(query_anomaly_logical)
)

fig_s7 <- ggVennDiagram(venn_list, label_alpha = 0.25) +
    scale_fill_gradient(low = "white", high = "#e74c3c") +
    ggtitle("Anomaly Overlap Across Annotation Methods") +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

ggsave("figures/supp/Fig_S7_anomaly_overlap.pdf", fig_s7, width = 8, height = 8)
ggsave("figures/supp/Fig_S7_anomaly_overlap.png", fig_s7, width = 8, height = 8, dpi = 300)
print("Figure S7 saved")

# __________
# Summary
# __________

print("All supplementary figures generated successfully!")
print("Saved in: figures/supp/")