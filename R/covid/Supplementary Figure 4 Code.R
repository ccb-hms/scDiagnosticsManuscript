# -----------------------------------------------
# COVID-19 PBMC - Supplementary Figure 4 Code
# -----------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scDiagnostics)
library(SingleR)

# -----------------------------------------------

# __________
# Load data
# __________

covid_data <- readRDS("data/covid/covid_data_sce.rds")
normal_data <- readRDS("data/covid/normal_data_sce.rds")

set.seed(0)

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

# SingleR
anomaly_singler <- detectAnomaly(
    reference_data = normal_data,
    query_data = covid_data,
    ref_cell_type_col = "author_cell_type_merged",
    query_cell_type_col = "singler_annotations_merged",
    cell_types = "CD14 mono",
    pc_subset = 1:3,
    n_tree = 500,
    anomaly_threshold = 0.5,
    assay_name = "logcounts", 
    max_cells_ref = NULL, 
    max_cells_query = NULL
)
singler_anomaly_logical <- anomaly_singler[["CD14 mono"]][["query_anomaly"]]

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
# Figure S4A: SingleR Delta Distribution
# ______________________________________

cd14_singler <- covid_data[, covid_data$singler_annotations_merged == "CD14 mono"]

scores_cd14 <- covid_data$singler_scores
rownames(scores_cd14) <- colnames(covid_data)
scores_cd14 <- scores_cd14[colnames(cd14_singler), ]

deltas <- apply(scores_cd14, 1, function(x) max(x) - median(x))

delta_data <- data.frame(
    Delta = deltas,
    Anomaly = factor(ifelse(singler_anomaly_logical[colnames(cd14_singler) %in% colnames(covid_data)], 
                     "Anomalous", "Non-anomalous"),
                     levels = c("Non-anomalous", "Anomalous"))
)

fig_s4a <- ggplot(delta_data, aes(x = Delta, fill = Anomaly)) +
    geom_density(alpha = 0.8, linewidth = 1, color = "white") +
    facet_wrap(~Anomaly, ncol = 1, labeller = labeller(Anomaly = c(
        "Anomalous" = "Anomalous CD14+ Monocytes",
        "Non-anomalous" = "Non-Anomalous CD14+ Monocytes"
    ))) +
    scale_fill_manual(
        values = c("Non-anomalous" = "#5B7C99", "Anomalous" = "#C1666B"),
        name = "Cell Status"
    ) +
    xlab("Delta Score") +
    ylab("Density") +
    ggtitle("SingleR") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 11, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 10, color = "#4B5563", family = "sans"),
        strip.text = element_text(size = 11, face = "bold", color = "#1F2937", family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.35),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.6, fill = NA),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA)
    )

# ______________________________________
# Figure S4B: Azimuth Confidence Scores
# ______________________________________

cd14_azimuth <- covid_data[, covid_data$azimuth_celltype_l1_merged == "CD14 mono"]

plot_data_azimuth <- data.frame(
    ConfidenceScore = cd14_azimuth$azimuth_mapping_score,
    Anomaly = factor(ifelse(azimuth_anomaly_logical[colnames(cd14_azimuth) %in% colnames(covid_data)], 
                     "Anomalous", "Non-anomalous"), 
                     levels = c("Non-anomalous", "Anomalous"))
)

fig_s4b <- ggplot(plot_data_azimuth, aes(x = ConfidenceScore, fill = Anomaly)) +
    geom_density(alpha = 0.8, linewidth = 1, color = "white") +
    facet_wrap(~Anomaly, ncol = 1, labeller = labeller(Anomaly = c(
        "Anomalous" = "Anomalous CD14+ Monocytes",
        "Non-anomalous" = "Non-Anomalous CD14+ Monocytes"
    ))) +
    scale_fill_manual(
        values = c("Non-anomalous" = "#5B7C99", "Anomalous" = "#C1666B"),
        name = "Cell Status"
    ) +
    xlab("Prediction Score") +
    ylab("Density") +
    xlim(0, 1) +
    ggtitle("Azimuth") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 11, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 10, color = "#4B5563", family = "sans"),
        strip.text = element_text(size = 11, face = "bold", color = "#1F2937", family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.35),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.6, fill = NA),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA)
    )

# _________________________________________
# Figure S4C: CellTypist Confidence Scores
# _________________________________________

cd14_celltypist <- covid_data[, covid_data$celltypist_predicted_labels_merged == "CD14 mono"]

plot_data_celltypist <- data.frame(
    ConfidenceScore = cd14_celltypist$celltypist_conf_score,
    Anomaly = factor(ifelse(celltypist_anomaly_logical[colnames(cd14_celltypist) %in% colnames(covid_data)], 
                     "Anomalous", "Non-anomalous"),
                     levels = c("Non-anomalous", "Anomalous"))
)

fig_s4c <- ggplot(plot_data_celltypist, aes(x = ConfidenceScore, fill = Anomaly)) +
    geom_density(alpha = 0.8, linewidth = 1, color = "white") +
    facet_wrap(~Anomaly, ncol = 1, labeller = labeller(Anomaly = c(
        "Anomalous" = "Anomalous CD14+ Monocytes",
        "Non-anomalous" = "Non-Anomalous CD14+ Monocytes"
    ))) +
    scale_fill_manual(
        values = c("Non-anomalous" = "#5B7C99", "Anomalous" = "#C1666B"),
        name = "Cell Status"
    ) +
    xlab("Confidence Score") +
    ylab("Density") +
    ggtitle("CellTypist") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 11, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 10, color = "#4B5563", family = "sans"),
        strip.text = element_text(size = 11, face = "bold", color = "#1F2937", family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.35),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.6, fill = NA),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA)
    )

# _______________________________________
# Figure S4D: scArches Confidence Scores
# _______________________________________

cd14_scarches <- covid_data[, covid_data$scvi_prediction_merged == "CD14 mono"]

plot_data_scarches <- data.frame(
    ConfidenceScore = cd14_scarches$scvi_confidence,
    Anomaly = factor(ifelse(scarches_anomaly_logical[colnames(cd14_scarches) %in% colnames(covid_data)], 
                     "Anomalous", "Non-anomalous"),
                     levels = c("Non-anomalous", "Anomalous"))
)

fig_s4d <- ggplot(plot_data_scarches, aes(x = ConfidenceScore, fill = Anomaly)) +
    geom_density(alpha = 0.8, linewidth = 1, color = "white") +
    facet_wrap(~Anomaly, ncol = 1, labeller = labeller(Anomaly = c(
        "Anomalous" = "Anomalous CD14+ Monocytes",
        "Non-anomalous" = "Non-Anomalous CD14+ Monocytes"
    ))) +
    scale_fill_manual(
        values = c("Non-anomalous" = "#5B7C99", "Anomalous" = "#C1666B"),
        name = "Cell Status"
    ) +
    xlab("Uncertainty Score") +
    ylab("Density") +
    xlim(0, 1) +
    ggtitle("scArches") +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "#1F2937", family = "sans"),
        axis.title = element_text(size = 11, face = "bold", color = "#374151", family = "sans"),
        axis.text = element_text(size = 10, color = "#4B5563", family = "sans"),
        strip.text = element_text(size = 11, face = "bold", color = "#1F2937", family = "sans"),
        legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.35),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#D1D5DB", linewidth = 0.6, fill = NA),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.background = element_rect(fill = "#FAFBFC", color = NA)
    )

# _________________
# Combine 2x2 grid
# _________________

fig_s4_combined <- plot_grid(fig_s4a, fig_s4b,
                              fig_s4c, fig_s4d,
                              nrow = 2, ncol = 2, labels = "AUTO",
                              label_size = 12, label_fontface = "bold")

ggsave("figures/supp/Fig_S4_confidence_scores.png", fig_s4_combined, width = 14, height = 12, dpi = 600)

# __________
# Summary
# __________

print("Supplementary Figure S4 complete!")
print("Saved in: figures/supp/")

# ___________________________________________________
# Table S4: Wilcoxon Tests - Anomaly vs Confidence
# ___________________________________________________

# Function to create table row
create_method_row <- function(covid_data, method_name, annotation_col, 
                              score_col, anomaly_logical, annotation_cell_type) {
    
    # Get CD14+ cells by name
    cd14_mask <- covid_data[[annotation_col]] == annotation_cell_type
    cd14_cell_names <- colnames(covid_data)[cd14_mask]  # ← FIXED: was cd14_indices
    
    # Extract scores by name
    scores <- covid_data[[score_col]]
    
    # Handle matrix case (SingleR) 
    if (is.matrix(scores)) {
        rownames(scores) <- colnames(covid_data)
        scores <- rowMeans(scores[cd14_cell_names, c("CD14_mono", "CD83_CD14_mono")])  
    } else {
        names(scores) <- colnames(covid_data)
        scores <- scores[cd14_cell_names]  
    }
    
    # Extract anomaly status by matching to cell names
    names(anomaly_logical) <- gsub("Query_", "", names(anomaly_logical))
    anomaly_status <- anomaly_logical[cd14_cell_names]
    
    # Split by anomaly status
    anom_scores <- scores[anomaly_status == TRUE]
    non_anom_scores <- scores[anomaly_status == FALSE]
    
    # Remove NAs
    anom_scores <- anom_scores[!is.na(anom_scores)]
    non_anom_scores <- non_anom_scores[!is.na(non_anom_scores)]
    
    # Check if we have enough data
    if (length(anom_scores) < 3 || length(non_anom_scores) < 3) {
        warning(paste(method_name, ": Not enough observations for Wilcoxon test"))
        return(NULL)
    }
    
    # Wilcoxon test
    wilcox_result <- wilcox.test(anom_scores, non_anom_scores)
    
    # Format p-value
    if (wilcox_result$p.value < 0.001) {
        p_val_str <- "<0.001"
    } else {
        p_val_str <- paste0("p=", format(round(wilcox_result$p.value, 4), nsmall = 4))
    }
    
    # Create row
    row <- data.frame(
        Method = method_name,
        N_Anomalous = length(anom_scores),
        N_Non_Anomalous = length(non_anom_scores),
        Median_Score_Anom = round(median(anom_scores, na.rm = TRUE), 2),
        Median_Score_Non_Anom = round(median(non_anom_scores, na.rm = TRUE), 2),
        Wilcoxon_p_value = p_val_str,
        stringsAsFactors = FALSE
    )
    
    return(row)
}

# Generate table data
table_s4_list <- list(
    # SingleR (uses CD14_mono with underscore)
    create_method_row(covid_data, "SingleR", "singler_annotations_merged", "singler_scores", 
                     singler_anomaly_logical, "CD14 mono"),
    
    # Azimuth
    create_method_row(covid_data, "Azimuth", "azimuth_celltype_l1_merged", "azimuth_mapping_score", 
                     azimuth_anomaly_logical, "CD14 mono"),
    
    # CellTypist
    create_method_row(covid_data, "CellTypist", "celltypist_predicted_labels_merged", "celltypist_conf_score", 
                     celltypist_anomaly_logical, "CD14 mono"),
    
    # scArches
    create_method_row(covid_data, "scArches", "scvi_prediction_merged", "scvi_confidence", 
                     scarches_anomaly_logical, "CD14 mono")
)

# Remove NULL entries
table_s4_data <- do.call(rbind, Filter(Negate(is.null), table_s4_list))
rownames(table_s4_data) <- NULL
print(table_s4_data)
