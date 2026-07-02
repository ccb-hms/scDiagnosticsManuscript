# ------------------------------------------
# ZEISEL BRAIN - Supplementary Figures 5
# ------------------------------------------

# Load libraries
library(scDiagnostics)
library(SingleCellExperiment)
library(scater)
library(SingleR)
library(scRNAseq)
library(ggplot2)
library(dplyr)
library(ggrepel)

# ------------------------------------------

# _____________________________
# 1. SETUP & DATA PREPARATION
# _____________________________

message("Loading and preparing Zeisel Brain Data...")
sce <- scRNAseq::ZeiselBrainData()
sce <- logNormCounts(sce)
sce$true_cell_type <- sce$level1class

set.seed(0)
split_idx <- sample(seq_len(ncol(sce)), size = floor(ncol(sce) * 0.7))
ref_sce <- scDiagnostics::processPCA(sce[, split_idx], n_hvgs = 1000)
query_sce <- sce[, -split_idx]

# Force a misannotation by removing 'pyramidal SS' from reference
ref_sce <- ref_sce[, ref_sce$true_cell_type != "pyramidal SS"]

message("Running SingleR annotation...")
pred <- SingleR(test = query_sce, ref = ref_sce, labels = ref_sce$true_cell_type)
query_sce$SingleR_annotation <- pred$labels

# Ground Truth Setup
target_cluster <- "pyramidal CA1"
query_target_cells <- query_sce[, query_sce$SingleR_annotation == target_cluster]
valid_cells <- query_target_cells$true_cell_type %in% c("pyramidal CA1", "pyramidal SS")
truth_labels <- query_target_cells$true_cell_type[valid_cells]
is_true_anomaly <- truth_labels == "pyramidal SS"

# __________________________
# 2. GRID SEARCH EXECUTION
# __________________________

message("Running Grid Search: Isolation Forest...")

pc_list <- list("1:3"=1:3, "1:5"=1:5, "1:10"=1:10)
hvg_list <- c(50, 100, 500)
thresholds <- list(
    list(name="Abs (0.5)", method="absolute", mad_val=2, abs_val=0.5),
    list(name="MAD (2x)", method="MAD", mad_val=2, abs_val=0.5),
    list(name="MAD (3x)", method="MAD", mad_val=3, abs_val=0.5)
)

results_if <- data.frame()

calc_metrics <- function(config_name, mode, param, thresh_name, preds) {
    TP <- sum(is_true_anomaly == TRUE & preds == TRUE, na.rm=TRUE)
    TN <- sum(is_true_anomaly == FALSE & preds == FALSE, na.rm=TRUE)
    FP <- sum(is_true_anomaly == FALSE & preds == TRUE, na.rm=TRUE)
    FN <- sum(is_true_anomaly == TRUE & preds == FALSE, na.rm=TRUE)
    data.frame(Config=config_name, Mode=mode, Parameter=param, Threshold=thresh_name, 
               Sensitivity = TP/(TP+FN), Specificity = TN/(TN+FP))
}

for (pc_name in names(pc_list)) {
    for (thresh in thresholds) {
        res <- suppressMessages(suppressWarnings(scDiagnostics::detectAnomaly(
            reference_data=ref_sce, query_data=query_sce, ref_cell_type_col="true_cell_type", query_cell_type_col="SingleR_annotation",
            cell_types=target_cluster, pc_subset=pc_list[[pc_name]], n_hvgs=1000, n_tree=300,
            threshold_method=thresh$method, mad_multiplier=thresh$mad_val, anomaly_threshold=thresh$abs_val
        )))
        preds <- res[[target_cluster]]$query_anomaly[valid_cells]
        results_if <- rbind(results_if, calc_metrics(paste0("PCA_", pc_name), "Global PCA", pc_name, thresh$name, preds))
    }
}

for (hvg in hvg_list) {
    for (thresh in thresholds) {
        res <- suppressMessages(suppressWarnings(scDiagnostics::detectAnomaly(
            reference_data=ref_sce, query_data=query_sce, ref_cell_type_col="true_cell_type", query_cell_type_col="SingleR_annotation",
            cell_types=target_cluster, pc_subset=NULL, n_hvgs=hvg, n_tree=300,
            threshold_method=thresh$method, mad_multiplier=thresh$mad_val, anomaly_threshold=thresh$abs_val
        )))
        preds <- res[[target_cluster]]$query_anomaly[valid_cells]
        results_if <- rbind(results_if, calc_metrics(paste0("HVG_", hvg), "Targeted HVGs", as.character(hvg), thresh$name, preds))
    }
}

# _________________
# 3. VISUALIZATION
# _________________

# Dynamically calculate limits so no point is left behind
min_x <- min(0.4, min(results_if$Specificity, na.rm = TRUE) - 0.05)
min_y <- min(0.4, min(results_if$Sensitivity, na.rm = TRUE) - 0.05)

pub_theme <- theme_bw(base_size = 13) +
    theme(
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 11),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
    )

p_if_sens <- ggplot(results_if, aes(x = Specificity, y = Sensitivity, color = Mode, shape = Threshold)) +
    annotate("rect", xmin = 0.8, xmax = Inf, ymin = 0.8, ymax = Inf, alpha = 0.12, fill = "forestgreen") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0.8, linetype = "dashed", color = "gray50") +
    geom_point(size = 4, alpha = 0.85) +
    ggrepel::geom_text_repel(aes(label = Parameter), size = 3.5, fontface = "bold", 
                             show.legend = FALSE, box.padding = 0.4, max.overlaps = Inf) +
    scale_color_manual(values = c("Global PCA" = "#C44E52", "Targeted HVGs" = "#4C72B0")) +
    coord_cartesian(xlim = c(min_x, 1.02), ylim = c(min_y, 1.02)) + 
    labs(x = "Specificity (True Negative Rate)", y = "Sensitivity (True Positive Rate)", 
         color = "Feature Space", shape = "Threshold") +
    pub_theme

print(p_if_sens)

ggsave("figures/supp/brain/Fig_S5_isolation_forest_tuning.png", p_if_sens, width = 15, height = 10, dpi = 600)

# ________
# Summary
# ________

print("Supplementary Figure S5 complete!")
print("Saved in: figures/supp/brain/")