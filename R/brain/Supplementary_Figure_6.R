# -------------------------------------------
# ZEISEL BRAIN - Supplementary Figures 6
# -------------------------------------------

# Load libraries
library(scDiagnostics)
library(SingleCellExperiment)
library(scater)
library(SingleR)
library(scRNAseq)
library(ggplot2)
library(dplyr)
library(ggrepel)

# -------------------------------------------

# ____________________________
# 1. SETUP & DATA PREPARATION
# _____________________________

message("Loading and preparing Zeisel Brain Data...")
sce <- scRNAseq::ZeiselBrainData()
sce <- logNormCounts(sce)
sce$true_cell_type <- sce$level1class

set.seed(0)
split_idx <- sample(seq_len(ncol(sce)), size = floor(ncol(sce) * 0.7))
ref_sce <- sce[, split_idx]
query_sce <- sce[, -split_idx]

ref_sce <- ref_sce[, ref_sce$true_cell_type != "pyramidal SS"]

message("Running SingleR annotation...")
pred <- SingleR(test = query_sce, ref = ref_sce, labels = ref_sce$true_cell_type)
query_sce$SingleR_annotation <- pred$labels

target_cluster <- "pyramidal CA1"
query_target_cells <- query_sce[, query_sce$SingleR_annotation == target_cluster]
valid_cells <- query_target_cells$true_cell_type %in% c("pyramidal CA1", "pyramidal SS")
truth_labels <- query_target_cells$true_cell_type[valid_cells]
is_true_anomaly <- truth_labels == "pyramidal SS"

# __________________________
# 2. GRID SEARCH EXECUTION
# __________________________

message("Running Grid Search: Reconstruction Error...")

pc_list_re <- list("1:3" = 1:3, "1:5" = 1:5, "1:10" = 1:10)
hvg_list_re <- c(50, 100, 500)
mad_list_re <- c(2, 3) 

results_re <- data.frame()

for (hvg in hvg_list_re) {
    for (pc_name in names(pc_list_re)) {
        pc <- pc_list_re[[pc_name]]
        if (max(pc) >= hvg) next 
        
        for (mad_val in mad_list_re) {
            
            res <- suppressMessages(suppressWarnings(scDiagnostics::calculateReconstructionError(
                reference_data=ref_sce, query_data=query_sce, 
                ref_cell_type_col="true_cell_type", query_cell_type_col="SingleR_annotation",
                cell_types=target_cluster, pc_subset=pc, n_hvgs=hvg, mad_multiplier=mad_val 
            )))
            
            preds <- res[[target_cluster]]$query_anomaly[valid_cells]
            
            TP <- sum(is_true_anomaly == TRUE & preds == TRUE, na.rm=TRUE)
            TN <- sum(is_true_anomaly == FALSE & preds == FALSE, na.rm=TRUE)
            FP <- sum(is_true_anomaly == FALSE & preds == TRUE, na.rm=TRUE)
            FN <- sum(is_true_anomaly == TRUE & preds == FALSE, na.rm=TRUE)
            
            results_re <- rbind(results_re, data.frame(
                HVGs = as.character(hvg), 
                PCs = pc_name, 
                MAD_Threshold = paste0("MAD (", mad_val, "x)"),
                Sensitivity = TP/(TP+FN), Specificity = TN/(TN+FP)
            ))
        }
    }
}

# Create the explicit "HVG (PC)" label
results_re$ConfigLabel <- paste0(results_re$HVGs, " (", results_re$PCs, ")")
results_re$HVGs <- factor(results_re$HVGs, levels = as.character(hvg_list_re))

# __________________
# 3. VISUALIZATION
# __________________

min_x <- min(0.4, min(results_re$Specificity, na.rm = TRUE) - 0.05)
min_y <- min(0.4, min(results_re$Sensitivity, na.rm = TRUE) - 0.05)

pub_theme <- theme_bw(base_size = 22) + 
    theme(
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold", size = 22), 
        axis.text = element_text(size = 18),                 
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 20), 
        legend.text = element_text(size = 18),                 
        legend.key.size = unit(1.2, "cm"),                     
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
    )

p_re_sens <- ggplot(results_re, aes(x = Specificity, y = Sensitivity, color = HVGs, shape = MAD_Threshold)) +
    annotate("rect", xmin = 0.8, xmax = Inf, ymin = 0.8, ymax = Inf, alpha = 0.12, fill = "forestgreen") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50", linewidth = 1) +
    geom_vline(xintercept = 0.8, linetype = "dashed", color = "gray50", linewidth = 1) +
    geom_point(size = 7, alpha = 0.85) + 
    ggrepel::geom_text_repel(aes(label = ConfigLabel), 
                             size = 6.5, 
                             fontface = "bold", 
                             show.legend = FALSE, box.padding = 0.6, max.overlaps = Inf) +
    scale_color_brewer(palette = "Dark2") +
    coord_cartesian(xlim = c(min_x, 1.02), ylim = c(min_y, 1.02)) + 
    labs(x = "Specificity (True Negative Rate)", y = "Sensitivity (True Positive Rate)", 
         color = "Local HVGs", shape = "Threshold") +
    pub_theme

print(p_re_sens)

ggsave("figures/supp/brain/Fig_S6_reconstruction_error_tuning.png", p_re_sens, width = 16, height = 11, dpi = 600)

# ________
# Summary
# ________

print("Supplementary Figure S6 complete!")
print("Saved in: figures/supp/brain/")