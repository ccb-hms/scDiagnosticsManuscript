# ---------------------------------------------
# ZEISEL BRAIN - Supplementary Figures 1-4
# ---------------------------------------------

library(scDiagnostics)
library(SingleCellExperiment)
library(scater)
library(SingleR)
library(scRNAseq)
library(pROC)     
library(PRROC)    
library(ggplot2)
library(dplyr)
library(patchwork)

# ---------------------------------------------

# ____________________________
# 1. SETUP, DATA PREP, & UMAP
# ____________________________

message("Loading and preparing Zeisel Brain Data...")
sce <- scRNAseq::ZeiselBrainData()
sce <- logNormCounts(sce)

# Use broad level1class for the entire analysis
sce$true_cell_type <- sce$level1class

set.seed(0)
split_idx <- sample(seq_len(ncol(sce)), size = floor(ncol(sce) * 0.7))
base_ref <- sce[, split_idx]
base_query <- sce[, -split_idx]

# --- Generate Professional UMAP for Panel A ---
message("Generating UMAP for visualization...")
sce_umap <- scater::runPCA(sce, ncomponents = 30) 
sce_umap <- scater::runUMAP(sce_umap, dimred = "PCA", n_neighbors = 50, min_dist = 0.5) 

umap_df <- data.frame(
    UMAP1 = reducedDim(sce_umap, "UMAP")[,1],
    UMAP2 = reducedDim(sce_umap, "UMAP")[,2],
    CellType = sce_umap$true_cell_type
)

# Highlight all test populations with clearer labels
umap_df$Highlight <- "Other Cell Types"
umap_df$Highlight[umap_df$CellType == "pyramidal CA1"] <- "Pyramidal CA1 (Base Neuronal Subtype)"
umap_df$Highlight[umap_df$CellType == "pyramidal SS"] <- "Pyramidal SS (Hidden Related Subtype)"
umap_df$Highlight[umap_df$CellType == "microglia"] <- "Microglia (Hidden Rare Population)"
umap_df$Highlight[umap_df$CellType == "astrocytes_ependymal"] <- "Astrocytes (Hidden Distinct Lineage)"

plot_order <- c("Other Cell Types", "Pyramidal CA1 (Base Neuronal Subtype)", 
                "Pyramidal SS (Hidden Related Subtype)", "Microglia (Hidden Rare Population)", 
                "Astrocytes (Hidden Distinct Lineage)")
umap_df$Highlight <- factor(umap_df$Highlight, levels = plot_order)
umap_df <- umap_df[order(umap_df$Highlight), ]

umap_colors <- c(
    "Other Cell Types" = "#E8E8E8",            
    "Pyramidal CA1 (Base Neuronal Subtype)" = "#4C72B0",        # Blue
    "Pyramidal SS (Hidden Related Subtype)" = "#1B9E77",        # Green
    "Microglia (Hidden Rare Population)" = "#E6AB02",           # Gold
    "Astrocytes (Hidden Distinct Lineage)" = "#D95F02"          # Burnt Orange
)

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Highlight)) +
    geom_point(size = 1.0, stroke = 0, alpha = 0.85) + 
    scale_color_manual(values = umap_colors) +
    theme_bw(base_size = 20) + 
    labs(title = "A. Ground Truth Manifold", 
         x = "UMAP1", y = "UMAP2") + 
    theme(
        panel.grid = element_blank(),     
        axis.text = element_blank(),      
        axis.ticks = element_blank(),   
        
        # FIX: Explicitly pull axis titles closer to the plot panel using negative margins
        axis.title.x = element_text(face = "bold", size = 18, margin = margin(t = -20)), 
        axis.title.y = element_text(face = "bold", size = 18, margin = margin(r = -20)),
        
        plot.title = element_text(face = "bold", size = 22, hjust = 0.5, margin = margin(b = 5)),
        plot.subtitle = element_text(size = 18, hjust = 0.5, margin = margin(b = 15)),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 16), 
        legend.margin = margin(t = -15, b = 0),
        legend.box.margin = margin(t = -20, b = 0, l = 0, r = 0) 
    ) +
    guides(color = guide_legend(override.aes = list(size = 6, alpha = 1), nrow = 3))

# _________________________
# 2. EVALUATION FUNCTIONS
# _________________________

compute_metrics <- function(scores, truth_logical) {
    if(length(scores) == 0 || length(truth_logical) == 0) return(c(AUROC=NA, AUPRC=NA))
    valid_idx <- !is.na(scores) & !is.na(truth_logical)
    scores <- scores[valid_idx]; truth_logical <- truth_logical[valid_idx]
    if(length(unique(truth_logical)) < 2) return(c(AUROC=NA, AUPRC=NA))
    roc_obj <- suppressMessages(roc(response = truth_logical, predictor = scores, quiet = TRUE))
    fg <- scores[truth_logical == TRUE]; bg <- scores[truth_logical == FALSE]
    pr_obj <- suppressWarnings(pr.curve(scores.class0 = fg, scores.class1 = bg, curve = FALSE))
    return(c(AUROC = as.numeric(auc(roc_obj)), AUPRC = pr_obj$auc.integral))
}

calc_sens_spec <- function(preds, truth_logical) {
    TP <- sum(truth_logical == TRUE & preds == TRUE, na.rm=TRUE)
    TN <- sum(truth_logical == FALSE & preds == FALSE, na.rm=TRUE)
    FP <- sum(truth_logical == FALSE & preds == TRUE, na.rm=TRUE)
    FN <- sum(truth_logical == TRUE & preds == FALSE, na.rm=TRUE)
    return(c(Sensitivity = TP/(TP+FN), Specificity = TN/(TN+FP)))
}

run_test <- function(ref_sce, query_sce, target_hidden_cluster) {
    pred <- SingleR(test = query_sce, ref = ref_sce, labels = ref_sce$true_cell_type)
    query_sce$SingleR_annotation <- pred$labels
    hidden_cells_predictions <- query_sce$SingleR_annotation[query_sce$true_cell_type == target_hidden_cluster]
    target_ref_cluster <- names(which.max(table(hidden_cells_predictions)))
    
    res_if <- suppressMessages(suppressWarnings(scDiagnostics::detectAnomaly(
        reference_data=ref_sce, query_data=query_sce, ref_cell_type_col="true_cell_type", query_cell_type_col="SingleR_annotation",
        cell_types=target_ref_cluster, pc_subset=NULL, n_hvgs=100, n_tree=500
    )))
    res_re <- suppressMessages(suppressWarnings(scDiagnostics::calculateReconstructionError(
        reference_data=ref_sce, query_data=query_sce, ref_cell_type_col="true_cell_type", query_cell_type_col="SingleR_annotation",
        cell_types=target_ref_cluster, pc_subset=1:5, n_hvgs=100
    )))
    
    query_target_cells <- query_sce[, query_sce$SingleR_annotation == target_ref_cluster]
    valid_cells <- query_target_cells$true_cell_type %in% c(target_ref_cluster, target_hidden_cluster)
    truth_logical <- query_target_cells$true_cell_type[valid_cells] == target_hidden_cluster
    
    if_scores <- res_if[[target_ref_cluster]]$query_anomaly_scores[valid_cells] 
    re_scores <- res_re[[target_ref_cluster]]$query_reconstruction_errors[valid_cells] 
    
    if_preds <- res_if[[target_ref_cluster]]$query_anomaly[valid_cells] 
    re_preds <- res_re[[target_ref_cluster]]$query_anomaly[valid_cells]
    
    m_if <- compute_metrics(if_scores, truth_logical)
    m_re <- compute_metrics(re_scores, truth_logical)
    
    ss_if <- calc_sens_spec(if_preds, truth_logical)
    ss_re <- calc_sens_spec(re_preds, truth_logical)
    
    return(data.frame(
        Method = c("Isolation_Forest", "Reconstruction_Error"), 
        AUROC = c(m_if["AUROC"], m_re["AUROC"]),
        AUPRC = c(m_if["AUPRC"], m_re["AUPRC"]),
        Sensitivity = c(ss_if["Sensitivity"], ss_re["Sensitivity"]),
        Specificity = c(ss_if["Specificity"], ss_re["Specificity"])
    ))
}

all_results <- list()

# __________________
# 3. RUN SCENARIOS
# __________________

message("\nRunning Baselines...")

# Distinct (Astrocytes)
res_A_dist <- run_test(base_ref[, base_ref$true_cell_type != "astrocytes_ependymal"], base_query, "astrocytes_ependymal")
res_A_dist$Test <- "Distinct\n(Astrocytes)"; res_A_dist$TestGroup <- "Baseline"; all_results[[length(all_results)+1]] <- res_A_dist

# Related (Pyr SS)
res_A_rel <- run_test(base_ref[, base_ref$true_cell_type != "pyramidal SS"], base_query, "pyramidal SS")
res_A_rel$Test <- "Related\n(Pyr SS)"; res_A_rel$TestGroup <- "Baseline"; all_results[[length(all_results)+1]] <- res_A_rel

# Rare (Microglia)
res_A_rare <- run_test(base_ref[, base_ref$true_cell_type != "microglia"], base_query, "microglia")
res_A_rare$Test <- "Rare\n(Microglia)"; res_A_rare$TestGroup <- "Baseline"; all_results[[length(all_results)+1]] <- res_A_rare

# Helper to run gradients dynamically
run_gradient_suite <- function(hidden_cluster, group_name, imb_levels = c(300, 100, 50, 20, 10)) {
    message(paste("Running Gradients for:", group_name))
    ref_grad <- base_ref[, base_ref$true_cell_type != hidden_cluster]
    
    pred <- SingleR(test = base_query, ref = ref_grad, labels = ref_grad$true_cell_type)
    hidden_preds <- pred$labels[base_query$true_cell_type == hidden_cluster]
    target_ref <- names(which.max(table(hidden_preds)))
    
    # Noise
    for(nl in c(0.0, 0.1, 0.2, 0.3, 0.4)) {
        ref_tmp <- ref_grad
        n_shuffle <- floor(ncol(ref_tmp) * nl)
        if(n_shuffle > 0) {
            shuffle_idx <- sample(seq_len(ncol(ref_tmp)), size = n_shuffle)
            ref_tmp$true_cell_type[shuffle_idx] <- sample(ref_tmp$true_cell_type[shuffle_idx])
        }
        res_tmp <- run_test(ref_tmp, base_query, hidden_cluster)
        res_tmp$Test <- "Noise"; res_tmp$X_Value <- as.character(nl); res_tmp$TestGroup <- group_name
        all_results[[length(all_results)+1]] <<- res_tmp
    }
    
    # Imbalance
    target_ref_idx <- which(ref_grad$true_cell_type == target_ref)
    actual_size <- length(target_ref_idx)
    max_clean_size <- floor(actual_size / 10) * 10
    valid_levels <- unique(sapply(imb_levels, function(x) min(x, max_clean_size)))
    valid_levels <- sort(valid_levels, decreasing = TRUE)
    
    for(cl in valid_levels) {
        ref_tmp <- ref_grad
        ca1_idx <- which(ref_tmp$true_cell_type == target_ref)
        
        keep_ca1 <- sample(ca1_idx, size = cl)
        ref_tmp <- ref_tmp[, -setdiff(ca1_idx, keep_ca1)]
        
        res_tmp <- run_test(ref_tmp, base_query, hidden_cluster)
        res_tmp$Test <- "Imbalance"
        res_tmp$X_Value <- as.character(cl)
        res_tmp$TestGroup <- group_name
        all_results[[length(all_results)+1]] <<- res_tmp
    }
    
    # Batch Effect
    for(bl in c(0.0, 0.5, 1.0, 1.5, 2.0)) {
        query_tmp <- base_query
        if(bl > 0) {
            genes_to_shift <- sample(seq_len(nrow(query_tmp)), size = floor(nrow(query_tmp) * 0.2))
            temp_counts <- as.matrix(logcounts(query_tmp))
            temp_counts[genes_to_shift, ] <- temp_counts[genes_to_shift, ] + bl
            logcounts(query_tmp) <- temp_counts
        }
        res_tmp <- run_test(ref_grad, query_tmp, hidden_cluster)
        res_tmp$Test <- "Batch"; res_tmp$X_Value <- as.character(bl); res_tmp$TestGroup <- group_name
        all_results[[length(all_results)+1]] <<- res_tmp
    }
}

run_gradient_suite("pyramidal SS", "Related", imb_levels = c(300, 100, 50, 20, 10))
run_gradient_suite("microglia", "Rare", imb_levels = c(300, 100, 50, 20, 10))

plot_data <- bind_rows(all_results)

# __________________________
# 4. PLOTTING PANELS B - H
# __________________________

my_theme <- theme_bw(base_size = 20) + theme(
    panel.grid.minor = element_blank(), 
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5, margin = margin(b=10)),
    axis.title = element_text(face = "bold", size = 18), 
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.title = element_blank(), 
    legend.text = element_text(size = 18),
    legend.margin = margin(t = 0, b = 0)
)

# Plot Generator Function 
generate_benchmark_figure <- function(df_plot, metric_col, metric_label, y_lims, pal_colors, pt_shape, hline_y = 0.5) {
    
    method_labels <- c("Isolation_Forest" = "Isolation Forest", "Reconstruction_Error" = "Reconstruction Error")
    
    # --- ROW 1: Baseline ---
    df_base <- df_plot |> filter(TestGroup == "Baseline")
    df_base$Test <- factor(df_base$Test, levels = c("Distinct\n(Astrocytes)", "Related\n(Pyr SS)", "Rare\n(Microglia)"))
    
    p_base <- ggplot(df_base, aes(x = Test, y = .data[[metric_col]], fill = Method)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
        scale_fill_manual(values = pal_colors, labels = method_labels) +
        coord_cartesian(ylim = y_lims) + 
        labs(title = paste0("B. Baseline Detection Accuracy (", metric_label, ")"), 
             x = "Missing Cell State", y = metric_label) + 
        my_theme + theme(legend.position = "bottom") 
        
    if(!is.na(hline_y)) { p_base <- p_base + geom_hline(yintercept = hline_y, linetype = "dashed", color = "gray50", linewidth = 1) }

    # --- ROW 2 & 3 GENERATOR ---
    create_line_plot <- function(df, test_type, title, x_label, is_factor=FALSE) {
        sub_df <- df |> filter(Test == test_type)
        if(is_factor) {
            numeric_vals <- sort(as.numeric(unique(sub_df$X_Value)), decreasing = TRUE)
            sub_df$X_Value <- factor(sub_df$X_Value, levels = as.character(numeric_vals))
        } else {
            sub_df$X_Value <- as.numeric(sub_df$X_Value)
            if(test_type == "Noise") sub_df$X_Value <- sub_df$X_Value * 100
        }
        
        p_line <- ggplot(sub_df, aes(x = X_Value, y = .data[[metric_col]], color = Method, group = Method)) +
            geom_line(linewidth = 1.5) + 
            geom_point(size = 4.0, shape = pt_shape, fill = "white", stroke = 1.5) + 
            scale_color_manual(values = pal_colors, labels = method_labels) +
            coord_cartesian(ylim = y_lims) + 
            labs(title = title, x = x_label, y = metric_label) + 
            my_theme
            
        if(!is.na(hline_y)) { p_line <- p_line + geom_hline(yintercept = hline_y, linetype = "dashed", color = "gray50", linewidth = 1) }
        return(p_line)
    }

    # Row 2 (Related)
    df_rel <- df_plot |> filter(TestGroup == "Related")
    pC <- create_line_plot(df_rel, "Noise", "C. Label Noise (Related)", "% of Labels Shuffled")
    pD <- create_line_plot(df_rel, "Imbalance", "D. Imbalance (Related)", "N Cells in Ref Cluster", TRUE)
    pE <- create_line_plot(df_rel, "Batch", "E. Batch Effect (Related)", "Mean Shift (+LogCounts)")

    # Row 3 (Rare)
    df_rare <- df_plot |> filter(TestGroup == "Rare")
    pF <- create_line_plot(df_rare, "Noise", "F. Label Noise (Rare)", "% of Labels Shuffled")
    pG <- create_line_plot(df_rare, "Imbalance", "G. Imbalance (Rare)", "N Cells in Ref Cluster", TRUE)
    pH <- create_line_plot(df_rare, "Batch", "H. Batch Effect (Rare)", "Mean Shift (+LogCounts)")

    # _______________________
    # ASSEMBLE NESTED PLOT
    # _______________________
    row1 <- (p_umap + theme(legend.position = "bottom")) | p_base
    row1 <- row1 + plot_layout(widths = c(1.3, 1))
    
    row_lines <- (pC | pD | pE) / (pF | pG | pH) + 
                 plot_layout(guides = "collect") & 
                 theme(legend.position = "bottom")

    final_plot <- row1 / row_lines + plot_layout(heights = c(1.2, 2))
    return(final_plot)
}

# --- Define metric-specific palettes & shapes ---
# SWAPPED COLORS: Isolation_Forest is now RED, Reconstruction_Error is now BLUE

# 1. AUROC (Standard Colors, Circles, Baseline 0.5)
pal_auroc <- c("Isolation_Forest" = "#C44E52", "Reconstruction_Error" = "#4C72B0")
p_auroc <- generate_benchmark_figure(plot_data, "AUROC", "AUROC", c(0.4, 1.0), pal_auroc, 21, hline_y = 0.5)

# 2. AUPRC (Distinct Colors, Diamonds, NO Baseline line)
pal_auprc <- c("Isolation_Forest" = "#A83C40", "Reconstruction_Error" = "#325C99")
p_auprc <- generate_benchmark_figure(plot_data, "AUPRC", "AUPRC", c(0.0, 1.0), pal_auprc, 23, hline_y = NA)

# 3. Sensitivity (Darker Colors, Triangles, Baseline 0.5)
pal_sens <- c("Isolation_Forest" = "#9E2A2B", "Reconstruction_Error" = "#1D4E89")
p_sens <- generate_benchmark_figure(plot_data, "Sensitivity", "Sensitivity", c(0.0, 1.0), pal_sens, 24, hline_y = 0.5)

# 4. Specificity (Lighter/Muted Colors, Squares, Baseline 0.5)
pal_spec <- c("Isolation_Forest" = "#D47A7B", "Reconstruction_Error" = "#7FA1C3")
p_spec <- generate_benchmark_figure(plot_data, "Specificity", "Specificity", c(0.0, 1.0), pal_spec, 22, hline_y = 0.5)

# Save figures 
ggsave("figures/supp/brain/Fig_S1_anomaly_detection_AUROC.png", p_auroc, width = 18, height = 18, dpi = 600)
ggsave("figures/supp/brain/Fig_S2_anomaly_detection_AUPRC.png", p_auprc, width = 18, height = 18, dpi = 600)
ggsave("figures/supp/brain/Fig_S3_anomaly_detection_Sensitivity.png", p_sens, width = 18, height = 18, dpi = 600)
ggsave("figures/supp/brain/Fig_S4_anomaly_detection_Specificity.png", p_spec, width = 18, height = 18, dpi = 600)

print("Supplementary Figures S1-4 complete!")
print("Saved in: figures/supp/brain/")

# __________________________
# 5. TEXT METRICS EXTRACTOR 
# __________________________

message("\n=== METRICS SUMMARY FOR MANUSCRIPT TEXT ===")
summary_table <- plot_data |>
    select(TestGroup, Test, X_Value, Method, AUROC, AUPRC, Sensitivity, Specificity) |>
    mutate(across(where(is.numeric), ~round(., 3))) 

message("\n--- Baseline Results ---")
print(summary_table |> filter(TestGroup == "Baseline") |> arrange(Test, Method))

message("\n--- Gradient Results (Related) ---")
print(summary_table |> filter(TestGroup == "Related") |> arrange(Test, as.numeric(X_Value), Method))

message("\n--- Gradient Results (Rare) ---")
print(summary_table |> filter(TestGroup == "Rare") |> arrange(Test, as.numeric(X_Value), Method))