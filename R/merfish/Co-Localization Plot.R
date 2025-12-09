# --------------------------------------------------------------------
# MERFISH - Spatial Validation: The "Spatial Shift" (Vertical Stack)
# Target Journal: Genome Biology
# --------------------------------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scDiagnostics)
library(FNN) 

# --- 1. Load Data ---
cat("\n--- Loading data ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")
healthy_data$tier2 <- as.character(healthy_data$tier2)
dss9_data$tier2 <- as.character(dss9_data$tier2)
set.seed(0)

# --- 2. Anomaly Detection ---
cat("\n--- Running Anomaly Detection ---\n")
is_fibro_lineage <- function(x) { grepl("^Fibro|^IAF", x) }
healthy_data$analysis_class <- ifelse(is_fibro_lineage(healthy_data$tier2), "Fibroblast_Lineage", "Other")
dss9_data$analysis_class <- ifelse(is_fibro_lineage(dss9_data$tier2), "Fibroblast_Lineage", "Other")

anomaly_output <- detectAnomaly(
    query_data = dss9_data,
    reference_data = healthy_data, 
    query_cell_type_col = "analysis_class", 
    ref_cell_type_col = "analysis_class", 
    cell_types = "Fibroblast_Lineage",
    pc_subset = 1:4, 
    anomaly_threshold = 0.5,
    max_cells_query = NULL,
    max_cells_ref = NULL
)

anomalous_barcodes <- unlist(lapply(anomaly_output, function(res) {
    if(is.null(res)) return(NULL)
    names(res$query_anomaly[res$query_anomaly == TRUE])
}))
anomalous_barcodes <- gsub("Query_", "", anomalous_barcodes)

# --- 3. Spatial Calculations: The 3 Key Niches ---
spatial_df <- data.frame(
    barcode = colnames(dss9_data),
    x = spatialCoords(dss9_data)[, "x"],
    y = spatialCoords(dss9_data)[, "y"],
    ground_truth_type = dss9_data$tier2,
    broad_class = dss9_data$analysis_class
)

# 1. Neutrophils (Inflammatory Target)
pattern_neutro <- "Neutrophil"
coords_neutro <- spatial_df %>% filter(grepl(pattern_neutro, ground_truth_type, ignore.case = TRUE)) %>% select(x, y)

# 2. Repair Epithelium (Damage Target)
pattern_repair <- "IAE|Repair associated|Epithelial \\(Clu\\+\\)"
coords_repair <- spatial_df %>% filter(grepl(pattern_repair, ground_truth_type, ignore.case = TRUE)) %>% select(x, y)

# 3. Stem Cell Niche (Homeostatic Source)
pattern_stem <- "Stem|TA"
coords_stem <- spatial_df %>% filter(grepl(pattern_stem, ground_truth_type, ignore.case = TRUE)) %>% select(x, y)

# Calculate Distances
fibro_df <- spatial_df %>% filter(broad_class == "Fibroblast_Lineage")

get_dists <- function(coords) {
    if(nrow(coords) == 0) return(rep(NA, nrow(fibro_df)))
    FNN::get.knnx(data = as.matrix(coords), query = as.matrix(fibro_df[, c("x", "y")]), k = 1)$nn.dist
}

fibro_df$dist_neutro <- get_dists(coords_neutro)
fibro_df$dist_repair <- get_dists(coords_repair)
fibro_df$dist_stem   <- get_dists(coords_stem)

fibro_df$prediction <- ifelse(fibro_df$barcode %in% anomalous_barcodes, "Anomalous", "Typical")

# --- 4. Prepare Data Structure ---
# A. Raw separated data
df_split <- fibro_df %>%
    select(barcode, prediction, dist_neutro, dist_repair, dist_stem)

# B. Duplicate for "All Fibroblasts" group
df_all <- df_split
df_all$prediction <- "All Fibroblasts"

# C. Combine and Label
combined_df <- rbind(df_split, df_all) %>%
    tidyr::pivot_longer(
        cols = c(dist_neutro, dist_repair, dist_stem), 
        names_to = "target_type", 
        values_to = "distance"
    ) %>%
    mutate(panel_label = case_when(
        target_type == "dist_neutro" ~ "Distance to Neutrophils",
        target_type == "dist_repair" ~ "Distance to Repair Epithelium",
        target_type == "dist_stem"   ~ "Distance to Stem/Crypt Base"
    )) %>%
    # Explicit ordering: Inflammation First, then Stem
    mutate(panel_label = factor(panel_label, levels = c(
        "Distance to Neutrophils", 
        "Distance to Repair Epithelium", 
        "Distance to Stem/Crypt Base"
    )))

# --- 5. Manual Smooth CDF Calculation ---
cat("\n--- Calculating Smooth Curves ---\n")

smooth_curves <- combined_df %>%
    group_by(panel_label, prediction) %>%
    reframe({
        d <- density(distance, from = 0, n = 512)
        cdf <- cumsum(d$y) / sum(d$y)
        data.frame(distance = d$x, cdf = cdf)
    })

smooth_curves$prediction <- factor(smooth_curves$prediction, 
                                   levels = c("Anomalous", "Typical", "All Fibroblasts"))

# --- 6. Smart Limits & Statistics ---

calc_limit <- function(type_name) {
    vals <- combined_df$distance[combined_df$panel_label == type_name]
    ceiling(quantile(vals, 0.90, na.rm=TRUE) / 50) * 50
}

limit_df <- data.frame(
    panel_label = unique(combined_df$panel_label),
    prediction = "Anomalous",
    cdf = 0.5
)
limit_df$distance <- sapply(limit_df$panel_label, calc_limit)

# Calculate Wilcoxon P-values
calc_pval <- function(col_name) {
    p <- wilcox.test(fibro_df[[col_name]] ~ fibro_df$prediction, alternative = "two.sided")$p.value
    format.pval(p, digits=2)
}

col_map <- list(
    "Distance to Neutrophils" = "dist_neutro",
    "Distance to Repair Epithelium" = "dist_repair",
    "Distance to Stem/Crypt Base" = "dist_stem"
)

ann_text <- limit_df %>%
    mutate(
        distance = distance * 0.65,
        cdf = 0.15,
        label = sapply(panel_label, function(x) {
            col <- col_map[[as.character(x)]]
            paste0("Wilcoxon (Diff)\np < ", calc_pval(col))
        })
    )

# --- 7. Generate Plot (VERTICAL STACK) ---
colors_pred <- c(
    "Anomalous"       = "#D9230F",  # Red
    "Typical"         = "#228B22",  # Green
    "All Fibroblasts" = "#606060"   # Dark Gray
)

p <- ggplot(smooth_curves, aes(x = distance, y = cdf, color = prediction)) +
    geom_line(linewidth = 1.2) +
    
    # CHANGED: ncol = 1 (Vertical Stack)
    facet_wrap(~panel_label, scales = "free_x", ncol = 1) +
    
    geom_blank(data = limit_df) +
    
    geom_text(data = ann_text, aes(label = label), 
              color = "black", size = 3.5, fontface = "italic", hjust=0.5) +
    
    scale_color_manual(name = "Population", values = colors_pred) +
    
    labs(
        title = "Co-Localization Plots",
        subtitle = "Anomalous and Typical Fibroblasts",
        x = expression(bold("Distance to Feature ("*mu*"m)")), 
        y = "Cumulative Proportion"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.major.y = element_line(color = "grey92"), 
        legend.position = "top",
        plot.title = element_text(face="bold", hjust=0.5),
        plot.subtitle = element_text(hjust=0.5, size=11)
    ) +
    coord_cartesian(expand = FALSE, ylim = c(0, 1.05))

print(p)

# CHANGED: Height/Width swapped to accommodate vertical stack
ggsave("figures/merfish/co-localization_analysis.png", plot = p, width = 5, height = 10, dpi = 600)
cat("✓ Final Vertical Spatial Plot saved.\n")
