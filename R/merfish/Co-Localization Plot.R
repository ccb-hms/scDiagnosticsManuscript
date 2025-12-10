# --------------------------------------------------
# MERFISH Mouse Colon IBD - Co-Localization Figure
# --------------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scDiagnostics)
library(FNN) 

# --------------------------------------------------

# __________
# Load Data
# __________

cat("\n--- Loading data ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

# Setting the seed
set.seed(0)

# __________________
# Anomaly Detection
# __________________

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
    anomaly_threshold = 0.5
)

anomalous_barcodes <- unlist(lapply(anomaly_output, function(res) {
    if(is.null(res)) return(NULL)
    names(res$query_anomaly[res$query_anomaly == TRUE])
}))
anomalous_barcodes <- gsub("Query_", "", anomalous_barcodes)

# ________________________________________
# Spatial Calculations: The 3 Key Niches
# ________________________________________

spatial_df <- data.frame(
    barcode = colnames(dss9_data),
    x = spatialCoords(dss9_data)[, "x"],
    y = spatialCoords(dss9_data)[, "y"],
    ground_truth_type = dss9_data$tier2,
    broad_class = dss9_data$analysis_class
)

# 1. Neutrophils (The Immune Trigger)
pattern_neutro <- "Neutrophil"
coords_neutro <- spatial_df |> filter(grepl(pattern_neutro, ground_truth_type, ignore.case = TRUE)) |> select(x, y)

# 2. Inflamed Epithelium (The Damage Target)
pattern_inflamed <- "IAE|Epithelial \\(Clu\\+\\)"
coords_inflamed <- spatial_df |> filter(grepl(pattern_inflamed, ground_truth_type, ignore.case = TRUE)) |> select(x, y)

# 3. Stem Cell Niche (The Homeostatic Source - Negative Control)
pattern_stem <- "Stem|TA"
coords_stem <- spatial_df |> filter(grepl(pattern_stem, ground_truth_type, ignore.case = TRUE)) |> select(x, y)

# Calculate Distances
fibro_df <- spatial_df |> filter(broad_class == "Fibroblast_Lineage")

get_dists <- function(coords) {
    if(nrow(coords) == 0) return(rep(NA, nrow(fibro_df)))
    FNN::get.knnx(data = as.matrix(coords), query = as.matrix(fibro_df[, c("x", "y")]), k = 1)$nn.dist
}

fibro_df$dist_neutro <- get_dists(coords_neutro)
fibro_df$dist_inflamed <- get_dists(coords_inflamed)
fibro_df$dist_stem   <- get_dists(coords_stem)

fibro_df$prediction <- ifelse(fibro_df$barcode %in% anomalous_barcodes, "Anomalous", "Typical")

# _______________________
# Prepare Data Structure 
# _______________________

df_split <- fibro_df |>
    select(barcode, prediction, dist_neutro, dist_inflamed, dist_stem)

df_all <- df_split
df_all$prediction <- "All Fibroblasts"

combined_df <- rbind(df_split, df_all) |>
    tidyr::pivot_longer(
        cols = c(dist_neutro, dist_inflamed, dist_stem), 
        names_to = "target_type", 
        values_to = "distance"
    ) |>
    mutate(panel_label = case_when(
        target_type == "dist_neutro" ~ "Dist. to Neutrophils\n(Acute Inflammation)",
        target_type == "dist_inflamed" ~ "Dist. to Inflamed Epithelium\n(Damage)",
        target_type == "dist_stem"   ~ "Dist. to Stem/Crypt\n(Homeostasis)"
    )) |>
    mutate(panel_label = factor(panel_label, levels = c(
        "Dist. to Neutrophils\n(Acute Inflammation)", 
        "Dist. to Inflamed Epithelium\n(Damage)", 
        "Dist. to Stem/Crypt\n(Homeostasis)"
    )))

# ______________________________
# Manual Smooth CDF Calculation 
# ______________________________

smooth_curves <- combined_df |>
    group_by(panel_label, prediction) |>
    reframe({
        d <- density(distance, from = 0, n = 512)
        cdf <- cumsum(d$y) / sum(d$y)
        data.frame(distance = d$x, cdf = cdf)
    })

smooth_curves$prediction <- factor(smooth_curves$prediction, 
                                   levels = c("Anomalous", "Typical", "All Fibroblasts"))

# ___________________________
# Smart Limits & Statistics 
# ___________________________

calc_limit <- function(type_name) {
    clean_name <- type_name 
    vals <- combined_df$distance[combined_df$panel_label == clean_name]
    ceiling(quantile(vals, 0.90, na.rm=TRUE) / 50) * 50
}

limit_df <- data.frame(
    panel_label = unique(combined_df$panel_label),
    prediction = "Anomalous",
    cdf = 0.5
)
limit_df$distance <- sapply(limit_df$panel_label, calc_limit)

calc_pval <- function(col_name) {
    p <- wilcox.test(fibro_df[[col_name]] ~ fibro_df$prediction, alternative = "two.sided")$p.value
    format.pval(p, digits=2)
}

col_map <- list(
    "Dist. to Neutrophils\n(Acute Inflammation)" = "dist_neutro",
    "Dist. to Inflamed Epithelium\n(Damage)" = "dist_inflamed",
    "Dist. to Stem/Crypt\n(Homeostasis)" = "dist_stem"
)

ann_text <- limit_df |>
    mutate(
        distance = distance * 0.65,
        cdf = 0.15,
        label = sapply(panel_label, function(x) {
            col <- col_map[[as.character(x)]]
            paste0("Wilcoxon\np < ", calc_pval(col))
        })
    )

# ______________
# Generate Plot
# ______________

colors_pred <- c(
    "Anomalous"       = "#D9230F", 
    "Typical"         = "#228B22", 
    "All Fibroblasts" = "#606060"
)

p <- ggplot(smooth_curves, aes(x = distance, y = cdf, color = prediction)) +
    geom_line(linewidth = 1) +
    
    # 3 Panels Horizontal
    facet_wrap(~panel_label, scales = "free_x", nrow = 1) + 
    
    geom_blank(data = limit_df) +
    geom_text(data = ann_text, aes(label = label), 
              color = "black", size = 3.5, fontface = "italic", hjust=0.5) +
    
    scale_color_manual(name = "Population", values = colors_pred) +
    
    labs(
        title = "Co-Localization Analysis: The Inflammatory Niche",
        subtitle = "Anomalous Fibroblasts are spatially coupled with Damage and Immunity",
        x = expression(bold("Distance to Feature ("*mu*"m)")), 
        y = "Cumulative Proportion"
    ) +
    
    theme_classic(base_size = 14) +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.major.y = element_line(color = "grey92"), 
        legend.position = "bottom", 
        plot.title = element_text(face="bold", hjust=0.5),
        plot.subtitle = element_text(hjust=0.5, size=11),
        axis.text.x = element_text(angle=0, hjust=0.5)
    ) +
    coord_cartesian(expand = FALSE, ylim = c(0, 1.05))

print(p)

# Adjusted width for 3 panels (was 16 for 4 panels)
ggsave("figures/merfish/co-localization_analysis.png", plot = p, width = 12, height = 8, dpi = 600)
cat("✓ Final 3-Panel Localization Plot Saved.\n")