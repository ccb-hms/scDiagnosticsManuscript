# -----------------------------------------------------
# MERFISH Mouse Colon IBD - Network Interaction Figure
# -----------------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(ggraph) 
library(tidygraph)
library(scDiagnostics) 
library(patchwork)

# -----------------------------------------------------

# __________
# Load Data
# __________

cat("--- Loading Data ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")

# Setting the seed
set.seed(0)

# __________________
# Anomaly Detection
# __________________

cat("--- Running Anomaly Detection ---\n")
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

anomalous_barcodes <- gsub("Query_", "", names(anomaly_output$Fibroblast_Lineage$query_anomaly[anomaly_output$Fibroblast_Lineage$query_anomaly == TRUE]))

# __________________
# Define Metadata
# __________________

meta <- colData(dss9_data) |> as.data.frame()
meta$barcode <- rownames(meta)
meta$node_group <- "Other"
meta$node_group[meta$analysis_class == "Fibroblast_Lineage"] <- "Typical Fibro" 
meta$node_group[meta$barcode %in% anomalous_barcodes] <- "Anomalous Fibro"
meta$node_group[grepl("Neutrophil", meta$tier2)] <- "Neutrophil"
meta$node_group[grepl("IAE|Repair associated|Epithelial \\(Clu\\+\\)", meta$tier2)] <- "IAE"
meta$node_group[grepl("Stem|TA", meta$tier2)] <- "Stem/Crypt Base"

target_groups <- c("Anomalous Fibro", "Typical Fibro", "Neutrophil", "IAE", "Stem/Crypt Base")
sce_subset <- dss9_data[, meta$node_group %in% target_groups]
sce_subset$node_group <- meta$node_group[meta$node_group %in% target_groups]

# ________________________________
# Curate Ligand-Receptor Database
# ________________________________
 
interactions_db <- data.frame(
    Ligand = c("Il11", "Il1b", "Tgfb1", "Cxcl5", "Rspo1", "Wnt5a"),
    Receptor = c("Il11ra", "Il1r1", "Tgfbr2", "Cxcr2", "Lgr5", "Fzd5"),
    Pair_Name = c("Il11-Il11ra", "Il1b-Il1r1", "Tgfb1-Tgfbr2", "Cxcl5-Cxcr2", "Rspo1-Lgr5", "Wnt5a-Fzd5")
)

# __________________
# Calculate Scores
# __________________

receivers <- c("Neutrophil", "IAE", "Stem/Crypt Base")
score_results <- list()

cat("\n--- CALCULATING SCORES ---\n")
for (rcv in receivers) {
    cells_rcv <- sce_subset[, sce_subset$node_group == rcv]
    for (i in 1:nrow(interactions_db)) {
        lig <- interactions_db$Ligand[i]
        rec <- interactions_db$Receptor[i]
        pair <- interactions_db$Pair_Name[i]
        
        if(!lig %in% rownames(sce_subset) | !rec %in% rownames(sce_subset)) next
        rec_expr <- mean(logcounts(cells_rcv)[rec, ])
        if(rec_expr < 0.00001) next 
        
        cells_typ <- sce_subset[, sce_subset$node_group == "Typical Fibro"]
        cells_anom <- sce_subset[, sce_subset$node_group == "Anomalous Fibro"]
        cells_all <- sce_subset[, sce_subset$node_group %in% c("Typical Fibro", "Anomalous Fibro")]
        
        s_typ  <- mean(logcounts(cells_typ)[lig, ]) * rec_expr * 100
        s_anom <- mean(logcounts(cells_anom)[lig, ]) * rec_expr * 100
        s_all  <- mean(logcounts(cells_all)[lig, ]) * rec_expr * 100
        
        score_results[[length(score_results)+1]] <- data.frame(Sender="Typical Fibro", Receiver=rcv, Pair=pair, Raw_Score=s_typ)
        score_results[[length(score_results)+1]] <- data.frame(Sender="Anomalous Fibro", Receiver=rcv, Pair=pair, Raw_Score=s_anom)
        score_results[[length(score_results)+1]] <- data.frame(Sender="All Fibroblasts", Receiver=rcv, Pair=pair, Raw_Score=s_all)
    }
}

# _________________________
# Quadratic Normalization
# _________________________

all_scores_df <- do.call(rbind, score_results)
all_scores_df <- all_scores_df |>
    group_by(Pair, Receiver) |>
    mutate(Max_Pair_Score = max(Raw_Score)) |>
    ungroup() |>
    mutate(
        Ratio = ifelse(Max_Pair_Score > 0, Raw_Score / Max_Pair_Score, 0),
        Plot_Alpha = Ratio^2,           # Strong contrast
        Plot_Width = 0.5 + (Ratio * 0.7) # Subtle thickness
    )

# __________________
# Plotting Function
# __________________

col_pal <- c("All Fibroblasts"="grey40", "Typical Fibro"="#228B22", "Anomalous Fibro"="#D9230F", 
             "Neutrophil"="#1F78B4", "IAE"="#FDBF6F", "Stem/Crypt Base"="#008B8B")

create_graph_object <- function(sender_name, title_text) {
    
    subset_edges <- all_scores_df |> filter(Sender == sender_name)
    
    unique_pairs <- unique(subset_edges$Pair)
    receivers_present <- unique(subset_edges$Receiver)
    nodes_df <- data.frame(name = unique(c(sender_name, unique_pairs, receivers_present)))
    
    nodes_df$col_type <- case_when(
        nodes_df$name == sender_name ~ "Sender",
        nodes_df$name %in% receivers ~ "Receiver",
        TRUE ~ "Pair"
    )
    
    nodes_df$x <- case_when(nodes_df$col_type=="Sender"~1, nodes_df$col_type=="Pair"~2, nodes_df$col_type=="Receiver"~3)
    
    # Manual Y Ordering
    top_group <- c("Neutrophil", "Anomalous Fibro", "All Fibroblasts", "Il11-Il11ra", "Tgfb1-Tgfbr2", "Il1b-Il1r1")
    mid_group <- c("IAE", "Cxcl5-Cxcr2")
    bot_group <- c("Stem/Crypt Base", "Typical Fibro", "Rspo1-Lgr5", "Wnt5a-Fzd5")
    
    nodes_df <- nodes_df |> 
        mutate(sort_val = case_when(name %in% top_group ~ 3, name %in% mid_group ~ 2, name %in% bot_group ~ 1, TRUE ~ 1.5)) |>
        arrange(sort_val) |> group_by(col_type) |> mutate(y = row_number()) |> mutate(y = y - mean(y)) |> ungroup()
    
    edges_to_plot <- data.frame()
    for(i in 1:nrow(subset_edges)) {
        s <- subset_edges$Sender[i]
        p <- subset_edges$Pair[i]
        r <- subset_edges$Receiver[i]
        w <- subset_edges$Plot_Width[i] 
        a <- subset_edges$Plot_Alpha[i] 
        
        edges_to_plot <- rbind(edges_to_plot, 
                               data.frame(from=s, to=p, width=w, alpha=a, color=sender_name),
                               data.frame(from=p, to=r, width=w, alpha=a, color=sender_name))
    }
    
    gr <- tbl_graph(nodes=nodes_df, edges=edges_to_plot)
    
    p <- ggraph(gr, layout="manual", x=nodes_df$x, y=nodes_df$y) +
        geom_edge_diagonal(aes(edge_width=width, edge_alpha=alpha, color=color), strength=1) +
        scale_edge_width_identity() +
        scale_edge_alpha_identity() +
        scale_edge_color_manual(values=col_pal, guide="none") +
        
        # Nodes: Enable Guide (Legend)
        geom_node_point(aes(color=name, size=col_type)) +
        scale_size_manual(values=c("Sender"=6, "Receiver"=6, "Pair"=0), guide="none") +
        scale_color_manual(values=col_pal, 
                           name = "Cell Type / State", # Legend Title
                           breaks = c("Typical Fibro", "Anomalous Fibro", "All Fibroblasts", 
                                      "Neutrophil", "IAE", "Stem/Crypt Base"),
                           guide = guide_legend(override.aes = list(size=5))) +
        
        # Keep Gene Pair Labels (Center)
        geom_node_label(aes(label=ifelse(col_type=="Pair", name, ""), filter=col_type=="Pair"), 
                        size=5,                  
                        fontface="bold.italic",  
                        fill="white", 
                        label.size=0.25) +
        
        theme_void() + 
        coord_cartesian(clip="off", xlim=c(0.5, 3.5)) +
        labs(title = title_text) +
        theme(
            plot.title = element_text(hjust=0.5, face="bold", size=14),
            legend.position = "right",
            legend.title = element_text(face="bold"),
            legend.text = element_text(size=10)
        )
    
    return(p)
}

# ________________
# Combine & Save
# ________________

p1 <- create_graph_object("All Fibroblasts", "A. All Fibroblasts (Avg)")
p2 <- create_graph_object("Typical Fibro", "B. Typical Fibro (Homeostatic)")
p3 <- create_graph_object("Anomalous Fibro", "C. Anomalous Fibro (Inflamed)")

# Collect guides ensures one unified legend
final_plot <- p1 + p2 + p3 + 
    plot_layout(ncol = 3, guides = "collect") + 
    plot_annotation(
        title = "Functional Switch: Homeostasis vs. Inflammation",
        subtitle = "Quadratic scaling of opacity highlights dominant signaling pathways.",
        theme = theme(plot.title = element_text(size=18, face="bold", hjust=0.5), 
                      plot.subtitle = element_text(size=12, hjust=0.5))
    )

print(final_plot)
ggsave("figures/merfish/network_analysis.png", plot=final_plot, width=26, height=8, dpi=600)
cat("\n✓ Saved clean network plot with legend to figures/merfish/network_analysis.png\n")