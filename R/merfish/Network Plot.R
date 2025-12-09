# --------------------------------------------------------------------
# MERFISH - Figure 5: Functional Validation of scDiagnostics
# Network with Controlled Layout + FULL SCORE PRINTING
# --------------------------------------------------------------------

library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(ggraph) 
library(tidygraph)
library(scDiagnostics) 

# --- 1. Load Data ---
cat("--- Loading Data ---\n")
healthy_data <- readRDS("data/merfish/healthy_data.rds")
dss9_data <- readRDS("data/merfish/dss9_data.rds")
healthy_data$tier2 <- as.character(healthy_data$tier2)
dss9_data$tier2 <- as.character(dss9_data$tier2)
set.seed(0)

# --- 2. Run Anomaly Detection ---
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

anomalous_barcodes <- unlist(lapply(anomaly_output, function(res) {
    if(is.null(res)) return(NULL)
    names(res$query_anomaly[res$query_anomaly == TRUE])
}))
anomalous_barcodes <- gsub("Query_", "", anomalous_barcodes)

# --- 3. Define Senders and Receivers ---
meta <- colData(dss9_data) %>% as.data.frame()
meta$barcode <- rownames(meta)
meta$node_group <- "Other"

# Senders (Algorithmic)
meta$node_group[meta$analysis_class == "Fibroblast_Lineage"] <- "Typical Fibro"
meta$node_group[meta$barcode %in% anomalous_barcodes] <- "Anomalous Fibro"

# Receivers (Spatial/Biological)
meta$node_group[grepl("Neutrophil", meta$tier2)] <- "Neutrophil"
meta$node_group[grepl("IAE|Repair associated|Epithelial \\(Clu\\+\\)", meta$tier2)] <- "Repair Epithelium"
meta$node_group[grepl("Stem|TA", meta$tier2)] <- "Stem/Crypt Base"

target_groups <- c("Anomalous Fibro", "Typical Fibro", 
                   "Neutrophil", "Repair Epithelium", "Stem/Crypt Base")
sce_subset <- dss9_data[, meta$node_group %in% target_groups]
sce_subset$node_group <- meta$node_group[meta$node_group %in% target_groups]

# --- 4. The Data-Driven Database ---
interactions_db <- data.frame(
    Ligand = c("Il11", "Il1b", "Cxcl1", "Cxcl5", "Cxcl10", "Rspo1", "Bmp2", "Wnt5a", "Wnt2b"),
    Receptor = c("Il11ra", "Il1r1", "Cxcr2", "Cxcr2", "Cxcr3", "Lgr5", "Bmpr1a", "Fzd5", "Fzd4"),
    Pair_Name = c("Il11-Il11ra", "Il1b-Il1r1", "Cxcl1-Cxcr2", "Cxcl5-Cxcr2", "Cxcl10-Cxcr3", "Rspo1-Lgr5", "Bmp2-Bmpr1a", "Wnt5a-Fzd5", "Wnt2b-Fzd4")
)

# --- 5. Logic, Reporting, AND PRINTING ---
edges_list <- list()
all_scores_list <- list() # New list to track every calculation
receivers <- c("Neutrophil", "Repair Epithelium", "Stem/Crypt Base")

cat("\n--- CALCULATING SCORES ---\n")

for (rcv in receivers) {
    cells_rcv <- sce_subset[, sce_subset$node_group == rcv]
    
    for (i in 1:nrow(interactions_db)) {
        lig <- interactions_db$Ligand[i]
        rec <- interactions_db$Receptor[i]
        pair <- interactions_db$Pair_Name[i]
        
        # Check if genes exist in dataset
        if(!lig %in% rownames(sce_subset) | !rec %in% rownames(sce_subset)) next
        
        # Check receptor expression
        rec_expr <- mean(logcounts(cells_rcv)[rec, ])
        if(rec_expr < 0.0001) next 
        
        # Calculate Scores
        cells_anom <- sce_subset[, sce_subset$node_group == "Anomalous Fibro"]
        cells_typ  <- sce_subset[, sce_subset$node_group == "Typical Fibro"]
        
        score_anom <- mean(logcounts(cells_anom)[lig, ]) * rec_expr * 100
        score_typ  <- mean(logcounts(cells_typ)[lig, ])  * rec_expr * 100
        
        # Determine Winner
        winner <- "None"
        if (score_anom > 1.5 * score_typ & score_anom > 0.005) {
            winner <- "Anomalous"
            # Add to Edge List for Plotting
            edges_list[[length(edges_list)+1]] <- data.frame(from="Anomalous Fibro", to=pair, color="Anomalous Fibro", weight=score_anom, type="Sender")
            edges_list[[length(edges_list)+1]] <- data.frame(from=pair, to=rcv, color="Anomalous Fibro", weight=score_anom, type="Receiver")
            
        } else if (score_typ > 1.5 * score_anom & score_typ > 0.005) {
            winner <- "Typical"
            # Add to Edge List for Plotting
            edges_list[[length(edges_list)+1]] <- data.frame(from="Typical Fibro", to=pair, color="Typical Fibro", weight=score_typ, type="Sender")
            edges_list[[length(edges_list)+1]] <- data.frame(from=pair, to=rcv, color="Typical Fibro", weight=score_typ, type="Receiver")
        }
        
        # --- SAVE SCORES FOR PRINTING ---
        all_scores_list[[length(all_scores_list)+1]] <- data.frame(
            Receiver = rcv,
            Pair = pair,
            Anom_Score = round(score_anom, 4),
            Typ_Score = round(score_typ, 4),
            Winner = winner
        )
    }
}

# --- PRINT THE SCORES ---
all_scores_df <- do.call(rbind, all_scores_list)

cat("\n=======================================================\n")
cat("          FULL INTERACTION SCORE REPORT                \n")
cat("=======================================================\n")
cat(sprintf("%-20s %-15s %-12s %-12s %-10s\n", "Receiver", "Pair", "Anom_Score", "Typ_Score", "Winner"))
cat("-------------------------------------------------------\n")
if (!is.null(all_scores_df) && nrow(all_scores_df) > 0) {
    # Sort by Receiver then by Winner to make it readable
    all_scores_df <- all_scores_df[order(all_scores_df$Receiver, all_scores_df$Winner), ]
    
    for(i in 1:nrow(all_scores_df)) {
        cat(sprintf("%-20s %-15s %-12.4f %-12.4f %-10s\n", 
                    all_scores_df$Receiver[i], 
                    all_scores_df$Pair[i], 
                    all_scores_df$Anom_Score[i], 
                    all_scores_df$Typ_Score[i], 
                    all_scores_df$Winner[i]))
    }
} else {
    cat("No valid interactions found (Receptor expression too low everywhere).\n")
}
cat("=======================================================\n\n")


# --- 6. Graphing with Manual Sorting ---
edges_df <- do.call(rbind, edges_list)

if (is.null(edges_df) || nrow(edges_df) == 0) stop("No significant interactions found to plot.")

all_nodes <- unique(c(edges_df$from, edges_df$to))
nodes_df <- data.frame(name = all_nodes)

nodes_df$col_type <- case_when(
    nodes_df$name %in% c("Typical Fibro", "Anomalous Fibro") ~ "Sender",
    nodes_df$name %in% receivers ~ "Receiver",
    TRUE ~ "Pair"
)

# A. X-Axis (Columns)
nodes_df$x <- case_when(nodes_df$col_type=="Sender"~1, nodes_df$col_type=="Pair"~2, nodes_df$col_type=="Receiver"~3)

# B. Y-Axis (Manual Ordering)
top_group <- c("Neutrophil", "Anomalous Fibro", "Cxcl1-Cxcr2", "Cxcl5-Cxcr2", "Cxcl10-Cxcr3", "Il11-Il11ra", "Il1b-Il1r1")
mid_group <- c("Repair Epithelium")
bot_group <- c("Stem/Crypt Base", "Typical Fibro", "Rspo1-Lgr5", "Bmp2-Bmpr1a", "Wnt5a-Fzd5", "Wnt2b-Fzd4")

nodes_df <- nodes_df %>% 
    mutate(sort_val = case_when(
        name %in% top_group ~ 3,
        name %in% mid_group ~ 2,
        name %in% bot_group ~ 1,
        TRUE ~ 1.5
    )) %>%
    arrange(sort_val) %>% 
    group_by(col_type) %>% 
    mutate(y = row_number()) %>% 
    mutate(y = y - mean(y)) %>% 
    ungroup()

col_pal <- c(
    "Anomalous Fibro" = "#D9230F", 
    "Typical Fibro"   = "#228B22", 
    "Neutrophil"      = "#1F78B4",        
    "Repair Epithelium" = "#FDBF6F",      
    "Stem/Crypt Base" = "#B2DF8A"         
)

gr <- tbl_graph(nodes=nodes_df, edges=edges_df)

p <- ggraph(gr, layout="manual", x=nodes_df$x, y=nodes_df$y) +
    geom_edge_diagonal(aes(edge_width=weight, color=color), alpha=0.7, strength=1) +
    scale_edge_width_continuous(range=c(0.5, 3), guide="none") +
    scale_edge_color_manual(values=col_pal, guide="none") +
    geom_node_point(aes(color=name, size=col_type)) +
    scale_size_manual(values=c("Sender"=8, "Receiver"=8, "Pair"=0), guide="none") +
    scale_color_manual(values=col_pal, na.value="grey50", guide="none") +
    geom_node_text(aes(label=ifelse(col_type!="Pair", name, ""), hjust=ifelse(x==1, 1.1, -0.1)), fontface="bold", size=5) +
    geom_node_label(aes(label=ifelse(col_type=="Pair", name, ""), filter=col_type=="Pair"), 
                    size=2.8, fontface="italic", fill="white", label.size=0.1) +
    theme_void() + coord_cartesian(clip="off", xlim=c(0.5, 3.5)) +
    labs(title="Algorithmic Identification of Pathogenic Signaling", 
         subtitle="scDiagnostics-identified Anomalous Fibroblasts drive Il11/Cxcl inflammation")

print(p)
ggsave("figures/merfish/network_analysis.png", plot=p, width=11, height=7)