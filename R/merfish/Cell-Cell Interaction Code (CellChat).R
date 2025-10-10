# -----------------------------------------------
# MERFISH Mouse Colon IBD - Figure 3 Cell-Cell Interaction
# -----------------------------------------------

# Load libraries
library(SingleCellExperiment)
library(CellChat)
library(patchwork)
library(ggplot2)

# -----------------------------------------------

# Load the processed DSS9 data object (which should have 'anomaly_status')
# If you don't have it in your environment, re-run the Figure 1/2 script to generate it.
dss9_data <- readRDS("data/merfish/dss9_data.rds")
# Re-run the anomaly detection if needed to add the 'anomaly_status' column.

# __________________________________________________________________
# 1. Prepare Data for CellChat
# __________________________________________________________________

# --- Step 1a: Create the "Interaction Group" Label ---
# This is the most crucial step. We combine lineage and anomaly status.

# We will focus on Fibroblasts, Epithelial cells, and key Immune cells.
# All other cell types will be grouped into "Other" to keep the analysis clean.
dss9_data$interaction_group <- as.character(dss9_data$tier2) # Start with the ground truth

# Create the specific "Anomalous" and "Typical" fibroblast groups
dss9_data$interaction_group[grepl("Fibro", dss9_data$tier2) & dss9_data$anomaly_status == "Anomalous"] <- "Fibroblast (Anomalous)"
dss9_data$interaction_group[grepl("Fibro", dss9_data$tier2) & dss9_data$anomaly_status == "Typical"] <- "Fibroblast (Typical)"

# Also create groups for the key receivers
dss9_data$interaction_group[grepl("IAF", dss9_data$tier2)] <- "Fibroblast (Anomalous)" # Group true IAFs with the anomalous ones
dss9_data$interaction_group[grepl("Macrophage", dss9_data$tier2)] <- "Macrophage"
dss9_data$interaction_group[grepl("Monocyte", dss9_data$tier2)] <- "Monocyte"
dss9_data$interaction_group[grepl("Neutrophil", dss9_data$tier2)] <- "Neutrophil"

# Simplify epithelial cells
dss9_data$interaction_group[grepl("Colonocytes|IAE|Stem", dss9_data$tier2)] <- "Epithelial"

# Group all other cell types into "Other"
main_groups <- c("Fibroblast (Anomalous)", "Fibroblast (Typical)", "Macrophage", "Monocyte", "Neutrophil", "Epithelial")
dss9_data$interaction_group[!dss9_data$interaction_group %in% main_groups] <- "Other"

# Make it a factor for CellChat
dss9_data$interaction_group <- factor(dss9_data$interaction_group)

cat("Created 'interaction_group' labels. Distribution:\n")
print(table(dss9_data$interaction_group))


# --- Step 1b: Create the CellChat Object ---
# CellChat needs the expression data and the cell metadata with our new labels.
# We will use the raw 'counts' matrix, as recommended by CellChat.
cellchat <- createCellChat(
    object = counts(dss9_data), 
    meta = as.data.frame(colData(dss9_data)), 
    group.by = "interaction_group"
)

# --- Step 1c: Set up the Ligand-Receptor Database ---
# We are using mouse data, so we load the mouse database.
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB

# Pre-process the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat) # Subset the expression data of signaling genes
# Future::plan("multisession", workers = 4) # Optional: for parallelization
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


# __________________________________________________________________
# 2. Run Cell-Cell Interaction Inference
# __________________________________________________________________

cat("\nInferring cell-cell communication network...\n")
# This is the main calculation step.
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10) # Filter out weak interactions
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cat("Inference complete.\n")


# __________________________________________________________________
# 3. Create the Visualization (Two-Panel Figure)
# __________________________________________________________________

# --- Set the groups we want to compare ---
group.new <- levels(cellchat@idents)
# We will compare signaling FROM the two fibroblast groups TO everyone else.
sender_groups <- c("Fibroblast (Anomalous)", "Fibroblast (Typical)")
receiver_groups <- levels(cellchat@idents)


# --- Panel A: The "Gained Interactions" Dot Plot ---
# This plot will show pathways that are stronger in the "Anomalous" group.

# We need to identify the significantly differential pathways
cellchat <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE, 
                    comparison = sender_groups)

# Get the differential results for visualization
net <- cellchat@net.sum

# Create the dot plot
panel_a <- netVisual_bubble(
    cellchat, 
    sources.use = "Fibroblast (Anomalous)", 
    targets.use = receiver_groups,
    comparison = sender_groups, 
    max.dataset = 2, # Show pathways stronger in "Typical" (dataset 2)
    title.name = "Pathways weaker in Anomalous Fibroblasts",
    angle.x = 45, 
    remove.isolate = TRUE
) + ggtitle("Pathways from Typical Fibroblasts")

panel_b <- netVisual_bubble(
    cellchat, 
    sources.use = "Fibroblast (Anomalous)", 
    targets.use = receiver_groups,
    comparison = sender_groups, 
    max.dataset = 1, # Show pathways stronger in "Anomalous" (dataset 1)
    title.name = "Pathways stronger in Anomalous Fibroblasts",
    angle.x = 45, 
    remove.isolate = TRUE
) + ggtitle("Pathways from Anomalous Fibroblasts")

# Combine the two dot plots for a full comparison
diff_dot_plot <- panel_a + panel_b


# --- Panel B: A "Spotlight" Circle Plot on the CXCL Pathway ---
# The CXCL pathway is key for recruiting neutrophils.

# Get all pathways
pathways.show <- cellchat@netP$pathways
# Create the plot, highlighting the CXCL pathway
panel_c <- netVisual_aggregate(
    cellchat, 
    signaling = "CXCL", 
    layout = "circle",
    sources.use = sender_groups,
    targets.use = "Neutrophil" # Focus on the key target
) + ggtitle("CXCL Signaling to Neutrophils")


# --- Combine and Save the Final Figure ---
# We will arrange the dot plots and the circle plot
final_interaction_figure <- diff_dot_plot / panel_c + 
    plot_layout(heights = c(2, 1.5)) +
    plot_annotation(
        title = "Figure 3: Anomalous fibroblasts adopt a pro-inflammatory communication profile",
        theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
    )

print(final_interaction_figure)

# Save the figure
ggsave("MERFISH_Figure3_Interactions.pdf", 
       plot = final_interaction_figure, 
       width = 12, 
       height = 10, 
       units = "in")