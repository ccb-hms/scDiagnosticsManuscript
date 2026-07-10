# ---------------------------------------------
# COVID-19 PBMC - Supplementary Figure 1 Code
# ---------------------------------------------

library(zellkonverter)
library(SingleCellExperiment)
library(HDF5Array)
library(Matrix)
library(scDiagnostics)
library(bench)
library(ggplot2)
library(dplyr)
library(patchwork)

# ---------------------------------------------

# _____________________________
# 1. SETUP & DELAYED LOADING
# _____________________________

input_path <- "data/covid/covid_data.h5ad"
message("Loading raw h5ad file: ", input_path)

# Load using HDF5-backed mode to save RAM!
sce_full <- readH5AD(input_path, use_hdf5 = TRUE, reader = "R")

# Minimal Cleaning (Filter out "respiratory system disorder")
keep_disease <- sce_full$disease %in% c("normal", "COVID-19")
sce_filtered <- sce_full[, keep_disease]

# Clear memory from the full object
rm(sce_full)
gc()

# Get total available pool (~624,000 cells)
available_indices <- seq_len(ncol(sce_filtered))
message(sprintf("Total available cells for benchmarking: %d", length(available_indices)))

# Sizes to benchmark
cell_sizes <- c(10000, 50000, 100000, 250000, 500000)
benchmark_results <- data.frame()

# ______________________
# 2. BENCHMARKING LOOP
# ______________________

set.seed(0)

for (n_cells in cell_sizes) {
    message(sprintf("\n=== Benchmarking at N = %d total cells ===", n_cells))
    
    # 1. Randomly sample the desired number of cells
    sampled_idx <- sample(available_indices, n_cells)
    
    # 2. Arbitrarily assign 30% to Ref, 70% to Query to mimic real-world mapping
    n_ref <- floor(n_cells * 0.3)
    ref_idx <- sampled_idx[1:n_ref]
    query_idx <- sampled_idx[(n_ref + 1):length(sampled_idx)]
    
    # 3. Realize into Memory as Sparse Matrices (Crucial for CPU timing)
    message("  > Realizing sparse matrices to RAM (No Metadata)...")
    counts_ref <- as(realize(assay(sce_filtered[, ref_idx], "X")), "CsparseMatrix")
    counts_query <- as(realize(assay(sce_filtered[, query_idx], "X")), "CsparseMatrix")
    
    # 4. Create CLEAN, stripped-down SCE objects (No old PCA/UMAP baggage)
    ref_sce <- SingleCellExperiment(
        assays = list(logcounts = counts_ref),
        colData = DataFrame(author_cell_type = sce_filtered$author_cell_type[ref_idx])
    )
    query_sce <- SingleCellExperiment(
        assays = list(logcounts = counts_query),
        colData = DataFrame(author_cell_type = sce_filtered$author_cell_type[query_idx])
    )
    
    # Identify a large cell type to test
    target_cell_type <- names(sort(table(ref_sce$author_cell_type), decreasing = TRUE))[1]
    
    # --- TASK 1: processPCA ---
    message("  > Task 1: processPCA...")
    bm_process <- bench::mark(
        processPCA = { ref_pca <- scDiagnostics::processPCA(ref_sce, n_hvgs = 1000) },
        iterations = 1, check = FALSE, memory = TRUE
    )
    ref_pca <- scDiagnostics::processPCA(ref_sce, n_hvgs = 1000)
    
    # --- TASK 2: projectPCA ---
    message("  > Task 2: projectPCA...")
    bm_project <- bench::mark(
        projectPCA = {
            scDiagnostics::projectPCA(query_data = query_sce, reference_data = ref_pca,
                                      query_cell_type_col = "author_cell_type", ref_cell_type_col = "author_cell_type",
                                      cell_types = target_cell_type, pc_subset = 1:5, max_cells_query = NULL, max_cells_ref = NULL)
        },
        iterations = 1, check = FALSE, memory = TRUE
    )
    
    # --- TASK 3: detectAnomaly ---
    message("  > Task 3: detectAnomaly...")
    bm_if <- bench::mark(
        detectAnomaly = {
            scDiagnostics::detectAnomaly(query_data = query_sce, reference_data = ref_pca, 
                                         query_cell_type_col = "author_cell_type", ref_cell_type_col = "author_cell_type", 
                                         cell_types = target_cell_type, pc_subset = NULL, n_hvgs = 1000, n_tree = 500,
                                         threshold_method = "MAD", mad_multiplier = 2)
        },
        iterations = 1, check = FALSE, memory = TRUE
    )
    
    # --- TASK 4: calculateReconstructionError ---
    message("  > Task 4: calculateReconstructionError...")
    bm_re <- bench::mark(
        calculateReconstructionError = {
            scDiagnostics::calculateReconstructionError(reference_data = ref_sce, query_data = query_sce,
                                                        ref_cell_type_col = "author_cell_type", query_cell_type_col = "author_cell_type",
                                                        cell_types = target_cell_type, pc_subset = 1:5, n_hvgs = 100, mad_multiplier = 2)
        },
        iterations = 1, check = FALSE, memory = TRUE
    )
    
    # --- TASK 5: calculateGeneShifts ---
    message("  > Task 5: calculateGeneShifts...")
    query_pca <- scDiagnostics::processPCA(query_sce, n_hvgs = 1000) # Precompute for anomaly check inside function
    bm_gs <- bench::mark(
        calculateGeneShifts = {
            scDiagnostics::calculateGeneShifts(query_data = query_pca, reference_data = ref_pca, 
                                               query_cell_type_col = "author_cell_type", ref_cell_type_col = "author_cell_type",
                                               cell_types = target_cell_type, pc_subset = 1:5, n_top_loadings = 30,
                                               p_value_threshold = 0.05, adjust_method = "fdr", 
                                               detect_anomalies = TRUE, anomaly_comparison = TRUE, mad_multiplier = 2)
        },
        iterations = 1, check = FALSE, memory = TRUE
    )
    
    # Append results
    benchmark_results <- rbind(benchmark_results, data.frame(
        Total_Cells = n_cells,
        Method = c("processPCA", "projectPCA", "detectAnomaly", "calculateReconstructionError", "calculateGeneShifts"),
        Time_Seconds = c(as.numeric(bm_process$median), as.numeric(bm_project$median), 
                         as.numeric(bm_if$median), as.numeric(bm_re$median), as.numeric(bm_gs$median)),
        Memory_GB = c(as.numeric(bm_process$mem_alloc)/1024^3, as.numeric(bm_project$mem_alloc)/1024^3, 
                      as.numeric(bm_if$mem_alloc)/1024^3, as.numeric(bm_re$mem_alloc)/1024^3, as.numeric(bm_gs$mem_alloc)/1024^3)
    ))
    
    # Free up RAM before the next massive loop
    rm(ref_sce, query_sce, counts_ref, counts_query, ref_pca, query_pca)
    gc()
}

print(benchmark_results)

# _____________________________
# 3. VISUALIZATION (Figure S1)
# _____________________________

my_colors <- c("processPCA" = "#7F8C8D", "projectPCA" = "#F39C12", 
               "detectAnomaly" = "#4C72B0", "calculateReconstructionError" = "#C44E52", "calculateGeneShifts" = "#9B59B6")

p_time <- ggplot(benchmark_results, aes(x = Total_Cells, y = Time_Seconds, color = Method)) +
    geom_line(linewidth = 1.2) + geom_point(size = 3) +
    scale_color_manual(values = my_colors) +
    scale_x_continuous(labels = scales::comma, breaks = cell_sizes) +
    theme_bw(base_size = 14) + 
    labs(title = "Execution Time Scaling", x = "Total Cells (Reference + Query)", y = "Time (Seconds)") +
    theme(legend.position = "none", panel.grid.minor = element_blank())

p_mem <- ggplot(benchmark_results, aes(x = Total_Cells, y = Memory_GB, color = Method)) +
    geom_line(linewidth = 1.2) + geom_point(size = 3) +
    scale_color_manual(values = my_colors) +
    scale_x_continuous(labels = scales::comma, breaks = cell_sizes) +
    theme_bw(base_size = 14) + 
    labs(title = "Memory Allocation Scaling", x = "Total Cells (Reference + Query)", y = "Peak RAM (GB)") +
    theme(panel.grid.minor = element_blank())

# Combine using patchwork (stacked vertically)
final_plot <- (p_time / p_mem) + 
    plot_layout(guides = "collect") & 
    theme(
        plot.title = element_text(face = "bold", size = 24, hjust = 0.5, margin = margin(b = 15)),
        axis.title = element_text(face = "bold", size = 22),
        axis.text = element_text(size = 18),
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.5, "cm"),
        legend.margin = margin(t = 10) 
    ) &
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) 

print(final_plot)

ggsave("figures/supp/covid/Fig_S1_computing_benchmark.png", final_plot, width = 14, height = 18, dpi = 600)

# ________
# Summary
# ________

print("Supplementary Figure S1 complete!")
print("Saved in: figures/supp/covid/")

