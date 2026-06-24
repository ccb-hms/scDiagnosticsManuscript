# -------------------------------------------------
# COVID-19 PBMC: Practical Scalability Benchmark 
# -------------------------------------------------

library(SingleCellExperiment)
library(scDiagnostics)
library(bench)
library(scran)
library(scater)

# -------------------------------------------------

# _________________
# 1. Load the Data
# _________________

message("Loading COVID-19 dataset...")
normal_data <- readRDS("data/covid/normal_data_sce.rds")
covid_data <- readRDS("data/covid/covid_data_sce.rds")

# Wipe PCA to simulate a completely raw start
reducedDims(normal_data) <- list()
reducedDims(covid_data) <- list()

message(sprintf("\nBenchmarking on Total Cells: %d (Ref: %d, Query: %d)", 
                ncol(normal_data) + ncol(covid_data), ncol(normal_data), ncol(covid_data)))

# ______________________
# 2. RUN THE BENCHMARKS
# ______________________

# Task 1: Process PCA (The Heavy SVD Lifting)
message("\nTask 1: processPCA (Global SVD on Reference)...")
bm_process <- bench::mark(
    processPCA = {
        normal_data_pca <- scDiagnostics::processPCA(normal_data, n_hvgs = 1000)
    },
    iterations = 1, check = FALSE, memory = TRUE
)

# Extract the computed PCA object to use in subsequent tasks
normal_data_pca <- scDiagnostics::processPCA(normal_data, n_hvgs = 1000)

# Task 2: Project PCA (Matrix Multiplication)
message("Task 2: projectPCA (Projecting Query onto Reference)...")
bm_project <- bench::mark(
    projectPCA = {
        scDiagnostics::projectPCA(
            query_data = covid_data, reference_data = normal_data_pca,
            query_cell_type_col = "azimuth_celltype_l1", ref_cell_type_col = "author_cell_type",
            cell_types = "CD14_mono", pc_subset = 1:5, assay_name = "logcounts",
            max_cells_query = NULL, max_cells_ref = NULL
        )
    },
    iterations = 1, check = FALSE, memory = TRUE
)

# Task 3: Isolation Forest Anomaly Detection
message("Task 3: detectAnomaly (Isolation Forest)...")
bm_if <- bench::mark(
    detectAnomaly = {
        scDiagnostics::detectAnomaly(
            query_data = covid_data, reference_data = normal_data_pca, 
            query_cell_type_col = "azimuth_celltype_l1", ref_cell_type_col = "author_cell_type", 
            cell_types = "CD14_mono", pc_subset = NULL, n_hvgs = 1000, n_tree = 500,
            threshold_method = "MAD", mad_multiplier = 2
        )
    },
    iterations = 1, check = FALSE, memory = TRUE
)

# Task 4: Reconstruction Error (Computes LOCAL PCA internally)
message("Task 4: calculateReconstructionError (Local PCA + SSE)...")
bm_re <- bench::mark(
    calculateReconstructionError = {
        scDiagnostics::calculateReconstructionError(
            reference_data = normal_data, query_data = covid_data,
            ref_cell_type_col = "author_cell_type", query_cell_type_col = "azimuth_celltype_l1",
            cell_types = "CD14_mono", pc_subset = 1:5, n_hvgs = 100, mad_multiplier = 2
        )
    },
    iterations = 1, check = FALSE, memory = TRUE
)

# Pre-compute Query PCA (Required for detectAnomaly inside GeneShifts)
message("Pre-computing Query PCA for Task 5...")
covid_data_pca <- scDiagnostics::processPCA(covid_data, n_hvgs = 1000)

# Task 5: Gene Shifts (Automated Top Loading Search + Differential Expression)
message("Task 5: calculateGeneShifts (Automated Search)...")
bm_gs <- bench::mark(
    calculateGeneShifts = {
        scDiagnostics::calculateGeneShifts(
            query_data = covid_data_pca, reference_data = normal_data_pca, 
            query_cell_type_col = "azimuth_celltype_l1", ref_cell_type_col = "author_cell_type",
            cell_types = "CD14_mono", pc_subset = 1:5, n_top_loadings = 30,
            assay_name = "logcounts", p_value_threshold = 0.05, adjust_method = "fdr", 
            detect_anomalies = TRUE, anomaly_comparison = TRUE, mad_multiplier = 2
        )
    },
    iterations = 1, check = FALSE, memory = TRUE
)

# ___________________
# 3. COMPILE RESULTS
# ___________________

benchmark_results <- data.frame(
    Function = c("processPCA", "projectPCA", "detectAnomaly", 
                 "calculateReconstructionError", "calculateGeneShifts"),
    Time_Seconds = c(as.numeric(bm_process$median), as.numeric(bm_project$median), 
                     as.numeric(bm_if$median), as.numeric(bm_re$median), as.numeric(bm_gs$median)),
    Memory_MB = c(as.numeric(bm_process$mem_alloc)/1024^2, as.numeric(bm_project$mem_alloc)/1024^2, 
                  as.numeric(bm_if$mem_alloc)/1024^2, as.numeric(bm_re$mem_alloc)/1024^2, 
                  as.numeric(bm_gs$mem_alloc)/1024^2)
)

message("\n=== MAIN MANUSCRIPT SCALABILITY METRICS ===")
print(benchmark_results)

