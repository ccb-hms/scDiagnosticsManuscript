# scDiagnostics: Manuscript Code and Analysis

This repository contains the full computational workflow for the manuscript "*scDiagnostics: diagnostic tools for assessing cell type annotation quality in scRNA-seq data*" (Christidis et al.), which introduces and demonstrates the `scDiagnostics` R/Bioconductor package.

The primary goal of this repository is to provide a transparent and reproducible demonstration of how `scDiagnostics` can be used to evaluate and diagnose automated cell type annotations from a variety of popular tools across different single-cell technologies and biological contexts.

[**Visit the project website for a full walkthrough of the analysis and results.**](https://ccb-hms.github.io/scDiagnosticsManuscript/)  

---

### Demonstration Overview

The utility of `scDiagnostics` is showcased through a comprehensive analysis of two distinct public datasets, as well as a controlled simulation.

1.  **COVID-19 PBMC scRNA-seq Data:** A dense, high-cell-count dataset used to diagnose annotations in a disease vs. healthy context.
2.  **MERFISH Mouse Colitis Spatial Transcriptomics Data:** A complex spatial dataset used to evaluate annotations where both cell state and spatial organization are key factors in a model of inflammatory bowel disease (IBD).
3.  **Simulated Data:** A `splatter`-based simulation used to illustrate core diagnostic principles in a controlled environment where the ground truth is perfectly known.

For each real-world dataset, we generate annotations using a suite of popular tools to provide a robust test bed for `scDiagnostics`:
*   **SingleR** (R)
*   **Azimuth** (R)
*   **CellTypist** (Python)
*   **scVI/scArches** (Python)

### Workflow and Repository Structure

The repository is structured to reflect a logical, multi-step analysis pipeline. The core workflow for each dataset is:
1.  **Data Preparation:** Initial import, QC, filtering, and normalization.
2.  **Annotation Generation:** Running SingleR, Azimuth, CellTypist and scVI/scArches to produce query labels.
3.  **Diagnostic Analysis:** Applying `scDiagnostics` to compare the tool-generated labels against the author's ground-truth annotations.

The main directories are organized as follows:

*   `/R`: Contains all R scripts, organized by dataset (`covid`, `merfish`) and function (`auxiliary`, `simulations`). This is the primary location for data processing, annotation generation (for R-based tools), and analysis code.
*   `/python`: Contains Python scripts for generating annotations with `CellTypist` and `scVI/scArches`.
*   `/data`: Stores key input files (e.g., `.h5ad`), intermediate objects, and the final annotation results from each tool.
*   `/website`: Contains `_quarto.yml`, `index.qmd`, etc.: Quarto files that generate the companion [project website](https://ccb-hms.github.io/scDiagnosticsManuscript/). 

### Reproducing the Analysis

To reproduce the analysis, follow these general steps:

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/ccb-hms/scDiagnosticsManuscript.git
    cd scDiagnosticsManuscript
    ```

2.  **Set up Environments:**
    The analysis requires both R and Python environments. Key R packages include `scDiagnostics`, `SpatialExperiment`, `SingleR`, and the `tidyverse`. The Python environment requires `celltypist`.

3.  **Execute the Workflow:**
    Run the data processing and annotation scripts found in `/R` and `/python`.

### Citation

If you use the code or data from this repository, we kindly ask that you cite our manuscript (details to be added upon publication).