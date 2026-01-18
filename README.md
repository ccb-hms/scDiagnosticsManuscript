# scDiagnostics: Manuscript Code and Analysis

This repository contains scripts, tutorials, and workflows for the manuscript "*scDiagnostics: systematic assessment of cell type annotation in single-cell transcriptomics data*" (Christidis et al.), which introduces the [scDiagnostics](https://bioconductor.org/packages/scDiagnostics) R/Bioconductor package.

The primary goal of this repository is to provide a transparent and reproducible demonstration of how [scDiagnostics](https://bioconductor.org/packages/scDiagnostics) can be used to assess automated cell type annotations from popular annotation tools across different single-cell technologies and experimental conditions.

**Visit the [project website](https://ccb-hms.github.io/scDiagnosticsManuscript) for a comprehensive walkthrough of the analysis and results.**

---

### Overview

We demonstrate the utility of `scDiagnostics` using a simulated and two different real-world single-cell datasets.

1.  **Simulated Data:** A simulation using [splatter](https://bioconductor.org/splatter) to illustrate common challenges of reference-based annotation transfer in a controlled environment with known ground truth composition,
2.  **COVID-19 PBMC scRNA-seq Data:** A sequencing-based single-cell dataset used to demonstrate how the package facilitates the discovery and characterization of a disease-associated cell state in COVID-19.
3.  **MERFISH Mouse Colitis Spatial Transcriptomics Data:** A spatially-resolved imaging-based single-cell dataset used to evaluate annotations where both cell state and spatial organization are key factors in a mouse model of colitis.


For both real-world datasets, we predict cell type labels using a selection of four popular annotation tools:
*   **[SingleR](https://bioconductor.org/packages/SingleR)** (R)
*   **[Azimuth](azimuth.hubmapconsortium.org)** (R)
*   **[CellTypist](https://www.celltypist.org)** (Python)
*   **[scVI](https://scvi-tools.org) / [scArches](https://docs.scarches.org/en/latest/)** (Python)

### Analysis Workflow and Repository Structure

The repository is structured to reflect a logical, multi-step analysis pipeline. The core workflow for each dataset is:
1.  **Data Preprocessing:** Initial import, QC, filtering, and normalization.
2.  **Cell Type Annotation:** Running SingleR, Azimuth, CellTypist and scVI/scArches to predict cell type labels.
3.  **Diagnostic Analysis:** Applying [scDiagnostics](https://bioconductor.org/packages/scDiagnostics) to compare the tool-generated labels against the author's ground-truth annotations.

The main directories are organized as follows:

*   `/R`: Contains R scripts, organized by dataset (`covid`, `merfish`) and function (`auxiliary`, `simulations`). This is the primary location for data preprocessing, cell type annotation (for R-based tools), and diagnostic analysis code.
*   `/python`: Contains Python scripts for predicting cell type annotations with `CellTypist` and `scVI/scArches`.
*   `/data`: Stores pre-computed input files (e.g., `.h5ad`), intermediate objects, and the obtained cell type annotations for each tool.
*   `/website`: Contains Quarto files for reproducible reporting that generate the companion [project website](https://ccb-hms.github.io/scDiagnosticsManuscript/). 

### Reproducing the Analysis

To reproduce the analysis, follow these general steps:

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/ccb-hms/scDiagnosticsManuscript.git
    cd scDiagnosticsManuscript
    ```

2.  **Set up Environments:**
    The analysis workflow has both R and Python dependencies. Required R packages include [scDiagnostics](https://bioconductor.org/packages/scDiagnostics), [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment), [SingleR](https://bioconductor.org/packages/SingleR), and the [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html). Required Python packages include [celltypist](https://pypi.org/project/celltypist/) and [scArches](https://pypi.org/project/scArches/).

3.  **Execute the Workflow:**
    Run the data processing and annotation scripts located in the `/R` and `/python` folders.

### Citation

If you use the code or data from this repository, we kindly ask that you cite our manuscript (details to be added upon publication).
