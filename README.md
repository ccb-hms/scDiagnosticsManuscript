# scDiagnostics Manuscript

Comprehensive tutorials, analysis code, and reproducible workflows demonstrating [`scDiagnostics`](https://bioconductor.org/packages/scDiagnostics) for systematic assessment of cell type annotation in single-cell transcriptomics data.

**Manuscript:** Christidis, A., Ghazi, A., Chawla, S., Turaga, N., Gentleman, R., & Geistlinger, L. scDiagnostics: systematic assessment of cell type annotation in single-cell transcriptomics data. *Submitted*. [**Read the preprint on bioRxiv**](https://www.biorxiv.org/content/10.64898/2026.01.29.701618v1).

**Analysis and Results:** [**Manuscript Website**](https://ccb-hms.github.io/scDiagnosticsManuscript/)

## Overview

We demonstrate scDiagnostics using a simulated and three different real-world single-cell datasets:

**1. Simulated single-cell data (splatter)**

- Synthetic data with known cell types and ground truth composition
- Use case: Illustration of common challenges of reference-based annotation transfer in a controlled setting

**2. Zeisel Mouse Brain Benchmarking**

- Mouse cortex and hippocampus scRNA-seq dataset (Zeisel et al., 2015)
- Use case: Systematic stress-testing of diagnostic sensitivity and specificity against label noise, class imbalance, and batch effects

**3. COVID-19 PBMC scRNA-seq**

- Single-cell RNA-seq data from severe COVID-19 patients and healthy controls
- Source: [CZI CELLxGENE](https://cellxgene.cziscience.com/collections/ddfad306-714d-4cc0-9985-d9072820c530) (Stephenson et al., 2021)
- Use case: Discovery and characterization of disease-associated immune cell states

**4. MERFISH Mouse Colitis**

- Spatial transcriptomics from a mouse model of inflammatory bowel disease
- Source: [MerfishData Bioconductor package](https://bioconductor.org/packages/MerfishData) (Cadinu et al., 2024)
- Use case: Spatial validation of annotation quality and disease-associated cell states

For each dataset, we predict cell type labels using four popular annotation tools:

- [**Azimuth**](https://github.com/satijalab/azimuth) — Weighted k-NN mapping
- [**SingleR**](https://bioconductor.org/packages/SingleR) — Correlation-based assignment
- [**CellTypist**](https://www.celltypist.org/) — Machine learning classifier
- [**scVI/scArches**](https://docs.scvi-tools.org/) — Deep learning with VAE (GPU-accelerated)

## Quick Start

### Installation

```{r}
source("R/covid/R_Package_Installation_Pipeline.R")
```

Or for MERFISH:

```{r}
source("R/merfish/R_Package_Installation_Pipeline.R")
```

### Download Data

All pre-processed datasets with annotations are available on [Zenodo](https://doi.org/10.5281/zenodo.18274941):

```{r}
source("data/downloadData.R")
downloadData()
```

This automated script downloads all four SingleCellExperiment/SpatialExperiment objects into your `data/covid/` and `data/merfish/` directories. For manual download, visit the [Zenodo repository](https://doi.org/10.5281/zenodo.18274941).

See detailed instructions: [Setup & Installation](setup.qmd), [Accessing Data](data-access.qmd)

## Documentation

Full tutorials and analysis code available at [https://ccb-hms.github.io/scDiagnosticsManuscript/](https://ccb-hms.github.io/scDiagnosticsManuscript/):

### Setup & Methods

Analysis environment setup, data retrieval, and reproducible analysis workflows:

- **[Setup & Installation](setup.qmd)** — Install R and Python dependencies (GPU recommended for scVI/scArches)
- **[Accessing Data](data-access.qmd)** — Download pre-processed datasets from Zenodo
- **[Cell type annotation](adding-annotations.qmd)** — Apply all four annotation methods to query data

### Tutorials & Workflows

Quick start, core functionality, and common analysis workflows:

- **[scDiagnostics Overview](overview-scdiagnostics.qmd)** — Introduction to diagnostic framework and key concepts
- **[Simulation Analysis](results-simulation.qmd)** — Demonstration of common challenges in reference-based annotation transfer
- **[Zeisel Brain Benchmarking](results-brain.qmd)** — Quantitative benchmarking of scDiagnostics under stress-tested single-cell scenarios
- **[COVID-19 Analysis](results-covid.qmd)** — Annotation assessment and anomaly detection in scRNA-seq data
- **[MERFISH Analysis](results-merfish.qmd)** — Annotation assessment and anomaly detection in spatial transcriptomics data
- **[Exploring Annotation Tool Diagnostics](annotation-diagnostics.qmd)** — Complementary aspects of scDiagnostics and built-in quality metrics from major annotation tools

## Installation

### R Packages

All required R packages are automatically installed by running:

```{r}
source("R/covid/R_Package_Installation_Pipeline.R")
source("R/merfish/R_Package_Installation_Pipeline.R")
```

### Python Environment

For GPU-accelerated scVI/scArches annotation:

```bash
conda env create -f environment-scvi.yml
conda activate scvi-env
```

See [Setup & Installation](setup.qmd) for detailed instructions.

## Citation

If you use this code, data, or analyses, please cite:

```bibtex
@article{christidis2024scDiagnostics,
  author = {Christidis, A. and Ghazi, A. and Chawla, S. and Turaga, N. and Gentleman, R. and Geistlinger, L.},
  title = {scDiagnostics: systematic assessment of cell type annotation in single-cell transcriptomics data},
  year = {2026},
  url = {https://www.biorxiv.org/content/10.64898/2026.01.29.701618v1},
  note = {Submitted}
}
```

## Repository

**Code and Scripts:** [github.com/ccb-hms/scDiagnosticsManuscript](https://github.com/ccb-hms/scDiagnosticsManuscript)

**Data:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18274942.svg)](https://doi.org/10.5281/zenodo.18274942)

## Contact

For questions or feedback, please [open an issue](https://github.com/ccb-hms/scDiagnosticsManuscript/issues) on GitHub.