# scDiagnostics Manuscript

Comprehensive tutorials, analysis code, and reproducible workflows demonstrating [`scDiagnostics`](https://bioconductor.org/packages/scDiagnostics) for systematic assessment of cell type annotation in single-cell transcriptomics data.

**Manuscript:** Christidis, A., Ghazi, A., Chawla, S., Turaga, N., Gentleman, R., & Geistlinger, L. scDiagnostics: systematic assessment of cell type annotation in single-cell transcriptomics data. *Submitted*.

**Website:** [https://ccb-hms.github.io/scDiagnosticsManuscript/](https://ccb-hms.github.io/scDiagnosticsManuscript/)

## Overview

We demonstrate scDiagnostics using two real-world single-cell datasets:

**1. COVID-19 PBMC scRNA-seq**

- COVID-19 PBMC scRNA-seq from Stephenson et al. (2021)
- 48,148 query cells (severe infection) + 23,201 reference cells (healthy controls)
- 10x Genomics scRNA-seq technology
- Demonstration of discovery and characterization of disease-associated CD14+ monocyte state

**2. MERFISH Mouse Colitis**

- Imaging-based spatially-resolved single-cell dataset from Cadinu et al. (2024)
- 29,040 query cells (DSS day 9) + 27,140 reference cells (healthy)
- MERFISH spatial transcriptomics technology
- Demonstration of discovery and characterization of disease-associated fibroblast state in mouse colitis model

For each dataset, we predict cell type labels using four popular annotation tools:

- [**Azimuth**](https://github.com/satijalab/azimuth) — Weighted k-NN mapping
- [**SingleR**](https://bioconductor.org/packages/SingleR) — Correlation-based assignment
- [**CellTypist**](https://www.celltypist.org/) — Machine learning classifier
- [**scVI/scArches**](https://docs.scvi-tools.org/) — Deep learning with VAE (GPU-accelerated)

## Quick Start

### Installation

```{r}
#| eval: false
source("R/covid/R Package Installation Pipeline.R")
```

Or for MERFISH:

```{r}
#| eval: false
source("R/merfish/R Package Installation Pipeline.R")
```

### Download Data

```{r}
#| eval: false
source("data/downloadData.R")
downloadData()
```

See documentation for details: [Setup & Installation](docs/setup.qmd), [Accessing Data](docs/data-access.qmd)

## Documentation

Full tutorials and analysis code available at [https://ccb-hms.github.io/scDiagnosticsManuscript/](https://ccb-hms.github.io/scDiagnosticsManuscript/):

- **[Setup & Installation](docs/setup.qmd)** — Install R and Python dependencies (GPU recommended for scVI/scArches)
- **[Accessing Data](docs/data-access.qmd)** — Download pre-processed datasets from Zenodo
- **[Cell type annotation](docs/adding-annotations.qmd)** — Apply all four annotation methods
- **[scDiagnostics Overview](docs/overview-scdiagnostics.qmd)** — Introduction to diagnostic framework
- **[COVID-19 Analysis](docs/results-covid.qmd)** — Annotation assessment in scRNA-seq
- **[MERFISH Analysis](docs/results-merfish.qmd)** — Annotation assessment in spatial transcriptomics
- **[Annotation Tool Diagnostics](docs/annotation-diagnostics.qmd)** — Complementary metrics across tools

## System Requirements

- **R:** Version 4.2 or later
- **Python:** Version 3.9 (for annotation pipelines)
- **GPU:** NVIDIA GPU with CUDA 11.8+ (required for scVI/scArches; tested on NVIDIA L40S on HMS O2 cluster)
- **Memory:** 16+ GB RAM (40+ GB for GPU annotation)
- **Disk:** ~50 GB for data + outputs

### Citation

If you use the code or data from this repository, we kindly ask that you cite our manuscript (details to be added upon publication).

```bibtex
@article{christidis2024scDiagnostics,
  author = {Christidis, A. and Ghazi, A. and Chawla, S. and Turaga, N. and Gentleman, R. and Geistlinger, L.},
  title = {scDiagnostics: systematic assessment of cell type annotation in single-cell transcriptomics data},
  year = {2024},
  note = {Submitted}
}
```

## Repository

**Code and Scripts:** [github.com/ccb-hms/scDiagnosticsManuscript](https://github.com/ccb-hms/scDiagnosticsManuscript)

**Data:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18274942.svg)](https://doi.org/10.5281/zenodo.18274942)

## Contact

For questions or feedback, please [open an issue](https://github.com/ccb-hms/scDiagnosticsManuscript/issues) on GitHub.








