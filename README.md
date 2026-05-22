# Single-nucleus RNA-seq Analysis Pipeline
## BINP52 - Master's Thesis

> A fully reproducible pipeline for single-nucleus RNA-seq quality control and downstream analysis of _Klf13_ knockout in the mouse brain, using Python and R.

> Author: Jyothishree Veerabhadran

> Programme: Bioinformatics MSc

> Institution: Lund University

> Year: 2026

---
**Table of Contents**
1. Project Overview
2. Repository Structure
3. Pipeline Steps
4. Setup & Installation
5. Usage
6. Data
7. Reproducibility
8. Troubleshooting
9. Acknowledgements
10. Citation
11. Contact
---

**Project Overview**

This repository contains a batch-processing pipeline for single-nucleus RNA sequencing (snRNA-seq) data, developed as part of the BINP52 Master's thesis. The pipeline covers all major steps from raw count matrix ingestion to cell type annotation and subcluster analysis, with a focus on:

- Ambient RNA correction using SoupX
- Doublet detection using scDblFinder and Scrublet
- Normalisation and feature selection using scran and Scanpy
- Dimensionality reduction, clustering, and annotation
- GABAergic neuron subclustering
  
The pipeline is designed for multi-sample batch processing with robust error handling and clear separation of each analysis step.

---

**Repository Structure**
```
binp52/
│
├── notebooks/                    # Jupyter notebooks (interactive analysis)
│   ├── 01_quality_filtering.ipynb
│   ├── 02_normalization.ipynb
│   ├── 03_concatenation.ipynb
│   ├── 04_feature_selection.ipynb
│   ├── 05_dimensionality_reduction.ipynb
│   ├── 06_clustering.ipynb
│   ├── 07_annotation.ipynb
│   └── 08_gaba_recluster.ipynb
│
│
├── Data/                         # Raw count matrices (.csv)
│
├── requirements.txt              # Python dependencies
├── install.R                     # R dependencies
├── .gitignore
├── LICENSE
├── CITATION.cff
└── README.md
```
---
**Pipeline Steps**

| Step | Notebook / Script | Description | Output |
|------|-------------------|-------------|--------|
| 1 | `quality_filtering` | Mitochondrial/ribosomal/Hb gene annotation, outlier detection, SoupX ambient RNA correction, doublet detection (scDblFinder + Scrublet) | Per-sample QC plots, filtered `.h5ad` |
| 2 | `normalization` | Scran-based normalisation per sample | Normalised `.h5ad` per sample |
| 3 | `concatenation` | Merge all per-sample objects into a single AnnData | Combined `.h5ad` |
| 4 | `feature_selection` | Highly variable gene selection | Feature-selected `.h5ad` |
| 5 | `dimensionality_reduction` | PCA and UMAP | Reduced-dimension `.h5ad` + UMAP plots |
| 6 | `clustering` | Leiden clustering on combined dataset | Cluster assignments + plots |
| 7 | `annotation` | Cell type annotation of clusters | Annotated `.h5ad` + summary tables |
| 8 | `gaba_recluster` | Subclustering of GABAergic neuron population | Subclustered `.h5ad` + plots |

---
**Setup & Installation**
**Requirements**

Python 3.10+
R 4.4.3
Conda (recommended) or pip

---
Python Setup
Option A — pip (simple):
```bash
pip install -r requirements.txt
```
Option B — Conda (recommended for full reproducibility):
```bash
conda create -n binp52 python=3.10
conda activate binp52
pip install -r requirements.txt
```
---
R Setup
Open R (version 4.4.3) and run:
```r
source("install.R")
```
This installs all required Bioconductor and CRAN packages.
> **Note:** `rpy2` and `anndata2ri` are Python packages that bridge Python and R. Make sure R is installed and accessible in your system PATH before installing them.
---

**Usage**

Run as Jupyter Notebooks (recommended)
```bash
conda activate binp52
jupyter notebook
```
Open the notebooks in order (`01_` → `08_`) from the `notebooks/` folder.

> Scripts must be run **in order**. Each step reads the output of the previous step.
---

**Data**

A sample raw data (`.csv` count matrices) is included in this repository.

Input format: Each sample is a raw count matrix (genes × cells) in `.csv` format, placed in `Data/`

This is how a raw count matrix typically looks like. Each row represents a gene, and each column represents a cell. The values are raw, unnormalized read counts.

| **Gene**  | **Cell1** | **Cell2** | **Cell3** | **Cell4** |
| --------- | --------- | --------- | --------- | --------- |
| **GeneA** | 12        | 0         | 0         | 10        |
| **GeneB** | 5         | 3         | 0         | 8         |
| **GeneC** | 0         | 0         | 0         | 2         |
| **GeneD** | 9         | 7         | 0         | 0         |

---
**Reproducibility**

- All random seeds are fixed in both Python and R
- All package versions are pinned in `requirements.txt` and `install.R`
- Input/output paths are fully parameterised within each script
- The pipeline skips samples with errors and logs them, ensuring batch runs complete

---
**Troubleshooting**

| Problem | Solution |
|---------|----------|
| Missing R packages | Run `source("install.R")` in R 4.4.3 |
| `rpy2` import error | Ensure R is in your system PATH and matches version 4.4.3 |
| Version mismatch warnings | Use a fresh conda environment with pinned versions from `requirements.txt` |
| Memory errors on large samples | Reduce batch size or increase available RAM; consider running per-sample steps on HPC |
| `.h5ad` file not found | Ensure previous pipeline step completed successfully |

---
**Acknowledgements**

 - [snRNA-seq Analysis Pipeline](https://www.sc-best-practices.org/)
 - [Cacoa Analysis](https://kkh.bric.ku.dk/xian/pipeline/CPHscPipe/cacoa.html)

---
**Citation**

If you use this pipeline or code in your work, please cite:
```
Veerabhadran, J. (2026). Single-nucleus RNA-seq Analysis Pipeline for KLF13 (BINP52).
GitHub. https://github.com/jveerabhadran/binp52
```
Or use the `CITATION.cff` file in this repository (GitHub renders a "Cite this repository" button automatically).

---
**License**

This project is licensed under the MIT License, see the LICENSE file for details.

---
**Contact**

For questions or issues, please open an issue on GitHub or contact:
Jyothishree Veerabhadran  
(itsmejyoee@gmail.com)
