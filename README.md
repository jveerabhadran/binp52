# Single-nucleus RNA-seq Analysis Pipeline

## BINP52: Masters Thesis

**Note: The Repository is still under construction and the README is not fully written or formatted.

This repository contains a fully reproducible, pipeline for single-nucleus RNA-seq quality control and downstream analysis using Python and R.  The pipeline is designed for batch-processing of multiple samples, robust error handling, and clear separation of each analysis step.

---
## Author

- [@Jyothishree Veerabhadran](https://github.com/jveerabhadran)
---

## Input Data
This is how a raw count matrix typically looks like. Each row represents a gene, and each column represents a cell. The values are raw, unnormalized read counts.

| **Gene**  | **Cell1** | **Cell2** | **Cell3** | **Cell4** |
| --------- | --------- | --------- | --------- | --------- |
| **GeneA** | 12        | 0         | 0         | 10        |
| **GeneB** | 5         | 3         | 0         | 8         |
| **GeneC** | 0         | 0         | 0         | 2         |
| **GeneD** | 9         | 7         | 0         | 0         |


## Processing Steps
The pipeline is organized into the following steps:

1. **Quality Filtering**  
   *Script:* `quality_filtering.py`  
   *Description:* Performs QC and filtering (mitochondrial/ribosomal/hb gene annotation, outlier detection, SoupX correction, doublet detection) on each raw `.csv` in `Data/` **separately**.  
   *Outputs:* Per-sample QC plots and filtered `.h5ad` files.

2. **Normalization**  
   *Script:* `normalization.py`  
   *Description:* Applies normalization to each filtered sample separately.  
   *Outputs:* Normalized `.h5ad` files per sample.

3. **Concatenation**  
   *Script:* `concatenate.py`  
   *Description:* Merges all per-sample `.h5ad` files into a single AnnData object for joint analysis.  
   *Outputs:* Combined `.h5ad` file.

4. **Feature Selection**  
   *Script:* `feature_selection.py`  
   *Description:* Selects highly variable genes for each sample.  
   *Outputs:* Feature-selected `.h5ad` for combined file.

5. **Dimensionality Reduction**  
   *Script:* `dimensionality_reduction.py`  
   *Description:* Performs PCA and UMAP on each sample.  
   *Outputs:* Reduced-dimension `.h5ad` files and plots for combined file.

6. **Clustering**  
   *Script:* `clustering.py`  
   *Description:* Performs clustering on the combined dataset.  
   *Outputs:* Cluster assignments and plots.

7. **Annotation**  
   *Script:* `annotation.py`  
   *Description:* Annotates clusters/cell types.  
   *Outputs:* Annotated `.h5ad` file and summary tables.

8. **Reclustering**  
   *Script:* `gaba_recluster.py`  
   *Description:* Reclusters GABAergic cell population.  
   *Outputs:* Subclustered data and plots.

---

## Directory Structure
## Setup
## Usage
## Reproducibility

- All random seeds are set in both Python and R.
- All package versions are pinned.
- Input/output is fully parameterized.
- The pipeline is robust to errors and logs issues per sample.

---

## Troubleshooting

- **Missing R packages:** Make sure all R dependencies are installed as above.
- **Version mismatch:** Check your Python and R versions match the requirements files.
- **Other errors:** Review the console output for error messages. The pipeline will skip files with errors and continue processing others.

---

    
## Acknowledgements

 - [snRNA-seq Analysis Pipeline](https://www.sc-best-practices.org/)
 - [Cacoa Analysis](https://kkh.bric.ku.dk/xian/pipeline/CPHscPipe/cacoa.html)

---
   ---
## Contact

For questions or issues, please open an issue on GitHub or contact the owner.
