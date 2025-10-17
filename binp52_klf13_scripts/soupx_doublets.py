import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation
import anndata2ri
import rpy2.robjects as ro

# Activate conversion between R and Python objects
anndata2ri.activate()

# Import required R libraries for SoupX and scDblFinder
ro.r('''
library(SoupX)
library(Seurat)
library(scater)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
''')

# Set up Scanpy settings
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)

# Load raw count matrix from CSV and convert to AnnData format
df = pd.read_csv("Data/124_2_M_WT_edda.csv", index_col=0)
adata = sc.AnnData(df.T)  # Transpose to match expected format (cells as rows, genes as columns)

# Ensure unique variable names (genes)
adata.var_names_make_unique()

# Print the total number of cells before filtering
print(f"Number of cells BEFORE filtering: {adata.n_obs}")

#==========Preprocessing==========#
# Annotate mitochondrial, ribosomal, and hemoglobin genes
adata.var["mt"] = adata.var_names.str.startswith("mt-")
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
adata.var["hb"] = adata.var_names.str.contains("^Hb[^(p)]")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
)

# Plot distribution of total counts per cell
sns.displot(adata.obs["total_counts"], bins=100, kde=False)
plt.title("Total Counts per Cell Distribution")
plt.xlabel("Total Counts")
plt.ylabel("Frequency")
plt.show()

# Violin plot for mitochondrial percentage
sc.pl.violin(adata, "pct_counts_mt", jitter=0.4, multi_panel=False)

# Scatter plot to explore total counts vs. gene counts
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

#==========Outlier Detection==========#
# Define outlier detection function based on MAD (Median Absolute Deviation)
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        M > np.median(M) + nmads * median_abs_deviation(M)
    )
    return outlier

# Identify outliers for total counts, gene counts, and mitochondrial percentage
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 3)  # Total counts outliers (3 MADs)
    | is_outlier(adata, "log1p_n_genes_by_counts", 3)  # Gene count outliers (3 MADs)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 3)  # Top 20 gene percentage outliers (3 MADs)
)

# Filter out low-quality cells based on outliers and mitochondrial percentage (>5%)
adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 5 # Hard threshold for mitochondrial percentage
)
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

print(f"Number of cells after filtering: {adata.n_obs}")

# Scatter plot to explore total counts vs. gene counts
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
plt.show()

#==========SoupX for Ambient RNA Removal===========#
# Preprocess AnnData for SoupX: Normalize and log-transform the data
adata_pp = adata.copy()
sc.pp.normalize_total(adata_pp, target_sum=1)
sc.pp.log1p(adata_pp)

# Perform PCA and Leiden clustering for SoupX groups
sc.pp.pca(adata_pp)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="soupx_groups")
soupx_groups = adata_pp.obs["soupx_groups"]

# Compute UMAP embedding
sc.tl.umap(adata_pp)

# Visualize Leiden clusters (SoupX groups) on UMAP
sc.pl.umap(
    adata_pp,
    color="soupx_groups",
    title="Leiden Clusters (SoupX Groups)",
    legend_loc='on data',  
    frameon=False,
    show=True
)
plt.show()

del adata_pp  # Clean up memory

# Prepare data for SoupX (transpose to match R's expected format)
cells = adata.obs_names
genes = adata.var_names
data = adata.X.T  # Transpose to genes x cells format
data_tod = data.copy()  # Placeholder for raw matrix if unavailable

# Run SoupX in R using rpy2 interface
ro.globalenv['data'] = data
ro.globalenv['data_tod'] = data_tod
ro.globalenv['genes'] = genes
ro.globalenv['cells'] = cells
ro.globalenv['soupx_groups'] = soupx_groups

# Run SoupX in R using rpy2 interface with TF-IDF adjustments and diagnostic plots enabled
ro.r('''
set.seed(123)
library(SoupX)

rownames(data) <- genes
colnames(data) <- cells

data <- as(data, "sparseMatrix")
data_tod <- as(data_tod, "sparseMatrix")

# Create SoupChannel object
sc <- SoupChannel(data_tod, data, calcSoupProfile=FALSE)

# Set soup profile based on raw data
soupProf <- data.frame(row.names=rownames(data), est=rowSums(data)/sum(data), counts=rowSums(data))
sc <- setSoupProfile(sc, soupProf)

# Add clustering information from Python (soupx_groups)
sc <- setClusters(sc, soupx_groups)

# (Not needed, coz of manual rho value setting)
# Estimate contamination fraction with adjusted TF-IDF parameters and diagnostic plots enabled
# sc <- autoEstCont(sc, tfidfMin=0.8, soupQuantile=0.85, doPlot=TRUE)

# Change for each file
sc <- setContaminationFraction(sc, 0.031)
     
# out <- adjustCounts(sc, roundToInt=TRUE)
out <- adjustCounts(sc, method="soupOnly", roundToInt=TRUE)

plotMarkerDistribution(sc)  # Diagnostic plot for marker distribution

''')

out = ro.globalenv['out']

# Save corrected counts matrix in AnnData layers and overwrite .X with corrected counts
adata.layers["counts"] = adata.X  # Original counts layer (pre-SoupX)
adata.layers["soupX_counts"] = out.T  # Corrected counts layer (post-SoupX)
adata.X = adata.layers["soupX_counts"]

print(f"Total number of genes: {adata.n_vars}")

sc.pp.filter_genes(adata, min_cells=5)  # Keep cells with at least 10 genes expressed
# sc.pp.filter_cells(adata, min_counts=500)  # Keep cells with at least 500 total counts

print(f"Number of genes after cell filter: {adata.n_vars}")

# Prepare the SoupX-corrected count matrix for scDblFinder (genes x cells format)
data_mat = adata.X.T  # Transpose to match R's expected format

#==========Doublet Detection===========#
# scdblFinder
# Pass the data matrix to the R environment and run scDblFinder
ro.globalenv['data_mat'] = data_mat

ro.r('''
set.seed(123)

library(SingleCellExperiment)
library(scDblFinder)

sce <- SingleCellExperiment(list(counts=data_mat))
sce <- scDblFinder(sce)

doublet_score <- sce$scDblFinder.score
doublet_class <- sce$scDblFinder.class
''')

# Retrieve scDblFinder results back into Python and store in AnnData .obs
doublet_score = ro.globalenv['doublet_score']
doublet_class = ro.globalenv['doublet_class']
adata.obs["scDblFinder_score"] = doublet_score
adata.obs["scDblFinder_class"] = doublet_class

print(adata.obs["scDblFinder_class"].value_counts())

# Compute PCA with reduced components
sc.pp.pca(adata, n_comps=50)

# Compute neighbors using PCA-reduced data
sc.pp.neighbors(adata, use_rep='X_pca')

# Run UMAP
sc.tl.umap(adata)

# Visualize UMAP embedding colored by scDblFinder classifications
sc.pl.umap(adata, color="scDblFinder_class")
plt.title("UMAP Embedding: Doublets vs Singlets")

# Scrublet
def run_scrublet(adata, expected_doublet_rate=0.08):
    np.random.seed(123) 
    counts_matrix = adata.layers["counts"].toarray() if hasattr(adata.layers["counts"], "toarray") else adata.layers["counts"]
    scrublet_instance = scr.Scrublet(counts_matrix, expected_doublet_rate=expected_doublet_rate)  # Transpose to cells x genes format

    doublet_scores, predicted_doublets = scrublet_instance.scrub_doublets()

    # Store results in AnnData object
    adata.obs['scrublet_score'] = doublet_scores
    adata.obs['scrublet_class'] = ['doublet' if x else 'singlet' for x in predicted_doublets]

    # Plot histogram of doublet scores using Scrublet's built-in method
    scrublet_instance.plot_histogram()

    print(f"Number of doublets detected by Scrublet: {sum(predicted_doublets)}")
    return adata

# Run Scrublet on the dataset with SoupX-corrected counts:
run_scrublet(adata)

#===========Save the Files for further processing============# 
# Save AnnData object with doublets included (both scDblFinder and Scrublet results)
adata.write("quality_control_with_doublets.h5ad")

# Remove doublets based on scDblFinder classification ("singlets" only)
adata_filtered_scDblFinder = adata[adata.obs["scDblFinder_class"] == "singlet"].copy()

# Save filtered dataset without doublets (based on scDblFinder classification)
adata_filtered_scDblFinder.write("quality_control_filtered_scDblFinder.h5ad")

print("Datasets saved with and without doublets.")