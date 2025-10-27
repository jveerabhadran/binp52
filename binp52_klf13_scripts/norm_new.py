import logging
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scanpy as sc
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
from scipy.sparse import issparse, csr_matrix

# Set verbosity and figure parameters
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

# Suppress R callbacks
rcb.logger.setLevel(logging.ERROR)

# Activate R-to-Python conversion
ro.pandas2ri.activate()
anndata2ri.activate()

# Load processed AnnData object
adata = sc.read(filename="148_6_F_KO_filtered_qc_singlets.h5ad")

# Raw Counts Distribution Plot
p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False)
plt.title("Raw Counts Distribution")
plt.show()

# Shifted Logarithm Normalization
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
logging.info("Applied Shifted Logarithm Normalization")

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Shifted logarithm")
plt.show()

# Scran Normalization
ro.r('''
library(scran)
library(BiocParallel)
''')

# Preliminary clustering
adata_pp = adata.copy()
sc.pp.normalize_total(adata_pp)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="groups")

data_mat = adata_pp.X.T
if issparse(data_mat):
    if data_mat.nnz > 2**31 - 1:
        data_mat = data_mat.tocoo()
    else:
        data_mat = data_mat.tocsc()
ro.globalenv["data_mat"] = data_mat
ro.globalenv["input_groups"] = adata_pp.obs["groups"]
del adata_pp

ro.r('''
size_factors = sizeFactors(
    computeSumFactors(
        SingleCellExperiment(
            list(counts=data_mat)), 
            clusters = input_groups,
            min.mean = 0.1,
            BPPARAM = MulticoreParam()
    )
)
''')

adata.obs["size_factors"] = ro.r("size_factors")

scran = adata.X / adata.obs["size_factors"].values[:, None]
scran_log1p = np.log1p(scran)
adata.layers["scran_normalization"] = csr_matrix(scran_log1p)
logging.info("Applied Scran Normalization")

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
sns.histplot(adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("log1p with Scran estimated size factors")
plt.show()

# Analytic Pearson Residuals Normalization
analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])
logging.info("Applied Analytic Pearson Residuals Normalization")

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
analytic_sum = np.array(adata.layers["analytic_pearson_residuals"].sum(1)).flatten()
sns.histplot(analytic_sum, bins=100, kde=False, ax=axes[1])
axes[1].set_title("Analytic Pearson residuals")
plt.show()

# Save the normalized AnnData object
adata.write("148_6_F_KO_qc_norm.h5ad")
logging.info("Saved normalized AnnData object.")
