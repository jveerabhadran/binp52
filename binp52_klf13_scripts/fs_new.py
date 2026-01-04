import anndata2ri
import matplotlib.pyplot as plt
import numpy as np
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scanpy as sc
import seaborn as sns
from rpy2.robjects.conversion import localconverter
import warnings
import logging
import copy
import scipy.sparse
from memory_profiler import profile


# Suppress all warnings (for cleaner output)
warnings.filterwarnings('ignore')


# Reproducibility and Logging
SEED = 123
np.random.seed(SEED)
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)
rcb.logger.setLevel(logging.ERROR)


# Directly set file paths and parameters
input_file = "combined_8_samples_batch_corrected.h5ad"
output_file = "feature_selected_scran.h5ad"


# Load dataset
adata = sc.read_h5ad(input_file)


# Integrate R scry package for feature selection
anndata2ri.activate()
with localconverter(anndata2ri.converter):
     ro.globalenv["adata"] = adata
     ro.r('''
         library(scry)
        sce = devianceFeatureSelection(adata, assay="X")
        binomial_deviance = rowData(sce)$binomial_deviance
    ''')
binomial_deviance = np.array(ro.r("binomial_deviance"))

# Select top N genes
idx = binomial_deviance.argsort()[-4000:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True


# Annotate AnnData
adata.var["highly_deviant"] = mask
adata.var["binomial_deviance"] = binomial_deviance


# Identify highly variable genes using Scanpy
sc.pp.highly_variable_genes(adata, layer="scran_normalization")


# Plot feature selection
plt.figure(figsize=(8, 6))
ax = sns.scatterplot(
     data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5,
     palette={True: "red", False: "grey"},
     legend=False,
)
plt.title("Feature Selection (Red: Highly Deviant)")
plt.show()


# Save the annotated AnnData
adata.write(output_file)


print("Feature selection complete.")
print(f"Saved to: {output_file}")