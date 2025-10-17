import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from anndata2ri import py2rpy, rpy2py
import matplotlib.pyplot as plt
import pandas as pd
from rpy2.robjects import pandas2ri

# Setup scanpy plotting and verbosity
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)
sns.set(style="whitegrid")

# Enable R <-> Python conversion
anndata2ri.activate()
pandas2ri.activate()

# 1) Load your CSV data: genes as rows, cells as columns
input_csv_file = "Data/96_6_M_KO_edda.csv" 
df = sc.read_csv(input_csv_file).T  # Transpose so cells x genes

adata = sc.AnnData(df)
adata.var_names_make_unique()

# 2) Annotate gene sets for QC 
adata.var["mt"] = adata.var_names.str.startswith("mt-")
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
adata.var["hb"] = adata.var_names.str.contains("^Hb[^(p)]")

# 3) Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

# Print Summary Statistics Before Filtering
print("Summary Statistics BEFORE Filtering:")
print(adata.obs[['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb', 'pct_counts_in_top_20_genes']].describe())

# 4) Visualize QC metrics Before Filtering (side by side)
fig, axs = plt.subplots(2, 3, figsize=(18, 10))

sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axs[0,0])
axs[0,0].set_title("Total Counts Before Filtering")
axs[0,0].set_xlim(0,70000)

sns.violinplot(y=adata.obs['pct_counts_mt'], ax=axs[0,1])
axs[0,1].set_title("Mitochondrial % Before Filtering")

sns.scatterplot(x=adata.obs["total_counts"], y=adata.obs["n_genes_by_counts"], hue=adata.obs["pct_counts_mt"],
                palette="magma", ax=axs[0,2])
axs[0,2].set_title("Counts vs Genes (Colored by Mito %) Before Filtering")

sns.violinplot(y=adata.obs['pct_counts_ribo'], ax=axs[1,0])
axs[1,0].set_title("Ribosomal % Before Filtering")

sns.violinplot(y=adata.obs['pct_counts_hb'], ax=axs[1,1])
axs[1,1].set_title("Hemoglobin % Before Filtering")

sns.violinplot(y=adata.obs['pct_counts_in_top_20_genes'], ax=axs[1,2])
axs[1,2].set_title("Top 20 Genes % Before Filtering")

plt.tight_layout()
plt.show()

# 5) MAD-based outlier detection function
def is_outlier(adata, metric:str, nmads:int):
    M = adata.obs[metric]
    return (M < np.median(M) - nmads * median_abs_deviation(M)) | (M > np.median(M) + nmads * median_abs_deviation(M))

# 6) Mark outliers
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > 0.05

print(f"Cells flagged as outliers: {adata.obs['outlier'].sum()}")
print(f"Cells flagged as mitochondrial outliers: {adata.obs['mt_outlier'].sum()}")

# 7) Filter low-quality cells
print(f"Cells before filtering: {adata.n_obs}")
adata_filtered = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
print(f"Cells after filtering: {adata_filtered.n_obs}")

# Print Summary Statistics After Filtering
print("\nSummary Statistics AFTER Filtering:")
print(adata_filtered.obs[['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb', 'pct_counts_in_top_20_genes']].describe())

# 8) Visualize QC metrics After Filtering (side by side)
fig, axs = plt.subplots(2, 3, figsize=(18, 10))

sns.histplot(adata_filtered.obs["total_counts"], bins=100, kde=False, ax=axs[0,0])
axs[0,0].set_title("Total Counts After Filtering")
axs[0,0].set_xlim(0, 70000)

sns.violinplot(y=adata_filtered.obs['pct_counts_mt'], ax=axs[0,1])
axs[0,1].set_title("Mitochondrial % After Filtering")

sns.scatterplot(x=adata_filtered.obs["total_counts"], y=adata_filtered.obs["n_genes_by_counts"], hue=adata_filtered.obs["pct_counts_mt"],
                palette="magma", ax=axs[0,2])
axs[0,2].set_title("Counts vs Genes (Colored by Mito %) After Filtering")

sns.violinplot(y=adata_filtered.obs['pct_counts_ribo'], ax=axs[1,0])
axs[1,0].set_title("Ribosomal % After Filtering")

sns.violinplot(y=adata_filtered.obs['pct_counts_hb'], ax=axs[1,1])
axs[1,1].set_title("Hemoglobin % After Filtering")

sns.violinplot(y=adata_filtered.obs['pct_counts_in_top_20_genes'], ax=axs[1,2])
axs[1,2].set_title("Top 20 Genes % After Filtering")

plt.tight_layout()
plt.show()

# 9) Run clustering for SoupX ambient RNA correction
# Run SoupX on raw counts before normalization
adata_filtered.layers["raw_counts"] = adata_filtered.X.copy()
# Create a copy to preserve original data
adata_pp = adata_filtered.copy()
# Normalize counts per cell to correct sequencing depth
sc.pp.normalize_per_cell(adata_pp)
# Log-transform the normalized counts
sc.pp.log1p(adata_pp)
# PCA for dimensionality reduction
sc.pp.pca(adata_pp)
# Compute neighborhood graph on PCA representation
sc.pp.neighbors(adata_pp)
# Run Leiden clustering, save cluster labels as 'soupx_groups'
sc.tl.leiden(adata_pp, key_added='soupx_groups')
# Compute UMAP embedding from neighborhood graph
sc.tl.umap(adata_pp)
# Now you can access cluster labels and plot UMAP colored by them
soupx_groups = adata_pp.obs["soupx_groups"]
# Plot UMAP embedding colored by the clusters found
sc.pl.umap(adata_pp, color="soupx_groups")
# Optional: clean up
del adata_pp

# Prepare count matrix for SoupX (raw counts)
data_tod = adata_filtered.layers["raw_counts"].T
cells = adata_filtered.obs_names
genes = adata_filtered.var_names
data = adata_filtered.layers["raw_counts"].T

ro.globalenv['data'] = data
ro.globalenv['data_tod'] = data_tod
ro.globalenv['genes'] = genes
ro.globalenv['cells'] = cells
ro.globalenv['soupx_groups'] = soupx_groups

# 10) Run SoupX correction in R with automatic contamination estimation
ro.r('''
set.seed(123)
library(SoupX)
library(Matrix)

rownames(data) <- genes
colnames(data) <- cells
data <- as(data, "sparseMatrix")
data_tod <- as(data_tod, "sparseMatrix")

sc <- SoupChannel(data_tod, data, calcSoupProfile = FALSE)

soupProf <- data.frame(row.names=rownames(data), est=rowSums(data)/sum(data), counts=rowSums(data))
sc <- setSoupProfile(sc, soupProf)
sc <- setClusters(sc, soupx_groups)
sc <- autoEstCont(sc, contaminationRange=c(0.01, 0.1), doPlot=TRUE)

cat("Estimated contamination fraction:", sc$conEstimate, "\\n")

out <- adjustCounts(sc, method="soupOnly", roundToInt=TRUE)
out <- as.matrix(out)
''')

out = ro.globalenv['out']

adata_filtered.layers["raw_counts"] = adata_filtered.X.copy()
adata_filtered.layers["soupX_counts"] = out.T

# Convert 'out' (SoupX output) matrix to DataFrame with correct gene names
out_df = pd.DataFrame(out.T, index=adata_filtered.obs_names, columns=adata_filtered.var_names)

# No need to subset columns here if you already used var_names of filtered object

# Assign back to layered matrix
adata_filtered.layers["soupX_counts"] = out_df.values

# Update X if you want to use corrected counts as main expression matrix
adata_filtered.X = adata_filtered.layers["soupX_counts"]

# adata_filtered is your AnnData object after QC and SoupX
# 'raw_counts' layer: before SoupX
# 'soupX_counts' layer: after SoupX

gene = "Klf13"

# Pre-SoupX
if 'raw_counts' in adata_filtered.layers and gene in adata_filtered.var_names:
    raw_idx = list(adata_filtered.var_names).index(gene)
    raw_klf13 = adata_filtered.layers["raw_counts"][:, raw_idx]
else:
    raw_klf13 = None

# Post-SoupX
if 'soupX_counts' in adata_filtered.layers and gene in adata_filtered.var_names:
    soupx_idx = list(adata_filtered.var_names).index(gene)
    soupx_klf13 = adata_filtered.layers["soupX_counts"][:, soupx_idx]
else:
    soupx_klf13 = None

# Plot comparison
if (raw_klf13 is not None) and (soupx_klf13 is not None):
    plt.figure(figsize=(10,4))
    plt.subplot(1,2,1)
    plt.hist(raw_klf13, bins=50, color="blue")
    plt.title("Klf13 Expression (Raw Counts)")
    plt.xlabel("Raw Counts")
    plt.ylabel("Cells")
    plt.subplot(1,2,2)
    plt.hist(soupx_klf13, bins=50, color="green")
    plt.title("Klf13 Expression (SoupX Corrected)")
    plt.xlabel("SoupX Counts")
    plt.ylabel("Cells")
    plt.tight_layout()
    plt.show()
else:
    print("Could not find Klf13 in layers or AnnData.")

# 11) Filter genes expressed in at least 10 cells on filtered data
print(f"Genes before filtering: {adata_filtered.n_vars}")
sc.pp.filter_genes(adata_filtered, min_cells=10)
print(f"Genes after filtering: {adata_filtered.n_vars}")


# 12) Run scDblFinder (R) doublet detection
data_mat = adata_filtered.X.T
ro.globalenv['data_mat'] = data_mat

ro.r('''
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)

set.seed(123)
sce <- SingleCellExperiment(list(counts=data_mat))
sce <- scDblFinder(sce)

doublet_score <- sce$scDblFinder.score
doublet_class <- sce$scDblFinder.class
''')

adata_filtered.obs["scDblFinder_score"] = np.array(ro.globalenv['doublet_score'])
adata_filtered.obs["scDblFinder_class"] = [str(x) for x in np.array(ro.globalenv['doublet_class'])]

print("scDblFinder doublet call counts:")
print(adata_filtered.obs.scDblFinder_class.value_counts())

print(adata_filtered.obs["scDblFinder_class"].value_counts())

# 13) Embedding for visualization
sc.pp.pca(adata_filtered, n_comps=50)
sc.pp.neighbors(adata_filtered, use_rep="X_pca")
sc.tl.umap(adata_filtered)
# Plot doublet calls on UMAP
sc.pl.umap(adata_filtered, color="scDblFinder_class", title="UMAP: Doublets vs Singlets", frameon=False)

# Save your final filtered and annotated AnnData object
adata_filtered.write("96_6_M_KO_filtered_qc_doublets.h5ad")

adata_singlets = adata_filtered[adata_filtered.obs["scDblFinder_class"] == "singlet"].copy()
adata_singlets.write("96_6_M_KO_filtered_qc_singlets.h5ad")

print("Files written and done.")




