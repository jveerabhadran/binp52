import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import silhouette_score
import pandas as pd
import anndata


# Load data
adata = anndata.read_h5ad("clustered_scran_annotated_ko.h5ad")


# Select GABAergic cells with expression > 0 in any marker gene
#gaba_marker_genes = ['Gad2', 'Gad1', 'Dlx1', 'Arx', 'Slc32a1', 'Dlx5', 'Dlx6', 'Lhx6', 'Lhx8']
gaba_marker_genes = ['Gad2', 'Gad1', 'Dlx1', 'Dlx2']
gaba_marker_genes = [g for g in gaba_marker_genes if g in adata.var_names]


X_gaba = adata[:, gaba_marker_genes].X
if hasattr(X_gaba, "A1"):
    mask = (X_gaba > 0).sum(axis=1).A1 > 0
else:
    mask = (X_gaba > 0).sum(axis=1) > 0


adata_gaba = adata[mask].copy()
print(f"Number of GABAergic cells: {adata_gaba.n_obs}")

adata_gaba = adata_gaba.copy()

# USE EXISTING LAYERS (no new normalization!)
#adata_gaba.X = adata_gaba.layers["scran_normalization"].copy()
#adata_gaba.X = np.nan_to_num(adata_gaba.X, nan=0.0, posinf=10.0, neginf=0.0)
#print(f"GABA subset: {adata_gaba.n_obs} cells, {adata_gaba.n_vars} genes")

# Scale + PCA (safe for sparse)
sc.pp.scale(adata_gaba, max_value=10, zero_center=False)
sc.pp.pca(adata_gaba, n_comps=50, svd_solver='arpack')
sc.pp.neighbors(adata_gaba, use_rep='X_pca', n_neighbors=15)

# Normalize and log transform
adata_gaba.layers['counts'] = adata_gaba.X.copy()
sc.pp.normalize_total(adata_gaba, target_sum=1e4)
sc.pp.log1p(adata_gaba)
adata_gaba.layers['log1p_norm_gaba'] = adata_gaba.X.copy()

adata_gaba.X = np.nan_to_num(adata_gaba.X, nan=0, posinf=10, neginf=-10)
print("NaNs cleaned safely!")

# GABA reclustering 
print(f"GABA subset: {adata_gaba.n_obs} cells, {adata_gaba.n_vars} genes")
print("Using GABA-specific processing")

#sc.pp.scale(adata_gaba, max_value=10)
#sc.pp.pca(adata_gaba, n_comps=50)
#sc.pp.neighbors(adata_gaba, use_rep='X_pca', n_neighbors=15)


sc.tl.leiden(adata_gaba, resolution=0.6, key_added='leiden_res0_6')
sc.tl.leiden(adata_gaba, resolution=0.7, key_added='leiden_res0_7')
sc.tl.leiden(adata_gaba, resolution=0.8, key_added='leiden_res0_8')
sc.tl.leiden(adata_gaba, resolution=0.9, key_added='leiden_res0_9')
sc.tl.leiden(adata_gaba, resolution=1.5, key_added='leiden_res1_5')


# Leiden clustering at multiple resolutions with silhouette scoring
resolutions = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0]
sil_scores = {}
leiden_keys = [f'leiden_{res}' for res in resolutions]


for res, key in zip(resolutions, leiden_keys):
    sc.tl.leiden(adata_gaba, resolution=res, key_added=key)
    labels = adata_gaba.obs[key].astype(int)
    sil_scores[res] = silhouette_score(adata_gaba.obsm['X_pca'], labels)
    print(f"Resolution {res} - Silhouette score: {sil_scores[res]:.3f}")


# Compute UMAP embedding
sc.tl.umap(adata_gaba, min_dist=0.3, spread=1.0)


# Step 2: Cell Cycle Marker Analysis
# Define marker genes for S phase and G2/M phase
s_genes = [
     "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1",
    "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2",
    "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1",
    "CHAF1B", "BRIP1", "E2F8"
]


g2m_genes = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
    "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2",
    "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1",
    "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2",
    "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3",
    "GAS2L3", "CBX5", "CENPA"
]


# Step 2: Cell Cycle Marker Analysis
s_genes_title = [gene.title() for gene in s_genes]
g2m_genes_title = [gene.title() for gene in g2m_genes]


s_genes_filtered = [gene for gene in s_genes_title if gene in adata_gaba.var_names]
g2m_genes_filtered = [gene for gene in g2m_genes_title if gene in adata_gaba.var_names]


print("Filtered S genes:", s_genes_filtered)
print("Filtered G2M genes:", g2m_genes_filtered)


if s_genes_filtered and g2m_genes_filtered:
    sc.tl.score_genes_cell_cycle(adata_gaba, s_genes=s_genes_filtered, g2m_genes=g2m_genes_filtered)
else:
    print("No valid S or G2M genes found for cell cycle scoring.")


sc.pl.umap(adata_gaba, color=["S_score", "G2M_score"], cmap="viridis", title="UMAP: Cell Cycle Scores")


adata_gaba.obs["cell_cycle_phase"] = np.where(
    adata_gaba.obs["S_score"] > adata_gaba.obs["G2M_score"], "S",
    np.where(adata_gaba.obs["G2M_score"] > adata_gaba.obs["S_score"], "G2M", "G1")
)


# Tabulate cell cycle phase assignment by cluster
cluster_phase_counts = adata_gaba.obs.groupby(['leiden_res0_7', 'cell_cycle_phase']).size().unstack(fill_value=0)


# Visualize UMAP, colored by phase
sc.pl.umap(adata_gaba, color="cell_cycle_phase", palette="Set2")


# Define lineage marker sets
marker_sets = {
    "MGE Progenitors": ["Nkx2-1", "Sp9", "Prdm16", "Lhx6"],
    "MGE Neuroblasts": ["Ackr3", "St18", "Maf", "Erbb4", "Sox6", "Foxp1"],
    "CGE Progenitors": ["Nr2f2", "Pax6", "Nr2f1", "Sp8", "Ascl1"],
    "CGE Neuroblasts": ["Prox1", "Calb2", "Lamp5", "Synpr"],
    "LGE Progenitors": ["Isl1", "Ebf1"],
    "LGE Neuroblasts": ["Tac1", "Ddah1"]
}


# Create a unique set of all markers that exist in your var_names (removes duplications)
# After HVG filtering:
# Ensure all_markers only contains genes still present in filtered adata_gaba.var_names 
all_markers = []
for genes in marker_sets.values():
    for g in genes:
        if g in adata_gaba.var_names:
            all_markers.append(g)


# Remove duplicates but preserve order (optional)
all_markers = list(dict.fromkeys(all_markers))


# Similarly, for lineage_markers:
selected_sets = ["CGE Progenitors", "MGE Progenitors", "LGE Progenitors", "CGE Neuroblasts", "MGE Neuroblasts", "LGE Neuroblasts"]


lineage_markers = []
for s in selected_sets:
    for g in marker_sets[s]:
        if g in adata_gaba.var_names:
            lineage_markers.append(g)


lineage_markers = list(dict.fromkeys(lineage_markers))


for lineage_name, markers in marker_sets.items():
    gene_list = [g for g in markers if g in adata_gaba.var_names]
    if gene_list:
        sc.tl.score_genes(adata_gaba, gene_list=gene_list, score_name=f"{lineage_name}_score")


sc.pl.umap(adata_gaba, color=[f"{ln}_score" for ln in marker_sets.keys()])


# 1. UMAP with only lineage marker set names as cluster labels (e.g. "MGE Neuroblast")


key_for_umap = leiden_keys[2] # using resolution 0.8


labels = adata_gaba.obs[key_for_umap].astype(str).values
groups = sorted(set(labels))


# Filter marker genes to those present in your dataset
marker_list = [g for g in all_markers if g in adata_gaba.var_names]


# Create a temporary AnnData object with all cells but only marker genes
adata_marker_subset = adata_gaba[:, marker_list].copy()


umap = adata_gaba.obsm['X_umap']


lineage_scores = {}
for group in groups:
    cells = labels == group  # boolean mask for all cells in adata_gaba
    scores = {}
    for lineage_name, marker_genes in marker_sets.items():
        marker_genes_in = [g for g in marker_genes if g in adata_marker_subset.var_names]
        if marker_genes_in:
            expr = adata_marker_subset[cells, marker_genes_in].X
            if hasattr(expr, "toarray"):
                expr = expr.toarray()
            marker_mean = np.mean(expr)
        else:
            marker_mean = 0
        scores[lineage_name] = marker_mean
    best_lineage = max(scores, key=scores.get)
    lineage_scores[group] = best_lineage.replace("Progenitors", "_prog").replace("Neuroblast", "_neuro")

lineage_annotations = {}
for cluster in groups:
    scores = {ln: adata_gaba.obs.loc[adata_gaba.obs[key_for_umap] == cluster, f"{ln}_score"].mean()
              for ln in marker_sets.keys()}
    best_lineage = max(scores, key=scores.get)
    lineage_annotations[cluster] = best_lineage
adata_gaba.obs['lineage_annotation'] = adata_gaba.obs[key_for_umap].map(lineage_annotations)


sc.pl.umap(adata_gaba, color='lineage_annotation')



for gene in all_markers:
    if gene in adata_gaba.var_names:
        sc.pl.umap(adata_gaba, color=gene, cmap='viridis', title=f"{gene} expression", show=True)



# 2. Differential Ranking and Plotting Top Markers per Cluster
# 'marker_list' is your custom list of genes (all_markers from earlier)
marker_list = [g for g in all_markers if g in adata_gaba.var_names]


# Create a temporary AnnData with only these marker genes
adata_marker_subset = adata[:, marker_list].copy()


# Run differential ranking only on marker genes subset
sc.tl.rank_genes_groups(
    adata_marker_subset, # Use subset with marker genes only
    groupby=key_for_umap, # clustering key, e.g., 'leiden_0.8'
    method='wilcoxon',
    n_genes=20
)


# Extract clusters/groups from the subset
groups = adata_marker_subset.obs[key_for_umap].cat.categories if hasattr(adata_marker_subset.obs[key_for_umap], 'cat') else adata_marker_subset.obs[key_for_umap].unique()


unique_top_genes = {}
used_genes = set()


# Find one unique top gene per cluster from DE results
for group in groups:
    gene_found = False
    for gene in adata_marker_subset.uns['rank_genes_groups']['names'][group]:
        if gene not in used_genes:
            unique_top_genes[group] = gene
            used_genes.add(gene)
            gene_found = True
            break
    if not gene_found:
        unique_top_genes[group] = adata_marker_subset.uns['rank_genes_groups']['names'][group][0]

print("Top marker genes per cluster:", unique_top_genes)



sc.tl.umap(
    adata_gaba,
    min_dist=0.2, # more separation; try 0.1–0.4
    spread=1.2,
)


# Plot UMAP with cluster colors (from full adata to preserve embeddings and clusters)
sc.pl.umap(adata_gaba, color=key_for_umap, size=30, legend_loc='on data')


# Plot UMAP again manually to add labels at cluster centroids
labels = adata_gaba.obs[key_for_umap].astype(str).values
umap = adata_gaba.obsm['X_umap']
groups = sorted(set(labels))
palette = sns.color_palette('tab20', n_colors=len(groups))


plt.figure(figsize=(8, 7))
for i, group in enumerate(groups):
    idx = labels == group
    plt.scatter(umap[idx, 0], umap[idx, 1], s=10, color=palette[i], label=f'Cluster {group}', alpha=0.7)
    centroid = umap[idx].mean(axis=0)
    label_text = f"{group}\n{unique_top_genes.get(group, '')}"
    plt.text(centroid[0], centroid[1], label_text, fontsize=10,
             ha='center', va='center', bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'))   

plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.title('UMAP with Cluster Labels and Top Markers')
plt.tight_layout()
plt.legend(loc='best', fontsize=8)
plt.show()



import matplotlib.pyplot as plt
import seaborn as sns


# Your marker sets
marker_sets = {
    "MGE Progenitors": ["Nkx2-1", "Sp9", "Prdm16", "Lhx6"],
    "MGE Neuroblasts": ["Ackr3", "St18", "Maf", "Erbb4", "Sox6", "Foxp1"],
    "CGE Progenitors": ["Nr2f2", "Pax6", "Sp8", "Ascl1"],
    "CGE Neuroblasts": ["Prox1", "Synpr"],
    "LGE Progenitors": ["Isl1", "Ebf1"],
    "LGE Neuroblasts": ["Tac1", "Ddah1"]
}


# Map cluster names to group short codes
cluster_group_codes = {
    "MGE Progenitors": "MP",
    "MGE Neuroblasts": "MN",
    "CGE Progenitors": "CP",
    "CGE Neuroblasts": "CN",
    "LGE Progenitors": "LP",
    "LGE Neuroblasts": "LN"
}


labels = adata_gaba.obs[key_for_umap].astype(str).values  # cluster numbers as strings
umap = adata_gaba.obsm['X_umap']
groups = sorted(set(labels), key=int)
palette = sns.color_palette('tab20', n_colors=len(groups))


# Reverse dictionary: markers to group code, to find group code by marker gene presence
marker_to_group_code = {}
for group_name, markers in marker_sets.items():
    gc = cluster_group_codes.get(group_name, "")
    for m in markers:
        marker_to_group_code[m] = gc


plt.figure(figsize=(8, 7))


for i, group in enumerate(groups):
    idx = labels == group
    plt.scatter(umap[idx, 0], umap[idx, 1], s=10, color=palette[i], label=f'Cluster {group}', alpha=0.7)

    centroid = umap[idx].mean(axis=0)

    # Suppose you have the marker gene expression average per cluster in this example:
    # Calculate mean expression of each marker gene for this cluster
    mean_exp = adata_gaba[idx].to_df().mean()

    # Find which group code has the highest expressed marker genes in this cluster
    group_code_scores = {}
    for gene, expr in mean_exp.items():
        gc = marker_to_group_code.get(gene, None)
        if gc:
            group_code_scores[gc] = group_code_scores.get(gc, 0) + expr

    # Determine dominant group code by max score
    if group_code_scores:
        dominant_gc = max(group_code_scores, key=group_code_scores.get)
    else:
        dominant_gc = ""

    # Select markers from dominant group
    group_name = None
    for name, code in cluster_group_codes.items():
        if code == dominant_gc:
            group_name = name
            break
    markers = marker_sets.get(group_name, []) if group_name else []

    marker_text = ', '.join(markers[:2])
    label_text = f"{group} ({dominant_gc})\n{marker_text}"
    
    plt.text(centroid[0], centroid[1], label_text, fontsize=10,
            ha='center', va='center', bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'))

plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.title('UMAP with Cluster Number and Group Code Based on Marker Expression')
plt.tight_layout()
plt.legend(loc='best', fontsize=8)
plt.show()


chosen_key = key_for_umap  # Assuming 'key_for_umap' holds the clustering key you want (e.g., leiden_0_7)


import matplotlib.pyplot as plt
import pandas as pd


cluster_key = key_for_umap # your clustering key, e.g., 'leiden_0.6'


# Create a contingency table (counts of cell cycle phases per cluster)
phase_cluster_counts = pd.crosstab(adata_gaba.obs[cluster_key], adata_gaba.obs['cell_cycle_phase'])


# Plot stacked bar plot
phase_cluster_counts.plot(kind='bar', stacked=True, figsize=(12, 6), colormap='Set2')
plt.ylabel('Number of Cells')
plt.title('Cell Cycle Phase Distribution per Cluster')
plt.xlabel('Cluster')
plt.xticks(rotation=45)
plt.legend(title='Cell Cycle Phase')
plt.tight_layout()
plt.show()


import matplotlib.pyplot as plt
import scanpy as sc


cluster_key = key_for_umap # e.g., 'leiden_0.6'
clusters = adata_gaba.obs[cluster_key].cat.categories if hasattr(adata_gaba.obs[cluster_key], 'cat') else adata_gaba.obs[cluster_key].unique()


n_clusters = len(clusters)
ncols = 4
nrows = (n_clusters + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), squeeze=False)


for i, cluster in enumerate(sorted(clusters)):
    ax = axes[i // ncols, i % ncols]
    # Subset adata to cells in the cluster
    cluster_mask = adata_gaba.obs[cluster_key] == cluster
    sc.pl.umap(
        adata_gaba[cluster_mask, :],
        color='cell_cycle_phase',
        palette='Set2',
        title=f'Cluster {cluster}',
        show=False,
        ax=ax,
        legend_loc='on data'
)

# Hide empty subplots if any
for j in range(i + 1, nrows * ncols):
    axes[j // ncols, j % ncols].axis('off')

plt.tight_layout()
plt.show()


import pandas as pd
import matplotlib.pyplot as plt


cluster_key = key_for_umap # your cluster key, e.g., 'leiden_0.6'
genes = all_markers  # list of marker genes you want to plot


# Calculate mean expression of each gene per cluster (using log-normalized layer)
# FIXED:
expr_df = pd.DataFrame(
    data=adata_gaba[:, genes].X.toarray() if hasattr(adata_gaba[:, genes].X, 'toarray') else adata_gaba[:, genes].X,
    columns=genes,
    index=adata_gaba.obs_names
)


expr_df[cluster_key] = adata_gaba.obs[cluster_key].values


# Group by cluster and calculate mean expression per gene
mean_expr_per_cluster = expr_df.groupby(cluster_key).mean()


# Plot stacked bar plot
mean_expr_per_cluster.plot(kind='bar', stacked=True, figsize=(14, 7), colormap='tab20')
plt.ylabel('Average Expression (log normalized)')
plt.title('Stacked Barplot: Marker Gene Expression by Cluster')
plt.xlabel('Cluster')
plt.xticks(rotation=45)
plt.legend(title='Genes', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()


# 1. Dotplot for lineage_markers grouped by chosen clustering resolution
sc.pl.dotplot(
    adata_gaba,
    unique_top_genes,
    groupby=chosen_key,
    layer=None, # Use default layer (X)
    standard_scale='var',
    title=f"Dotplot of Markers",
    swap_axes=True,
    show=True
)


# 2. Stacked bar plot for proportions of WT and KO cells in each cluster


# Count cells per cluster and condition
cluster_condition_counts = adata_gaba.obs.groupby([chosen_key, 'condition']).size().unstack(fill_value=0)


# Calculate proportions (normalizing each cluster total to 1)
cluster_condition_proportions = cluster_condition_counts.div(cluster_condition_counts.sum(axis=1), axis=0)


plt.figure(figsize=(12, 6))
cluster_condition_proportions.plot(
    kind='bar',
    width=0.8,
    colormap='Paired'
)
plt.ylabel('Proportion of Cells')
plt.title(f'Proportion of WT vs HT Cells per Cluster')
plt.xticks(rotation=45, ha='right')
plt.legend(title='Condition')
plt.tight_layout()
plt.show()


import matplotlib.pyplot as plt

# Count cells per cluster and condition
cluster_condition_counts = (
    adata_gaba.obs
    .groupby([chosen_key, 'condition'])
    .size()
    .unstack(fill_value=0)
)

# Calculate proportions within each cluster
cluster_condition_proportions = cluster_condition_counts.div(
    cluster_condition_counts.sum(axis=1),
    axis=0
)

# Optional: set cluster order if they are numeric strings
cluster_condition_proportions.index = cluster_condition_proportions.index.astype(str)
cluster_order = sorted(cluster_condition_proportions.index, key=lambda x: int(x))
cluster_condition_proportions = cluster_condition_proportions.loc[cluster_order]

# Rename legend labels exactly how you want
cluster_condition_proportions = cluster_condition_proportions.rename(
    columns={
        "WT": "WT",
        "HT": r"$\it{Klf13}$ +/-",
        "Klf13 +/-": r"$\it{Klf13}$ +/-"
    }
)

# Choose soft colors
colors = ["#9ecae1", "#fdae6b"]   # light blue, light orange

fig, ax = plt.subplots(figsize=(12, 6))

cluster_condition_proportions.plot(
    kind="bar",
    stacked=True,           # this makes them on top of each other
    width=0.8,
    color=colors,
    edgecolor="none",       # remove outlines
    linewidth=0,
    ax=ax
)

ax.set_ylabel("Proportion of Cells")
ax.set_xlabel("Cluster")
ax.set_title("Proportion of WT vs Klf13 +/- Cells per Cluster")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

# Clean legend
ax.legend(
    title="Condition",
    frameon=False
)

# Optional: remove top/right spines for cleaner look
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.show()


import matplotlib.pyplot as plt

# Count cells per cluster and condition
cluster_condition_counts = (
    adata_gaba.obs
    .groupby([chosen_key, 'condition'])
    .size()
    .unstack(fill_value=0)
)

# Calculate proportions within each cluster
cluster_condition_proportions = cluster_condition_counts.div(
    cluster_condition_counts.sum(axis=1),
    axis=0
)

# Optional: numeric cluster order
cluster_condition_proportions.index = cluster_condition_proportions.index.astype(str)
cluster_order = sorted(cluster_condition_proportions.index, key=lambda x: int(x))
cluster_condition_proportions = cluster_condition_proportions.loc[cluster_order]

# Rename legend labels
cluster_condition_proportions = cluster_condition_proportions.rename(
    columns={
        "WT": "WT",
        "KO": r"$\it{Klf13}$ -/-",
        "Klf13 -/-": r"$\it{Klf13}$ -/-"
    }
)

# Soft colors
colors = ["#9ecae1", "#fdd0a2"]   # light blue, light peach

fig, ax = plt.subplots(figsize=(8, 4.5))

cluster_condition_proportions.plot(
    kind="bar",
    stacked=True,
    width=0.75,
    color=colors,
    edgecolor="none",
    linewidth=0,
    ax=ax
)

ax.set_ylabel("Proportion of Cells")
ax.set_xlabel("Cluster")
ax.set_title("WT and Klf13 -/- proportions by cluster", pad=10)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

# Clean plot appearance
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(False)

# Make plotting area smaller and centered
ax.set_position([0.12, 0.2, 0.55, 0.65])

# Put legend away from plot
ax.legend(
    title="Condition",
    frameon=False,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5)
)

plt.show()

import matplotlib.pyplot as plt

# cluster_names = {
#     "0": "mge_prog_prdm16",
#     "1": "cge_prog_nr2f2",
#     "2": "mge_neu_sox6",
#     "3": "lge_neu_ddah1",
#     "4": "mge_neu_st18",
#     "5": "lge_prog_ebf1",
#     "6": "lge_prog_isl1",
#     "7": "mge_prog_lhx6",
#     "8": "mge_neu_erbb4",
#     "9": "mge_prog_sp9",
#     "10": "cge_prog_sp8",
#     "11": "cge_prog_pax6",
#     "12": "mge_prog_nkx21",
#     "13": "cge_neu_synpr",
#     "14": "cge_prog_ascl1",
# }

import matplotlib.pyplot as plt

# New cluster names (0–13)
cluster_names = {
    "0":  "cge_prog_nr2f1",
    "1":  "mge_prog_lhx6",
    "2":  "lge_neu_ddah1",
    "3":  "cge_prog_sp8",
    "4":  "mge_prog_sp9",
    "5":  "cge_prog_ascl1",
    "6":  "mge_prog_prdm16",
    "7":  "mge_neu_errb4",
    "8":  "lge_prog_isl1",
    "9":  "mge_neu_foxp1",
    "10": "cge_prog_nr2f2",
    "11": "cge_neu_synpr",
    "12": "cge_prog_pax6",
    "13": "lge_prog_ebf1",   # keep this key around for renaming
}

# Count cells per cluster and condition
cluster_condition_counts = (
    adata_gaba.obs
    .groupby([chosen_key, 'condition'])
    .size()
    .unstack(fill_value=0)
)

# Convert cluster IDs to string
cluster_condition_counts.index = cluster_condition_counts.index.astype(str)

# Remove cluster "13"
cluster_condition_counts = cluster_condition_counts.loc[
    cluster_condition_counts.index != "13"
]

# Order clusters numerically (0,1,...,12 only, since 13 is gone)
cluster_order = [str(i) for i in range(14) if str(i) in cluster_condition_counts.index]
cluster_condition_counts = cluster_condition_counts.loc[cluster_order]

# Rename cluster numbers to cluster names (only for clusters that remain)
cluster_condition_counts = cluster_condition_counts.rename(index=cluster_names)

# Calculate proportions within each cluster
cluster_condition_proportions = cluster_condition_counts.div(
    cluster_condition_counts.sum(axis=1),
    axis=0
)

# Rename legend labels (WT and Klf13+/-)
cluster_condition_proportions = cluster_condition_proportions.rename(
    columns={
        "WT": "WT",
        "HT": r"$\it{Klf13}$ +/-",
        "Klf13 +/-": r"$\it{Klf13}$ +/-"
    }
)

# Soft colors
colors = ["#9ecae1", "#fdd0a2"]

fig, ax = plt.subplots(figsize=(8, 4.5))

cluster_condition_proportions.plot(
    kind="bar",
    stacked=True,
    width=0.75,
    color=colors,
    edgecolor="none",
    linewidth=0,
    ax=ax
)

ax.set_ylabel("Proportion of Cells")
ax.set_xlabel("Cluster")
ax.set_title(r"WT vs $\it{Klf13}$ +/- proportions by cluster", pad=10)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

# Clean look
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(False)

# Make the plotting area smaller and centered
ax.set_position([0.12, 0.2, 0.55, 0.65])

# Put legend away from the plot
ax.legend(
    title="Condition",
    frameon=False,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5)
)

plt.show()

# Count cells per cluster and condition
cluster_condition_counts = (
    adata_gaba.obs
    .groupby([chosen_key, 'condition'])
    .size()
    .unstack(fill_value=0)
)

# Convert cluster IDs to string
cluster_condition_counts.index = cluster_condition_counts.index.astype(str)

# Remove cluster 14
cluster_condition_counts = cluster_condition_counts.loc[
    cluster_condition_counts.index != "14"
]

# Order clusters numerically
cluster_order = [str(i) for i in range(15) if i != 14 and str(i) in cluster_condition_counts.index]
cluster_condition_counts = cluster_condition_counts.loc[cluster_order]

# Rename cluster numbers to cluster names
cluster_condition_counts = cluster_condition_counts.rename(index=cluster_names)

# Calculate proportions within each cluster
cluster_condition_proportions = cluster_condition_counts.div(
    cluster_condition_counts.sum(axis=1),
    axis=0
)

# Rename legend labels
cluster_condition_proportions = cluster_condition_proportions.rename(
    columns={
        "WT": "WT",
        "KO": r"$\it{Klf13}$ -/-",
        "Klf13 -/-": r"$\it{Klf13}$ -/-"
    }
)

# Soft colors
colors = ["#9ecae1", "#fdd0a2"]

fig, ax = plt.subplots(figsize=(8, 4.5))

cluster_condition_proportions.plot(
    kind="bar",
    stacked=True,
    width=0.75,
    color=colors,
    edgecolor="none",
    linewidth=0,
    ax=ax
)

ax.set_ylabel("Proportion of Cells")
ax.set_xlabel("Cluster")
ax.set_title(r"WT vs $\it{Klf13}$ -/- proportions by cluster", pad=10)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

# Clean look
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(False)

# Make the plotting area smaller and centered
ax.set_position([0.12, 0.2, 0.55, 0.65])

# Put legend away from the plot
ax.legend(
    title="Condition",
    frameon=False,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5)
)

plt.show()


#====================================================================================

import pandas as pd
import matplotlib.pyplot as plt

cluster_key = key_for_umap      # e.g. 'leiden_res0_5'
file_key    = "batch"           # per-animal / per-file ID

# Count cells per (cluster, batch)
cluster_file_counts = pd.crosstab(
    adata_gaba.obs[cluster_key],
    adata_gaba.obs[file_key]
)

# Absolute counts
cluster_file_counts.plot(
    kind="bar", stacked=True, figsize=(12, 6), colormap="tab20"
)
plt.ylabel("Number of cells")
plt.xlabel("Cluster")
plt.title("Cells from each batch/animal per cluster")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Batch / animal", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.show()

adata_gaba.write("gaba_reclust_scran_new_05.h5ad")
print("All done!")
