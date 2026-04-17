from zarr import group
import scanpy as sc
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt


figdir = "01_subset_genes/outputs/figures"
os.makedirs(figdir,exist_ok=True)
sc.settings.figdir=figdir
sc.set_figure_params(dpi=500,dpi_save=500)
# read panel metadata and object

panel = pd.read_csv("data/Xenium_mMulti_v1_metadata_annotations.csv")
genes = panel["Gene"].tolist()
adata = sc.read_h5ad("data/atlas_light.h5ad")

# ----- subset adata to genes in panel -------

genes_intersection = [g for g in genes if g in adata.var_names]
len(genes_intersection) # 375 genes
len(genes) # 379 genes

adata = adata[:,genes_intersection].copy()
sc.pp.calculate_qc_metrics(adata,percent_top=[50, 100, 200],inplace=True)

# ---- re-process and visualize umap -----

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)

# plot variance per pc

sc.pl.pca_variance_ratio(adata,save="")


# correlation library size ~ PC1

sns.scatterplot(
    x = adata.obs["total_counts"],
    y = adata.obsm["X_pca"][:,0])

fig = plt.gcf()
fig.savefig(os.path.join(figdir,"PC1vsLibsize.png"))



# calculate umap

sc.pp.neighbors(adata,n_pcs=10)
sc.tl.umap(adata, n_components=10)
sc.pl.umap(adata, color=["sample_id","condition","study"])


sc.pl.umap(adata, color=["log1p_total_counts","total_counts"],size=4,save="_qc.png")

sc.pl.umap(adata, color=["cell_type_level0","cell_type_level1"],size=3,ncols=1,save="_cts.png")


# DGE

sc.tl.rank_genes_groups(adata, groupby="cell_type_level0",n_genes = 7,key_added="dge_subset_ct0")
sc.pl.rank_genes_groups_dotplot(adata, groupby="cell_type_level0",n_genes = 7,key="dge_subset_ct0",save="_dotplot_ct0.png")


# save adata
os.makedirs("01_subset_genes/outputs/files")
sc.settings.writedir="01_subset_genes/outputs/files"
sc.write(adata=adata,filename="adatasubset",ext="h5ad")


# remove low count cells

adata = adata[adata.obs["total_counts"] >= 10].copy()
sc.pl.umap(adata, color=["log1p_total_counts","total_counts"],size=4,save="_qc_nolowcells.png")

sc.pp.pca(adata,key_added="pca_nolow")

sc.pl.pca_variance_ratio(adata,save="_nolow")

sc.pp.neighbors(adata,n_pcs=10,use_rep="pca_nolow",key_added="nn_nolow")
sc.tl.umap(adata, n_components=10,neighbors_key="nn_nolow",key_added="umap_nolow")
sc.pl.embedding(adata,basis = "umap_nolow", color=["log1p_total_counts","total_counts"],size=4,save="_qc_umapnolow.png")
sc.pl.embedding(adata,basis="umap_nolow",color="log1p_total_counts")


sc.pl.embedding(adata,basis="umap_nolow",color="cell_type_level0",mask_obs=[ct for ct in adata.obs["cell_type_level0"].values.unique()])
