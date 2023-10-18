logging.info("Running scran normalization.")
scran = importr("scran")
SingleCellExperiment = importr("SingleCellExperiment")
BiocParallel = importr("BiocParallel")
adata_pp = adata.copy()
sc.pp.filter_genes(adata_pp, min_cells=20)
sc.pp.normalize_total(adata_pp)
sc.pp.log1p(adata_pp)
sc.pp.highly_variable_genes(adata_pp, n_top_genes=3000)
adata_pp = adata_pp[:, adata_pp.var.highly_variable]
sc.pp.pca(adata_pp)
sc.pp.neighbors(adata_pp, n_neighbors=initial_clust_n_neighbors)
sc.tl.umap(adata_pp)
sc.tl.leiden(adata_pp, key_added="scran_groups", resolution=initial_clust_resolution)
gene_markers = [x for x in gene_markers if x in adata_pp.var_names]
with plt.rc_context():
    sc.pl.umap(adata_pp, color=["scran_groups"])
    plt.savefig(os.path.join(outdir_path, "scran_groups_umap.png"))
    plt.close()
with plt.rc_context():
    sc.pl.dotplot(
        adata_pp,
        groupby="scran_groups",
        var_names=gene_markers,
        standard_scale="var",
    )
    plt.savefig(os.path.join(outdir_path, "scran_marker_genes.png"))
    plt.close()
with plt.rc_context():
    sc.pl.umap(
        adata_pp,
        color=gene_markers,
        vmin=0,
        vmax="p99", 
        sort_order=False, 
        frameon=False,
        cmap="Reds", 
    )
    plt.savefig(os.path.join(outdir_path, "scran_marker_genes_umap.png"))
    plt.close()

# convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
data_mat = adata_pp.X.T
if issparse(data_mat):
    if data_mat.nnz > 2**31 - 1:
        data_mat = data_mat.tocoo()
    else:
        data_mat = data_mat.tocsc()
ro.globalenv["data_mat"] = data_mat
ro.globalenv["input_groups"] = adata_pp.obs["scran_groups"]
adata.obs["scran_groups"] = adata_pp.obs["scran_groups"]

# Manage memory
del adata_pp

# Run scran normalization
ro.r('''sce <- SingleCellExperiment(assays=list(counts=data_mat))
        sce <- computeSumFactors(sce, cluster=input_groups, min.mean=0.1, BPPARAM=MulticoreParam())
        sf <- sizeFactors(sce)
    ''')
size_factors = ro.globalenv["sf"]

# Add the scran normalization to the adata object
adata.obs["size_factors"] = size_factors
scran = adata.X / adata.obs["size_factors"].values[:, None]
adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))

# Plot the side-by-side distribution of total counts before and after normalization
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(
    adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1]
)
axes[1].set_title("log1p with Scran estimated size factors")
plt.savefig(os.path.join(outdir_path, f"scran_normalization_distribution.png"))
plt.close()