import sys
sys.path.append("/cellar/users/aklie/opt/igvf-ucsd/single_cell_utilities/data_processing")
from normalizations import proportional_filtering_normalize, log1p_normalize
proportional_filtering_normalize(
    adata=adata,
    layer_out="depth_normalization",
)
log1p_normalize(
    adata=adata,
    layer_in="depth_normalization",
    layer_out="depth_normalization",
)
proportional_filtering_normalize(
    adata=adata,
    layer_in="depth_normalization",
    layer_out="depth_normalization",
)

# Plot the side-by-side distribution of total counts before and after normalization
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(
    adata.layers["depth_normalization"].sum(1), bins=100, kde=False, ax=axes[1]
)
axes[1].set_title("Depth normalization")
plt.savefig(os.path.join(outdir_path, f"depth_normalization_distribution.png"))
plt.close()