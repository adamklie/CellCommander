import logging
import os

import matplotlib
import seaborn as sns
from anndata import AnnData

matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

logger = logging.getLogger("cellcommander")


def plot_against_raw(
    adata: AnnData,
    outdir_path: str,
    layer_key: Optional[str] = None,
    obsm_key: Optional[str] = None,
):
    # Plot the side-by-side distribution of total counts before and after normalization
    _, axes = plt.subplots(1, 2, figsize=(10, 5))
    if "total_counts" not in adata.obs.columns:
        adata.obs["total_counts"] = adata.X.sum(1)
    p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
    axes[0].set_title("Total raw counts across cells")
    norm_mat = adata.layers[layer_key] if layer_key is not None else adata.obsm[obsm_key]
    key = layer_key if layer_key is not None else obsm_key
    p2 = sns.histplot(norm_mat.sum(1), bins=100, kde=False, ax=axes[1])
    axes[1].set_title(f"Total {key} counts across cells")
    plt.savefig(os.path.join(outdir_path, f"total_counts_{key}.png"))
    plt.close()
    