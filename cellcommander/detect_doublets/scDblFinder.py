import logging
import os

import matplotlib
import scanpy as sc
import anndata2ri
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr
from anndata import AnnData
anndata2ri.activate()
anndata2ri.scipy2ri.activate()
matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

logger = logging.getLogger("cellcommander")


def run_scDblFinder(
    adata: AnnData,
    atac_params: bool,
    random_state: int,
):
    logger.info("Importing necessary R modules for scDblFinder.")
    scDblFinder = importr("scDblFinder")
    SingleCellExperiment = importr("SingleCellExperiment")
    logger.info("Running scDblFinder for doublet detection.")
    data_mat = adata.X.T.astype("float32")  # needed to avoid error in anndata2ri v1.1)
    ro.r("set.seed({})".format(random_state))
    if atac_params:
        logger.info("Setting scDblFinder parameters to defaults for scATAC-seq data.")
        ro.globalenv["data_mat"] = data_mat
        ro.r('''sce <- scDblFinder(SingleCellExperiment(list(counts=data_mat)), clusters=TRUE, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")
                scDblFinder_df <- as.data.frame(colData(sce))
             ''')
        scDblFinder_df = ro.globalenv["scDblFinder_df"]
    else:
        logger.info("Setting scDblFinder parameters to defaults for scRNA-seq data.")
        scDblFinder_obj = scDblFinder.scDblFinder(SingleCellExperiment.SingleCellExperiment(assays=ro.ListVector({"counts": data_mat})))
        scDblFinder_df = anndata2ri.rpy2py(scDblFinder_obj.do_slot("colData"))
    adata.obs["scDblFinder_doublet_score"] = scDblFinder_df["scDblFinder.score"].values
    adata.obs["scDblFinder_doublet_class"] = scDblFinder_df["scDblFinder.class"].values
    adata.obs["scDblFinder_predicted_doublet"] = (adata.obs["scDblFinder_doublet_class"] == "doublet")


def plot_scDblFinder(
    adata: AnnData,
    outdir_path: str,
):
    """Plot the scDblFinder scores as histogram"""
    logger.info("Generating scDblFinder score plots.")
    with plt.rc_context():
        _, ax = plt.subplots(figsize=(4, 4))
        ax.hist(adata.obs["scDblFinder_doublet_score"], bins=100)
        ax.set_xlabel("scDblFinder score")
        ax.set_ylabel("Number of cells")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir_path, "scDblFinder_score_distribution.png"))
        plt.close()


def save_scDblFinder(
    adata: AnnData,
    outdir_path: str,
):
    # Add in the scDblFinder scores
    logger.info("Saving scDblFinder scores and predicted doublet barcodes")
    doublet_df = adata.obs[['scDblFinder_doublet_score', 'scDblFinder_doublet_class', 'scDblFinder_predicted_doublet']]
    doublet_df.to_csv(os.path.join(outdir_path, "scDblFinder_doublet_scores.csv"))
    doublet_bcs = doublet_df[doublet_df["scDblFinder_doublet_class"] == "doublet"].index.tolist()
    with open(os.path.join(outdir_path, "scDblFinder_predicted_doublets.txt"), "w") as f:
        for bc in doublet_bcs:
            f.write(bc + "\n")


def scDblFinder_recipe(
    adata: AnnData,
    outdir_path: str,
    random_state: int,
):
    # Run scDblFinder
    run_scDblFinder(
        adata=adata,
        random_state=random_state,
    )

    # Plot scDblFinder scores
    plot_scDblFinder(
        adata=adata,
        outdir_path=outdir_path,
    )

    # Save scDblFinder scores
    save_scDblFinder(
        adata=adata,
        outdir_path=outdir_path,
    )
    