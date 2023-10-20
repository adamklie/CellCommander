import logging
import os

import matplotlib
import numpy as np
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


def run_amulet(
    adata: AnnData,
    outdir_path: str,
):
    # Load packages
    logger.info("Importing necessary R modules for AMULET.")
    scDoubletFinder = importr("scDblFinder")
    SingleCellExperiment = importr("SingleCellExperiment")
    GenomicRanges = importr("GenomicRanges")
    rtracklayer = importr("rtracklayer")
    
    # Set up a GRanges objects of repeat elements, mitochondrial genes and sex chromosomes we want to exclude
    logger.info("Running AMULET for doublet detection.")
    ro.globalenv["frag_path"] = adata.uns["files"]["fragments"]
    ro.r('''repeats <- import('/cellar/users/aklie/opt/AMULET_SourceAndScripts/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed')
         otherChroms <- GRanges(c("chrM","chrX","chrY","MT"),IRanges(1L,width=10^8)) # check which chromosome notation you are using c("M", "X", "Y", "MT")
         toExclude <- suppressWarnings(c(repeats, otherChroms))
         amulet_result <- amulet(frag_path, regionsToExclude=toExclude)
        ''')
    amulet_result = ro.globalenv["amulet_result"]
    
    # Save raw output
    amulet_result.to_csv(os.path.join(outdir_path, "amulet_raw_output.csv"))

     # Merge with barcodes in adata.obs
    amulet_merged_result = adata.obs.merge(amulet_result, left_index=True, right_index=True, how="left")

    # Add columns to adata.obs
    adata.obs["amulet_pVal"] = amulet_merged_result["p.value"]
    adata.obs["amulet_qVal"] = amulet_merged_result["q.value"]
    adata.obs["amulet_predicted_doublet"] = adata.obs["amulet_qVal"] < 0.05
    adata.obs["amulet_negLog10qVal"] = -1 * np.log10(adata.obs["amulet_qVal"])
    adata.obs["amulet_doublet_score"] = adata.obs["amulet_negLog10qVal"]

    # Some are NA, likely cuz they have low counts
    adata.obs["amulet_na"] = adata.obs["amulet_qVal"].isna()


def plot_amulet(
    adata: AnnData,
    outdir_path: str,
):
    """Plot the AMULET scores as histogram"""
    logger.info("Generating AMULET score plots.")
    with plt.rc_context():
        _, ax = plt.subplots(figsize=(4, 4))
        ax.hist(adata.obs["amulet_doublet_score"].clip(upper=10), bins=100)
        ax.set_xlabel("AMULET score")
        ax.set_ylabel("Number of cells")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir_path, "amulet_score_distribution.png"))
        plt.close()


def save_amulet(
    adata: AnnData,
    outdir_path: str,
):
    doublet_df = adata.obs[['amulet_pVal', 'amulet_qVal', 'amulet_negLog10qVal', 'amulet_predicted_doublet', 'amulet_na']]
    doublet_df.to_csv(os.path.join(outdir_path, "amulet_scores.csv"))
    doublet_bcs = doublet_df[doublet_df["amulet_predicted_doublet"] == True].index.tolist()
    with open(os.path.join(outdir_path, "amulet_predicted_doublets.txt"), "w") as f:
        for bc in doublet_bcs:
            f.write(bc + "\n")


def amulet_recipe(
    adata: AnnData,
    outdir_path: str,
):
    # Run AMULET
    run_amulet(
        adata=adata,
        outdir_path=outdir_path,
    )

    # Plot AMULET scores
    plot_amulet(
        adata=adata,
        outdir_path=outdir_path,
    )

    # Save AMULET scores
    save_amulet(
        adata=adata,
        outdir_path=outdir_path,
    )
