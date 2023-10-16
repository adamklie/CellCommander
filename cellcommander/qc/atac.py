import logging
import os

import matplotlib
import numpy as np
import scanpy as sc
from muon import atac as ac
import seaborn as sns
from anndata import AnnData
from cellcommander.qc import consts

matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')
from cellcommander.qc.utils import is_outlier

logger = logging.getLogger("cellcommander")


def plot_nucleosome_signal():
    sns.histplot(atac.obs, x="nucleosome_signal")
    plt.title("Distribution of the nucleome signal")
    plt.show()

def adata_qc(
    adata: AnnData,
    n_for_ns_calc: int = 1e4,
    n_top_features: int = 20,
    total_counts_bins: int = 100,
    total_counts_nmads: int = 3,
    n_features_by_counts_nmads: int = 3,
    pct_counts_in_top_features_nmads: int = 3,
    pct_counts_mt_nmads: int = 3,
    pct_counts_mt_threshold: float = 10,
    n_features_high_threshold: Optional[int] = None,
    n_features_low_threshold: Optional[int] = None,
    min_cells: Optional[int] = None,
) -> AnnData:

    # If we have more than a threshold number cells, filter cells with less than 500 total counts
    if adata.n_obs > consts.DEFAULT_MAX_CELL_COUNT:
        logger.info(
            f"Filtering cells with less than 500 counts due to high initial cell count of {adata.n_obs}"
        )
        sc.pp.filter_cells(adata, min_counts=500)
        logger.info(f"Number of cells after filtering of cells with less than 20 genes: {adata.n_obs}")

    # We need to calculate these if the doublet detection command wasn't run
    if "log_total_fragment_counts" not in adata.obs:

        # Calculate general qc metrics using scanpy
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

        # Rename columns
        adata.obs.rename(
        columns={
                "n_genes_by_counts": "n_features_per_cell",
                "total_counts": "total_fragment_counts",
            },
            inplace=True,
        )

        # log-transform total counts and add as column
        adata.obs["log_total_fragment_counts"] = np.log10(adata.obs["total_fragment_counts"])

    # Calculate nucleosome signal
    ac.tl.nucleosome_signal(adata, n=n_for_ns_calc)

    

    # Add group labels for above and below the nucleosome signal threshold
    nuc_signal_threshold = 2
    atac.obs["nuc_signal_filter"] = [
        "NS_FAIL" if ns > nuc_signal_threshold else "NS_PASS"
        for ns in atac.obs["nucleosome_signal"]
    ]

    # Print number cells not passing nucleosome signal threshold
    atac.obs["nuc_signal_filter"].value_counts()

    # Plot fragment size distribution
    p1 = ac.pl.fragment_histogram(
        atac[atac.obs["nuc_signal_filter"] == "NS_PASS"], region="chr1:1-2000000"
    )

    p2 = ac.pl.fragment_histogram(
        atac[atac.obs["nuc_signal_filter"] == "NS_FAIL"], region="chr1:1-2000000"
    )

    tss = ac.tl.tss_enrichment(mdata, n_tss=10000, random_state=13)
    atac.obs["tss_score"] = tss.obs["tss_score"]

    fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

    p1 = sns.histplot(atac.obs, x="tss_score", ax=axs[0])
    p1.set_title("Full range")

    p2 = sns.histplot(
        atac.obs,
        x="tss_score",
        binrange=(0, atac.obs["tss_score"].quantile(0.995)),
        ax=axs[1],
    )
    p2.set_title("Up to 99.5% percentile")

    plt.suptitle("Distribution of the TSS score")

    plt.tight_layout()
    plt.show()

    tss_threshold = 1.5
    tss.obs["tss_filter"] = [
        "TSS_FAIL" if score < tss_threshold else "TSS_PASS"
        for score in atac.obs["tss_score"]
    ]

    # Print number cells not passing nucleosome signal threshold
    tss.obs["tss_filter"].value_counts()

    # Temporarily set different color palette
    sns.set_palette(palette="Set1")
    ac.pl.tss_enrichment(tss, color="tss_filter")
    # reset color palette
    sns.set_palette(palette="tab10")