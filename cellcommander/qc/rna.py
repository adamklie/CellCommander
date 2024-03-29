import logging
import os

import matplotlib
import numpy as np
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from cellcommander.qc import consts

matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')
from cellcommander.qc.utils import is_outlier

logger = logging.getLogger("cellcommander")


def rna_qc_triplet_plot(
    adata: AnnData,
    outdir_path: str,
    output_prefix: str,
    total_counts_bins: int = 100,
):
    with plt.rc_context():
        sns.displot(adata.obs["total_counts"], bins=total_counts_bins, kde=False)
        plt.savefig(
            os.path.join(outdir_path, f"{output_prefix}_total_counts_distribution.png")
        )
        plt.close()

        sc.pl.violin(adata, "pct_counts_mt")
        plt.savefig(
            os.path.join(outdir_path, f"{output_prefix}_pct_counts_mt_violin.png")
        )
        plt.close()

        sc.pl.violin(adata, "pct_counts_ribo")
        plt.savefig(
            os.path.join(outdir_path, f"{output_prefix}_pct_counts_ribo_violin.png")
        )
        plt.close()

        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
        plt.savefig(
            os.path.join(outdir_path, f"{output_prefix}_total_counts_vs_n_genes.png")
        )
        plt.close()


def rna_qc(
    adata: AnnData,
    n_top_genes: int = 20,
    total_counts_nmads: int = 3,
    n_genes_by_counts_nmads: int = 3,
    pct_counts_in_top_genes_nmads: int = 3,
    pct_counts_mt_nmads: int = 3,
    pct_counts_mt_hi: float = 10,
    pct_counts_ribo_nmads: int = 3,
    pct_counts_ribo_hi: float = 1,
    n_genes_hi: Optional[int] = None,
    n_genes_low: Optional[int] = None,
) -> None:
    """Perform QC and filtering on RNA data.

    Args:
        adata: Annotated data matrix.
        n_top_genes: Number of top genes to calculate the percentage of total counts in.
        total_counts_nmads: Number of MADs to use for determining outliers for total counts.
        n_genes_by_counts_nmads: Number of MADs to use for determining outliers for number of genes by counts.
        pct_counts_in_top_genes_nmads: Number of MADs to use for determining outliers for percentage of counts in top genes.
        pct_counts_mt_nmads: Number of MADs to use for determining outliers for percentage of counts in mitochondrial genes.
        pct_counts_mt_hi: Threshold for determining outliers for percentage of counts in mitochondrial genes.
        n_genes_hi: Threshold for determining outliers for number of genes by counts.
        n_genes_low: Threshold for determining outliers for number of genes by counts.

    Returns:
        Annotated data matrix with QC metrics calculated and stored in `adata.obs`.
        Importantly, it will add 'outlier' and 'mt_outlier' columns to `adata.obs` that when
        True indicate that the cell is an outlier based on the QC metrics.
        This function does not filter the data and modifies the input `adata` object in place.
    """

    # If we have more than a threshold number cells, filter cells with less than 20 genes
    if adata.n_obs > consts.DEFAULT_MAX_CELL_COUNT_FOR_PRELIM_FILTER:
        logger.info(
            f"Filtering cells with less than 20 genes due to high initial cell count of {adata.n_obs}"
        )
        sc.pp.filter_cells(adata, min_genes=20)
        logger.info(
            f"Number of cells after filtering of cells with less than 20 genes: {adata.n_obs}"
        )

    # Calculate QC metrics
    logger.info("Calculating QC metrics")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        inplace=True,
        percent_top=[n_top_genes],
        log1p=True,
    )

    # Define outliers either by MADs or by user-defined thresholds
    if n_genes_hi is not None or n_genes_low is not None:
        if n_genes_hi is None:
            n_genes_hi = np.inf
        if n_genes_low is None:
            n_genes_low = -np.inf
        logger.info(
            "Defining general feature-based outliers based on user-defined thresholds "
            f"of < {n_genes_low} and > {n_genes_hi}"
        )
        adata.obs["outlier"] = (adata.obs["n_genes_by_counts"] < n_genes_low) | (
            adata.obs["n_genes_by_counts"] > n_genes_hi
        )
    else:
        logger.info(
            "Defining general feature-based outliers based on "
            f"{total_counts_nmads} MADs for total counts, "
            f"{n_genes_by_counts_nmads} MADs for number of genes by counts, "
            f" and {pct_counts_in_top_genes_nmads} MADs for percentage of counts in top genes"
        )
        adata.obs["outlier"] = (
            is_outlier(adata, "log1p_total_counts", total_counts_nmads)
            | is_outlier(adata, "log1p_n_genes_by_counts", n_genes_by_counts_nmads)
            | is_outlier(
                adata, "pct_counts_in_top_20_genes", pct_counts_in_top_genes_nmads
            )
        )

    # Define mt outliers by either MADs or a user-defined threshold
    if n_genes_hi is not None or n_genes_low is not None:
        logger.info(
            "Defining mitochondrial feature outliers based on a user-defined threshold "
            f"of {pct_counts_mt_hi}% of total counts"
        )
        adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > pct_counts_mt_hi
    else:
        logger.info(
            "Defining mitochondrial feature outliers based on "
            f"{pct_counts_mt_nmads} MADs for the percentage of counts in mitochondrial features "
            f"and a user-defined threshold of {pct_counts_mt_hi}% of total counts"
        )
        adata.obs["mt_outlier"] = is_outlier(
            adata, "pct_counts_mt", pct_counts_mt_nmads
        ) | (adata.obs["pct_counts_mt"] > pct_counts_mt_hi)

    # Define ribo outliers by either MADs or a user-defined threshold
    if n_genes_hi is not None or n_genes_low is not None:
        logger.info(
            "Defining ribosomal feature outliers based on a user-defined threshold "
            f"of {pct_counts_ribo_hi}% of total counts"
        )
        adata.obs["ribo_outlier"] = adata.obs["pct_counts_ribo"] > pct_counts_ribo_hi
    else:
        logger.info(
            "Defining ribosomal feature outliers based on "
            f"{pct_counts_ribo_nmads} MADs for the percentage of counts in ribosomal features "
            f"and a user-defined threshold of {pct_counts_ribo_hi}% of total counts"
        )
        adata.obs["ribo_outlier"] = is_outlier(
            adata, "pct_counts_ribo", pct_counts_ribo_nmads
        ) | (adata.obs["pct_counts_ribo"] > pct_counts_ribo_hi)

    return


def rna_outlier_filter(
    adata: AnnData,
    outlier_cols: Iterable[str] = ["outlier", "mt_outlier"],
) -> Tuple[AnnData, Iterable[str]]:
    """Filter RNA data based on QC metrics.

    Args:
        adata: Annotated data matrix.
        outlier_cols: Columns in `adata.obs` to use for filtering.
            Defaults to ["outlier", "mt_outlier"].

    Returns:
        Filtered annotated data matrix and the filtered barcodes as a tuple.
    """
    # Actually filter the data
    logger.info(f"Number of cells pre-filtering: {adata.n_obs}")
    # filter_mask = (~adata.obs[outlier_cols].any(axis=1))
    filter_mask = (~adata.obs.outlier) & (~adata.obs.mt_outlier) & (~adata.obs.ribo_outlier)
    filtered_bc = adata.obs[~filter_mask].index
    adata = adata[filter_mask].copy()
    logger.info(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
    return adata, filtered_bc
