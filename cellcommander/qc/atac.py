import logging
import os
import sys
import traceback

import matplotlib
import numpy as np
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from cellcommander.qc import consts
from muon import atac as ac

matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')
from cellcommander.qc.utils import is_outlier

logger = logging.getLogger("cellcommander")


def atac_qc_triplet_plot(
    adata: AnnData,
    outdir_path: str,
    output_prefix: str,
    total_counts_bins: int = 100,
    total_counts_low: int = None,
    total_counts_hi: int = None,
    tss_low: float = None,
    tss_hi: float = None,
):
    with plt.rc_context():
        sns.displot(adata.obs["total_counts"], bins=total_counts_bins, kde=False)
        if total_counts_low:
            plt.axvline(total_counts_low, color="red", linestyle="--")
        if total_counts_hi:
            plt.axvline(total_counts_hi, color="red", linestyle="--")
        plt.savefig(
            os.path.join(outdir_path, f"{output_prefix}_total_counts_distribution.png")
        )
        plt.close()
        
        p1 = sc.pl.violin(adata, "tss_score", show=False)
        if tss_low:
            p1.axhline(y=tss_low, c="red")
        if tss_hi:
            p1.axhline(y=tss_hi, c="red")
        plt.savefig(os.path.join(outdir_path, f"{output_prefix}_tss_score_violin.png"))
        plt.close()

        p2 = sc.pl.scatter(
            adata, "total_counts", "n_features_by_counts", color="tss_score", show=False
        )
        if total_counts_low:
            p2.axvline(x=total_counts_low, c="red")  # Add horizontal line
        if total_counts_hi:
            p2.axvline(x=total_counts_hi, c="red")  # Add vertical line
        plt.savefig(
            os.path.join(outdir_path, f"{output_prefix}_total_counts_vs_n_features.png")
        )
        plt.close()


def plot_nucleosome_signal(
    adata: AnnData,
    outdir_path: str,
    output_prefix: str,
    ns_hi: int = 2,
):
    try:
        _, axs = plt.subplots(1, 2, figsize=(7, 3.5))
        sc.pl.violin(adata, "nucleosome_signal", show=False, ax=axs[0])
        axs[0].axhline(y=ns_hi, c="red")
        sns.histplot(adata.obs, x="nucleosome_signal", ax=axs[1])
        axs[1].axvline(x=ns_hi, color="red", linestyle="--")
        plt.suptitle("Distribution of the nucleome signal")
        plt.savefig(
            os.path.join(
                outdir_path, f"{output_prefix}_nucleosome_signal_distribution.png"
            )
        )
        plt.close()

    except Exception as e:
        logger.error(f"Error plotting nucleosome signal: {e}")
        traceback.logger.info_exc(file=sys.stdout)


def plot_per_barcode_tsse(
    adata: AnnData,
    outdir_path: str,
    output_prefix: str,
    tss_low: int = 1.5,
    tss_hi: int = 200,
):
    try:
        # TSS score distribution plot
        _, axs = plt.subplots(1, 3, figsize=(10.5, 3.5))

        # Violin plot
        p1 = sc.pl.violin(adata, "tss_score", show=False, ax=axs[0])
        p1.set_ylim(0, 200)
        p1.axhline(y=tss_hi, c="red")

        # Histogram plots
        p2 = sns.histplot(adata.obs, x="tss_score", ax=axs[1])
        p2.axvline(tss_low, color="red", linestyle="--")
        p2.set_title("Full range")

        p3 = sns.histplot(
            adata.obs,
            x="tss_score",
            binrange=(0, adata.obs["tss_score"].quantile(0.995)),
            ax=axs[2],
        )
        p3.axvline(tss_low, color="red", linestyle="--")
        p3.axvline(tss_hi, color="red", linestyle="--")
        p3.set_title("Up to 99.5% percentile")

        plt.suptitle("Distribution of the TSS score")
        plt.tight_layout()
        plt.savefig(
            os.path.join(outdir_path, f"{output_prefix}_tss_score_distribution.png")
        )
        plt.close()

    except Exception as e:
        logger.error(f"Error plotting tss enrichment: {e}")
        traceback.logger.info_exc(file=sys.stdout)


def plot_tsse_pass_fail(
    adata: AnnData,
    outdir_path: str,
    output_prefix: str,
):
    sns.set_palette(palette="Set1")
    ac.pl.tss_enrichment(adata, color="tss_filter")
    sns.set_palette(palette="tab10")
    plt.savefig(os.path.join(outdir_path, f"{output_prefix}_tss_score_pass_fail.png"))
    plt.close()


def plot_tsse_logcounts_density(
    adata: AnnData,
    outdir_path: str,
    output_prefix: str,
    total_count_low: int = 1000,
    total_count_hi: int = 50000,
    tss_low: int = 1.5,
):
    # Plot log counts vs tss score# Scatter plot & histograms
    g = sns.jointplot(
        data=adata[(adata.obs["tss_score"] < 20)].obs,
        x="log_total_counts",
        y="tss_score",
        color="black",
        marker=".",
    )
    # Density plot including lines
    g.plot_joint(sns.kdeplot, fill=True, cmap="Blues", zorder=1, alpha=0.75)
    g.plot_joint(sns.kdeplot, color="black", zorder=2, alpha=0.75)

    # Grab the axis
    ax = g.ax_joint

    # Lines thresholds
    ax.axvline(x=np.log10(total_count_hi), c="red", zorder=1)
    ax.axvline(x=np.log10(total_count_low), c="red", zorder=1)
    ax.axhline(y=tss_low, c="red")

    # Save
    plt.savefig(
        os.path.join(outdir_path, f"{output_prefix}_log_total_counts_vs_tss_score.png")
    )
    plt.close()


def get_tss_annotation():
    try:
        import pybiomart as pbm
    except ImportError:
        raise ImportError(
            "Please install pybiomart to use the tss_enrichment function (pip install pybiomart) or use the --no-tss flag."
        )
    # Annotation has to contain columns: Chromosome, Start, End
    dataset = pbm.Dataset(name="hsapiens_gene_ensembl", host="http://www.ensembl.org")
    annot = dataset.query(
        attributes=[
            "chromosome_name",
            "transcription_start_site",
            "strand",
            "external_gene_name",
            "transcript_biotype",
        ]
    )
    filter = annot["Chromosome/scaffold name"].str.contains("CHR|GL|JH|MT")
    annot = annot[~filter]
    annot["Chromosome/scaffold name"] = annot["Chromosome/scaffold name"].str.replace(
        r"(\b\S)", r"chr\1"
    )
    annot.columns = ["Chromosome", "Start", "Strand", "Gene", "Transcript_type"]
    annot = annot[annot.Transcript_type == "protein_coding"]
    annot = annot[
        annot.Chromosome.isin(
            ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
        )
    ]
    return annot


def atac_qc(
    adata: AnnData,
    n_top_features: int = 20,
    pct_counts_in_top_features_nmads: int = 3,
    n_for_ns_calc: int = 1e4,
    ns_nmads: int = 3,
    ns_hi: Optional[int] = None,
    n_tss: int = 1e4,
    tss_annot: Optional[str] = None,
    tss_nmads: int = 3,
    tss_low: Optional[int] = None,
    tss_hi: Optional[int] = None,
    total_counts_nmads: int = 5,
    total_counts_low: Optional[int] = None,
    total_counts_hi: Optional[int] = None,
    n_features_nmads: int = 5,
    n_features_low: Optional[int] = None,
    random_state: int = 13,
) -> AnnData:
    # If we have more than a threshold number cells, filter cells with less than 500 counts
    if adata.n_obs > consts.DEFAULT_MAX_CELL_COUNT_FOR_PRELIM_FILTER:
        logger.info(
            f"Filtering cells with less than 500 fragments due to high initial cell count of {adata.n_obs}"
        )
        sc.pp.filter_cells(adata, min_counts=500)
        logger.info(
            f"Number of cells after filtering of cells with less than 500 fragments: {adata.n_obs}"
        )

    # We need to calculate these if the doublet detection command wasn't run
    if "log_total_counts" not in adata.obs:
        logger.info("Calculating QC metrics")
        # Calculate general qc metrics using scanpy
        sc.pp.calculate_qc_metrics(
            adata, percent_top=[n_top_features], log1p=False, inplace=True
        )

        # Rename columns
        adata.obs.rename(
            columns={
                "n_genes_by_counts": "n_features_by_counts",
                f"pct_counts_in_top_{n_top_features}_genes": "pct_counts_in_top_20_features",
            },
            inplace=True,
        )

        # log-transform total counts and add as column
        adata.obs["log_total_counts"] = np.log10(adata.obs["total_counts"])
    
    # Calculate nucleosome signal
    logger.info(f"Calculating nucleosome signal with {n_for_ns_calc} fragments per cell")
    ac.tl.nucleosome_signal(adata, n=adata.n_obs * n_for_ns_calc)

    # Add group labels for above and below the nucleosome signal threshold
    adata.obs["nuc_signal_filter"] = [
        "NS_FAIL" if ns > ns_hi else "NS_PASS" for ns in adata.obs["nucleosome_signal"]
    ]

    # Calculate tss enrichment
    if tss_annot is None:
        logger.info("No TSS annotation provided, using protein coding genes from default Ensembl gene annotation")
        tss_annot = get_tss_annotation()
    logger.info(f"Calculating TSS enrichment with {n_tss} TSSs")
    tss = ac.tl.tss_enrichment(
        adata, features=tss_annot, n_tss=n_tss, random_state=random_state
    )
    adata.obs["tss_score"] = tss.obs["tss_score"]

    # Create a filter column for TSS score
    tss.obs["tss_filter"] = [
        "TSS_FAIL" if score < tss_low else "TSS_PASS" for score in tss.obs["tss_score"]
    ]

    # Define outliers either by MADs or by user-defined thresholds
    if total_counts_hi is not None or total_counts_low is not None:
        if total_counts_hi is None:
            total_counts_hi = np.inf
        if total_counts_low is None:
            total_counts_low = -np.inf
        logger.info(
            "Defining general feature-based outliers based on user-defined thresholds "
            f"of < {total_counts_low} and > {total_counts_hi} for total fragment counts "
            f"and < {n_features_low} for number of features",
        )
        adata.obs["outlier"] = (
            (adata.obs["total_counts"] < total_counts_low)
            | (adata.obs["total_counts"] > total_counts_hi)
            | (adata.obs["n_features_by_counts"] < n_features_low)
        )
    else:
        logger.info(
            "Defining general feature-based outliers based on MADs"
            f"of > {total_counts_nmads} MADs from the median for total fragment counts, "
            f"> {n_features_nmads} MADs from the median for number of features, "
            f"and > {pct_counts_in_top_features_nmads} MADs from the median for percent counts in top features"
        )
        adata.obs["outlier"] = (
            is_outlier(adata, "total_counts", total_counts_nmads)
            | is_outlier(adata, "n_features_by_counts", n_features_nmads)
            | is_outlier(
                adata, f"pct_counts_in_top_{n_top_features}_features", pct_counts_in_top_features_nmads
            )
        )

    # Define ATAC outliers by either MADs or a user-defined threshold
    if ns_hi is not None or tss_low is not None or tss_hi is not None:
        logger.info(
            "Defining ATAC specific QC outliers based on user-defined thresholds "
            f"of > {ns_hi} for nucleosome signal, < {tss_low} for TSS enrichment, and > {tss_hi} for TSS enrichment"
        )
        adata.obs["atac_outlier"] = (
            (adata.obs["nucleosome_signal"] > ns_hi)
            | (adata.obs["tss_score"] < tss_low)
            | (adata.obs["tss_score"] > tss_hi)
        )
    else:
        logger.info(
            "Defining ATAC specific QC outliers based on MADs "
            f"of > {ns_nmads} MADs from the median for nucleosome signal and "
            f"> {tss_nmads} MADs from the median for TSS enrichment"
        )
        adata.obs["atac_outlier"] = (
            is_outlier(adata, "nucleosome_signal", ns_nmads)
        ) | (is_outlier(adata, "tss_score", tss_nmads))

    return adata


def atac_outlier_filter(
    adata: AnnData,
    outlier_cols: Iterable[str] = ["outlier", "atac_outlier"],
) -> Tuple[AnnData, Iterable[str]]:
    """Filter RNA data based on QC metrics.

    Args:
        adata: Annotated data matrix.
        outlier_cols: Columns in `adata.obs` to use for filtering.
            Defaults to ["outlier", "atac_outlier"]

    Returns:
        Filtered annotated data matrix and the filtered barcodes as a tuple.
    """
    # Actually filter the data
    logger.info(f"Number of cells pre-filtering: {adata.n_obs}")
    # filter_mask = (~adata.obs[outlier_cols].any(axis=1))
    filter_mask = (~adata.obs.outlier) & (~adata.obs.atac_outlier)
    filtered_bc = adata.obs[~filter_mask].index
    adata = adata[filter_mask].copy()
    logger.info(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
    return adata, filtered_bc
