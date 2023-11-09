import logging
import os

import matplotlib
import scanpy.external as sce
from anndata import AnnData

matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

logger = logging.getLogger("cellcommander")


def run_scrublet(
    adata: AnnData,
    random_state: Optional[int] = 1234,
):
    # Run scrublet
    logger.info("Running Scrublet for doublet detection.")
    sce.pp.scrublet(
        adata,
        adata_sim=None,
        random_state=random_state
    )


def plot_scrublet(
    adata: AnnData,
    outdir_path: str,
):
    # Plot scrublet scores
    logger.info("Generating Scrublet score plots.")
    with plt.rc_context():
        sce.pl.scrublet_score_distribution(adata, show=False)
        plt.savefig(os.path.join(outdir_path, "scrublet_score_distribution.png"))
        plt.close()


def save_scrublet(
    adata: AnnData,
    outdir_path: str,
):
    # Save scrublet scores
    logger.info("Saving Scrublet scores and predicted doublet barcodes")
    adata.obs.rename(columns={"doublet_score": "scrublet_doublet_score"}, inplace=True)
    adata.obs.rename(columns={"predicted_doublet": "scrublet_predicted_doublet"}, inplace=True)
    doublet_df = adata.obs[['scrublet_doublet_score', 'scrublet_predicted_doublet']]
    doublet_df.to_csv(os.path.join(outdir_path, "scrublet_doublet_scores.csv"))
    doublet_bcs = doublet_df[doublet_df["scrublet_predicted_doublet"] == True].index.tolist()
    with open(os.path.join(outdir_path, "scrublet_predicted_doublets.txt"), "w") as f:
        for bc in doublet_bcs:
            f.write(bc + "\n")


def scrublet_recipe(
    adata: AnnData,
    outdir_path: str,
    random_state: int,
):
    # Run scrublet
    run_scrublet(
        adata=adata,
        random_state=random_state,
    )

    # Plot scrublet scores
    plot_scrublet(
        adata=adata,
        outdir_path=outdir_path,
    )

    # Save scrublet scores
    save_scrublet(
        adata=adata,
        outdir_path=outdir_path,
)