import logging
import os

import matplotlib
from muon import atac as ac
from anndata import AnnData
from scipy.io import mmwrite

matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

from cellcommander.normalize.utils import plot_against_raw
logger = logging.getLogger("cellcommander")


def run_tfidf(
    adata: AnnData,
    scale_factor: float = 1e4,
    to_layer: Optional[str] = "tfidf_norm"
):
    ac.pp.tfidf(adata, scale_factor=scale_factor, to_layer=to_layer)


def plot_tfidf(
    adata: AnnData,
    outdir_path: str,
    layer_key: Optional[str] = "tfidf_norm",
):
    plot_against_raw(adata, outdir_path, layer_key=layer_key)


def save_tfidf(
    adata: AnnData,
    outdir_path: str,
    layer_key: Optional[str] = "tfidf_norm",
):
    logger.info(f"Saving tfidf normalized data to {os.path.join(outdir_path, 'tfidf_norm.h5ad')}")
    X = adata.layers[layer_key]
    mmwrite(os.path.join(outdir_path, "tfidf_norm.mtx"), X)
    adata.obs.index.to_series().to_csv(os.path.join(outdir_path, "barcodes.tsv"), sep="\t", index=False, header=False)
    adata.var.index.to_series().to_csv(os.path.join(outdir_path, "features.tsv"), sep="\t", index=False, header=False)


def tfidf_recipe(
    adata: AnnData,
    outdir_path: str,
    scale_factor: float = 1e4,
    layer_key: Optional[str] = "tfidf_norm",
):
    logger.info("Running tfidf normalization.")
    run_tfidf(adata, scale_factor=scale_factor, to_layer=layer_key)
    plot_tfidf(adata, outdir_path, layer_key=layer_key)
    save_tfidf(adata, outdir_path, layer_key=layer_key)
