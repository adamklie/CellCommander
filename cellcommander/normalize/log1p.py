import logging
import os

import matplotlib
import numpy as np
import scanpy as sc
from anndata import AnnData
from scipy.io import mmwrite

matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

from cellcommander.normalize.utils import plot_against_raw

logger = logging.getLogger("cellcommander")


def run_log1p(
    adata: AnnData,
    target_sum: Optional[int] = None,
    to_layer: Optional[str] = "log1p_norm",
):
    logger.info("Running shifted logarithm normalization.")
    scales_counts = sc.pp.normalize_total(adata, target_sum=target_sum, inplace=False)
    adata.layers[to_layer] = sc.pp.log1p(scales_counts["X"], copy=True)


def plot_log1p(
    adata: AnnData,
    outdir_path: str,
    layer_key: Optional[str] = "log1p_norm",
):
    logger.info(f"Plotting log1p normalized data against raw data.")
    plot_against_raw(adata, outdir_path, layer_key=layer_key)


def save_log1p(
    adata: AnnData,
    outdir_path: str,
    layer_key: Optional[str] = "log1p_norm",
):
    logger.info(f"Saving log1p normalized data to {os.path.join(outdir_path)}")
    X = adata.layers[layer_key]
    mmwrite(os.path.join(outdir_path, "mtx.mtx"), X)
    adata.obs.index.to_series().to_csv(os.path.join(outdir_path, "barcodes.tsv"), sep="\t", index=False, header=False)
    adata.var.index.to_series().to_csv(os.path.join(outdir_path, "features.tsv"), sep="\t", index=False, header=False)


def log1p_recipe(
    adata: AnnData,
    outdir_path: str,
    layer_key: Optional[str] = "log1p_norm",
    save_normalized_mtx: bool = False,
):
    logger.info("Running log1p normalization.")
    run_log1p(adata, to_layer=layer_key)
    plot_log1p(adata, outdir_path, layer_key=layer_key)
    if save_normalized_mtx:
        if not os.path.exists(os.path.join(outdir_path, "log1p_norm")):
            os.makedirs(os.path.join(outdir_path, "log1p_norm"))
        save_log1p(adata, os.path.join(outdir_path, "log1p_norm"), layer_key=layer_key)
