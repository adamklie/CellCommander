""" Feature selection strategies that can all be run with simple ScanPy commands. """

import logging
import os

import numpy as np
import matplotlib
import scanpy as sc
from anndata import AnnData

matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

logger = logging.getLogger("cellcommander")


def run_seurat(
    adata: AnnData,
    layer: Optional[str] = None,
    min_mean: Optional[float] = 0.125,
    max_mean: Optional[float] = 3,
    min_disp: Optional[float] = 0.5,
    max_disp: Optional[float] = np.inf,
    n_bins: int = 20,
    key_added: Optional[str] = "highly_variable",
):
    logger.info("Running Seurat feature selection with ScanPy implementation (sc.pp.highly_variable_genes).")
    sc.pp.highly_variable_genes(
        adata=adata,
        layer=layer,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        max_disp=max_disp,
        n_bins=n_bins,
        flavor="seurat",
    )
    if key_added != "highly_variable":
        adata.var[key_added] = adata.var["highly_variable"].copy()
        adata.var.drop("highly_variable", axis=1, inplace=True)


def run_seurat_v3(
    adata: AnnData,
    layer: Optional[str] = None,
    n_top_genes: int = 3000,
    span: float = 0.3,
    key_added: Optional[str] = "highly_variable",
):
    logger.info("Running Seurat v3 feature selection with ScanPy implementation (sc.pp.highly_variable_genes).")
    sc.pp.highly_variable_genes(
        adata=adata,
        layer=layer,
        n_top_genes=n_top_genes,
        span=span,
        flavor="seurat_v3",
    )
    if key_added != "highly_variable":
        adata.var[key_added] = adata.var["highly_variable"].copy()
        adata.var.drop("highly_variable", axis=1, inplace=True)


def run_cell_ranger(
    adata: AnnData,
    min_mean: Optional[float] = None,
    max_mean: Optional[float] = None,
    min_disp: Optional[float] = None,
    max_disp: Optional[float] = None,
    n_bins: int = 20,
    key_added: Optional[str] = "highly_variable",
):
    logger.info("Running Cell Ranger feature selection with ScanPy implementation (sc.pp.highly_variable_genes).")
    sc.pp.highly_variable_genes(
        adata=adata,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        max_disp=max_disp,
        n_bins=n_bins,
        flavor="cell_ranger",
    )
    if key_added != "highly_variable":
        adata.var[key_added] = adata.var["highly_variable"].copy()
        adata.var.drop("highly_variable", axis=1, inplace=True)
