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


def run_scanpy_default(
    adata: AnnData,
    layer: Optional[str] = None,
    obsm_key: Optional[str] = None,
    n_comps: int = 50,
    random_state: Optional[int] = 1234,
):
    if layer is not None:
        assert layer in adata.layers.keys(), f"Layer {layer} not found in AnnData object."
        assert obsm_key is None, "Cannot specify both layer and obsm_key."
        logger.info(f"Running PCA on layer {layer} (ScanPy default)")
        sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack", use_highly_variable=True, zero_center=True, dtype="float32", copy=False, random_state=random_state)
    elif obsm_key is not None:
        assert obsm_key in adata.obsm.keys(), f"Key {obsm_key} not found in AnnData object."
        logger.info(f"Running PCA on obsm key {obsm_key} (ScanPy default), note that you cannot use highly variable genes with this method.")
        adata.obsm["X_pca"] = sc.pp.pca(adata.obsm[obsm_key].values, n_comps=n_comps, use_highly_variable=False, return_info=False, random_state=random_state)
    else:
        logger.info("Must specify either layer or obsm_key to run ScanPy default PCA.")