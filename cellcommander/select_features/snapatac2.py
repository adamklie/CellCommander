""" Feature selection strategies that can all be run with simple ScanPy commands. """

import logging
import os

import numpy as np
import matplotlib
import snapatac2 as snap
from anndata import AnnData

matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

logger = logging.getLogger("cellcommander")


def run_snapatac2(
    adata: AnnData,
    num_features: int = 50000,
):
    logger.info("Running Seurat feature selection with ScanPy implementation (sc.pp.highly_variable_genes).")
    snap.pp.select_features(adata, n_features=num_features)
