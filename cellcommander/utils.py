
import logging
import os

import matplotlib
from anndata import AnnData

matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

logger = logging.getLogger("cellcommander")


def describe_anndata(
    adata: AnnData,
):
    """Log important aspects of an AnnData object in a visually pleasing way (e.g. how many cells, features, etc.)"""
    logger.info(f"AnnData object has {adata.n_obs} cells and {adata.n_vars} features")
