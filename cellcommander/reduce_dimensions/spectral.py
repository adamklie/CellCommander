""" Feature selection strategies that can all be run with simple snapatac2 commands. """

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


def run_snapatac2_spectral(
    adata: AnnData,
    features_key: Optional[str] = "highly_variable",
    n_comps: int = 50,
    random_state: Optional[int] = 1234,
):
    logger.info("Running snapatac2 default dimensionality reduction.")
    logger.info("This assumes that '.X' is the layer to use.")
    snap.tl.spectral(adata, n_comps=n_comps, features=features_key, random_state=random_state)
