import logging
import os

import matplotlib
import numpy as np
import scanpy as sc
import anndata2ri
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr
from anndata import AnnData
anndata2ri.activate()
anndata2ri.scipy2ri.activate()
matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

logger = logging.getLogger("cellcommander")


def run_deviance(
    adata: AnnData,
):
    adata_ = sc.AnnData(adata.X.copy())
    adata_.obs_names = adata.obs_names.copy()
    adata_.var_names = adata.var_names.copy()
    adata_.X = adata_.X.astype(np.float32)
    
    # Set up R imports
    scry = importr("scry")

    # Put the data into R
    ro.globalenv["adata_"] = adata_

    ro.r('''sce = devianceFeatureSelection(adata_, assay="X")
            binomial_deviance = rowData(sce)$binomial_deviance
         ''')
    binomial_deviance = ro.globalenv["binomial_deviance"].T

    # Annotate highly deviant genes
    idx = binomial_deviance.argsort()[-4000:]
    mask = np.zeros(adata.var_names.shape, dtype=bool)
    mask[idx] = True
    adata.var["highly_deviant"] = mask
    adata.var["binomial_deviance"] = binomial_deviance