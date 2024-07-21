import logging
import os

import matplotlib
import numpy as np
import scanpy as sc
from muon import atac as ac
import anndata2ri
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr
from anndata import AnnData
import pandas as pd
anndata2ri.activate()
anndata2ri.scipy2ri.activate()
matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple, Union, List

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

logger = logging.getLogger("cellcommander")


def run_harmonyR(
    adata: AnnData,
    obsm_key: str = "X_pca",
    corrected_obsm_key: str = "X_harmony",
    vars_to_correct: List = ["sample"],
    max_iter_harmony: int = 10,
    random_state: int = 1234,
):
    # R imports
    harmony = importr("harmony")
    magrittr = importr("magrittr")

    # Grab the dim reduction to correct
    n_pcs = adata.obsm[obsm_key].shape[1]
    dim_reduce = pd.DataFrame(adata.obsm[obsm_key], columns=['PC{}'.format(i) for i in range(1,n_pcs+1)], index=adata.obs.index)
    obs = adata.obs

    # Toss into R's global environment
    ro.globalenv["dim_reduce"] = dim_reduce
    ro.globalenv["obs"] = obs
    ro.globalenv["vars_to_correct"] = ro.StrVector(vars_to_correct)
    ro.globalenv["max_iter_harmony"] = max_iter_harmony
    ro.globalenv["random_state"] = random_state

     # Run 
    ro.r('''set.seed(random_state)
            harmonized <- HarmonyMatrix(dim_reduce, obs, vars_use=vars_to_correct, do_pca=FALSE, max.iter.harmony=max_iter_harmony)
            harmonized <- data.frame(harmonized)
        ''')
    
    # Get the corrected matrix
    harmony_embedding_df = pd.DataFrame(ro.r['harmonized'], index=['PC{}'.format(i) for i in range(1,n_pcs+1)], columns=adata.obs.index).T

    # Add the corrected matrix
    adata.obsm[corrected_obsm_key] = harmony_embedding_df.loc[adata.obs.index].values


def run_harmony(
    adata: AnnData,
    obsm_key: str = "X_pca",
    corrected_obsm_key: str = "X_harmony",
    vars_to_correct: List = ["sample"],
    nclust: int = 100,
    max_iter_harmony: int = 10,
    random_state: int = 1234,
):
    from harmonypy import run_harmony
    ho = run_harmony(
        adata.obsm[obsm_key],
        meta_data=adata.obs,
        vars_use=vars_to_correct,
        nclust=nclust,
        max_iter_harmony=max_iter_harmony,
        random_state=random_state,
    )
    adata.obsm[corrected_obsm_key] = ho.Z_corr.T