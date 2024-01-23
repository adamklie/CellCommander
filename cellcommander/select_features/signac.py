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


def run_signac(
    adata: AnnData,
    layer: Optional[str] = None,
    key_added: Optional[str] = "highly_variable",
    min_cutoff: Optional[str] = "q0",
):
    # Import libraries
    logger.info("Importing necessary R libraries for Signac feature selection")
    Seurat = importr("Seurat")
    SeuratObject = importr("SeuratObject")
    Matrix = importr("Matrix")
    Signac = importr("Signac")

    # Prepare data for Seurat
    if layer is None:
        data_mat = adata.X.astype("float32").T
    else:
        data_mat = adata.layers[layer].T.astype("float32")
    cell_names = adata.obs_names
    gene_names = adata.var_names

    # Load into global environment
    ro.globalenv["data_mat"] = data_mat
    ro.globalenv["cell_names"] = cell_names
    ro.globalenv["gene_names"] = gene_names
    ro.globalenv["min_cutoff"] = min_cutoff

    ro.r('''mtx = Matrix(data_mat, sparse = TRUE)
        rownames(mtx) = gene_names
        colnames(mtx) = cell_names
        sobj = CreateSeuratObject(counts = mtx, assay = "ATAC")
        sobj <- FindTopFeatures(sobj, min.cutoff=min_cutoff)
        variable_features = VariableFeatures(object = sobj)
        ''')

    # Get results into Python
    variable_features = ro.globalenv["variable_features"]
    logger.info(f"Added key {key_added} to adata.var with {len(variable_features)} features total")
    adata.var[key_added] = adata.var.index.isin(variable_features)

    