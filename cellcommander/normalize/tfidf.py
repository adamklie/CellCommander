import logging
import os

import matplotlib
import numpy as np
import scanpy as sc
import anndata2ri
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr
from muon import atac as ac
from anndata import AnnData
from scipy.io import mmwrite
anndata2ri.activate()
anndata2ri.scipy2ri.activate()
matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

from cellcommander.normalize.utils import plot_against_raw

logger = logging.getLogger("cellcommander")


def run_tfidf(
    adata: AnnData,
    scale_factor: float = 1e4,
    to_layer: Optional[str] = "tfidf_norm"
):
    ac.pp.tfidf(adata, scale_factor=scale_factor, to_layer=to_layer)


def run_signac_tfidf(
    adata: AnnData,
    scale_factor: float = 1e4,
    to_layer: Optional[str] = "seurat_tfidf_norm"
):
    
    # Set up R imports
    Seurat = importr("Seurat")
    SeuratObject = importr("SeuratObject")
    Matrix = importr("Matrix")
    Signac = importr("Signac")

    # Prepare data for Seurat
    data_mat = adata.X.T.astype("float32")
    cell_names = adata.obs_names
    gene_names = adata.var_names

    # Load into global environment
    ro.globalenv["data_mat"] = data_mat
    ro.globalenv["cell_names"] = cell_names
    ro.globalenv["gene_names"] = gene_names
    ro.globalenv["scale_factor"] = scale_factor

    # Run
    ro.r('''mtx = Matrix(data_mat, sparse = TRUE)
            rownames(mtx) = gene_names
            colnames(mtx) = cell_names
            sobj = CreateSeuratObject(counts = mtx, assay = "ATAC")
            sobj = RunTFIDF(sobj, assay = "ATAC", scale.factor = scale_factor)
            data = GetAssayData(object = sobj, assay = "ATAC", slot = "data")
            ''')

    # Get data
    data = ro.globalenv["data"]
    adata.layers[to_layer] = data.T
         

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
    logger.info(f"Saving tfidf normalized data to {os.path.join(outdir_path)}")
    X = adata.layers[layer_key]
    mmwrite(os.path.join(outdir_path, "mtx.mtx"), X)
    adata.obs.index.to_series().to_csv(os.path.join(outdir_path, "barcodes.tsv"), sep="\t", index=False, header=False)
    adata.var.index.to_series().to_csv(os.path.join(outdir_path, "features.tsv"), sep="\t", index=False, header=False)


def tfidf_recipe(
    adata: AnnData,
    outdir_path: str,
    scale_factor: float = 1e4,
    layer_key: Optional[str] = "tfidf_norm",
    save_normalized_mtx: bool = False,
):
    logger.info("Running tfidf normalization.")
    run_tfidf(adata, scale_factor=scale_factor, to_layer=layer_key)
    plot_tfidf(adata, outdir_path, layer_key=layer_key)
    if save_normalized_mtx:
        if not os.path.exists(os.path.join(outdir_path, "tfidf_norm")):
            os.makedirs(os.path.join(outdir_path, "tfidf_norm"))
        save_tfidf(adata, os.path.join(outdir_path, "tfidf_norm"), layer_key=layer_key)
