import logging
import os

import matplotlib
import numpy as np
import pandas as pd
import anndata2ri
import scanpy as sc
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr
from anndata import AnnData
anndata2ri.activate()
anndata2ri.scipy2ri.activate()
matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

from cellcommander.normalize.utils import plot_against_raw
logger = logging.getLogger("cellcommander")


def run_sctransform(
    adata,
    filter_genes: bool = False,
):
    logger.info("Importing necessary R modules for running SCTransform with Seurat.")
    Seurat = importr("Seurat")
    SeuratObject = importr("SeuratObject")
    Matrix = importr("Matrix")
    
    # Make a copy of the data to run SCTransform on
    adata_pp = adata.copy()

    # I've seen better results when I filter out genes with low counts
    if filter_genes:
        logger.info("Filtering genes with low counts before running SCTransform.")
        logger.info(f"Total number of genes: {adata_pp.n_vars}")
        sc.pp.filter_genes(adata_pp, min_cells=20)
        logger.info(f"Number of genes after cell filter: {adata_pp.n_vars}")

    # Prepare data for Seurat
    data_mat = adata_pp.X.T.astype("float32")
    cell_names = adata_pp.obs_names
    gene_names = adata_pp.var_names

    # Load into global environment
    ro.globalenv["data_mat"] = data_mat
    ro.globalenv["cell_names"] = cell_names
    ro.globalenv["gene_names"] = gene_names

    # Run SCTransform
    logger.info("Running SCTransform using Seurat and rpy2.")
    ro.r('''mtx = Matrix(data_mat, sparse = TRUE)
            rownames(mtx) = gene_names
            colnames(mtx) = cell_names
            sobj = CreateSeuratObject(counts = mtx, assay = "RNA")
            sobj = SCTransform(sobj, verbose = FALSE, method = "glmGamPoi")
            counts = GetAssayData(object = sobj, assay = "SCT", slot = "counts")
            data = GetAssayData(object = sobj, assay = "SCT", slot = "data")
            scale_data = GetAssayData(object = sobj, assay = "SCT", slot = "scale.data")
            variable_genes = VariableFeatures(object = sobj)
        ''')
    
    # Get results into Python
    sct_counts = ro.globalenv["counts"]
    sct_data = ro.globalenv["data"]
    sct_scale_data = ro.globalenv["scale_data"]
    sct_variable_genes = ro.globalenv["variable_genes"]

    # Reformat the scale data
    scale_data_df = pd.DataFrame(data=sct_scale_data.T, index=cell_names, columns=sct_variable_genes)

    # Add the sctransform normalization to the adata object
    adata.obsm["sctransform_scale_data"] = scale_data_df
    if filter_genes:
        adata.layers["sctransform_corrected_counts"] = sct_counts.T
        adata.layers["sctransform_corrected_log1p_counts"] = sct_data.T
    else:
        logger.info("We can't add the corrected counts to the adata object if we didn't filter because of a shape mismatch")
    adata.var["sctransform_genes"] = adata.var.index.isin(sct_variable_genes)


def plot_sctransform(
    adata: AnnData,
    outdir_path: str,
):
    logger.info(f"Plotting SCTransform normalized data against raw data.")
    plot_against_raw(adata, outdir_path, obsm_key="sctransform_scale_data")


def save_sctransform(
    adata: AnnData,
    outdir_path: str,
):
    logger.info(f"Saving SCTransform normalized data to {os.path.join(outdir_path, 'sctransform')}")
    X = adata.obsm["sctransform_scale_data"]
    X.to_csv(os.path.join(outdir_path, "sctransform_scale_data.csv"))    
    X.index.to_series().to_csv(os.path.join(outdir_path, "barcodes.tsv"), sep="\t", index=False, header=False)
    X.columns.to_series().to_csv(os.path.join(outdir_path, "features.tsv"), sep="\t", index=False, header=False)


def sctransform_recipe(
    adata: AnnData,
    outdir_path: str,
    filter_genes: bool = False,
    save_normalized_mtx: bool = False,
):
    run_sctransform(adata, filter_genes=filter_genes)
    plot_sctransform(adata, outdir_path)
    if save_normalized_mtx:
        if not os.path.exists(os.path.join(outdir_path, "sctransform")):
            os.makedirs(os.path.join(outdir_path, "sctransform"))
        save_sctransform(adata, os.path.join(outdir_path, "sctransform"))
