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
anndata2ri.activate()
anndata2ri.scipy2ri.activate()
matplotlib.use("Agg")
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

logger = logging.getLogger("cellcommander")


def run_seurat_default(
    adata: AnnData,
    counts_key: Optional[str] = None,
    data_key: Optional[str] = None,
    scale_data_key: Optional[str] = None,
    var_genes_column: Optional[str] = "highly_variable",
    n_comps: int = 50,
    random_state: Optional[int] = 1234,
):
    # R imports
    logger.info("Importing necessary R modules for running Seurat.")
    Seurat = importr("Seurat")
    SeuratObject = importr("SeuratObject")
    Matrix = importr("Matrix")

    # Grab out of adata
    logger.info("Preparing data for Seurat.")
    if counts_key is None:
        counts = adata.X.T.astype("float32")
    else:
        if counts_key in adata.layers.keys():
            counts = adata.layers[counts_key].T.astype("float32")
        elif counts_key in adata.obsm.keys():
            counts = adata.obsm[counts_key].values.T.astype("float32")
        else:
            raise ValueError("counts_key not found in adata layers or obsm.")
    logger.info(f"Counts matrix shape: {counts.shape}")
    if data_key is None:
        data = adata.X.T.astype("float32")
    else:
        if data_key in adata.layers.keys():
            data = adata.layers[data_key].T.astype("float32")
        elif data_key in adata.obsm.keys():
            data = adata.obsm[data_key].values.T.astype("float32")
        else:
            raise ValueError("data_key not found in adata layers or obsm.")
    logger.info(f"Data matrix shape: {data.shape}")
    if scale_data_key is None:
        scale_data = adata.X.T.astype("float32")
    else:
        if scale_data_key in adata.layers.keys():
            scale_data = adata.layers[scale_data_key].T.astype("float32")
        elif scale_data_key in adata.obsm.keys():
            scale_data = adata.obsm[scale_data_key].values.T.astype("float32")
        else:
            raise ValueError("scale_data_key not found in adata layers or obsm.")
    logger.info(f"Scale data matrix shape: {scale_data.shape}")
    cell_names = adata.obs_names
    gene_names = adata.var_names
    var_gene_names = adata.var[adata.var[var_genes_column]].index.values
    logger.info(f"Number of variable genes: {len(var_gene_names)}")

    # Toss into R's global environment
    ro.globalenv["counts"] = counts
    ro.globalenv["data"] = data
    ro.globalenv["scale_data"] = scale_data
    ro.globalenv["cell_names"] = cell_names
    ro.globalenv["gene_names"] = gene_names
    ro.globalenv["var_gene_names"] = var_gene_names
    ro.globalenv["random_state"] = random_state
    ro.globalenv["pca_n_comps"] = n_comps

    # Run 
    logger.info("Running Seurat dimensionality reduction.")
    ro.r('''set.seed(random_state)
            counts_mtx = Matrix(counts, sparse = TRUE)
            rownames(counts_mtx) = gene_names
            colnames(counts_mtx) = cell_names
            mtx = Matrix(data, sparse = TRUE)
            rownames(mtx) = gene_names
            colnames(mtx) = cell_names
            rownames(scale_data) = var_gene_names
            colnames(scale_data) = cell_names
            sobj = CreateSeuratObject(assay="RNA", counts=counts_mtx, data=mtx, min.cells=0, min.features=0)
            sobj@assays$RNA@scale.data <- scale_data
            VariableFeatures(sobj[["RNA"]]) <- var_gene_names
            sobj = RunPCA(sobj, verbose = FALSE)
            sobj = RunUMAP(sobj, dims=1:pca_n_comps, reduction.name='umap.rna', reduction.key='rnaUMAP_', seed.use=random_state)
            pca_rna = Embeddings(sobj, reduction = "pca")
            umap_rna = Embeddings(sobj, reduction = "umap.rna")
            ''')
    
    # Grab the embedding to use
    pca_rna = ro.globalenv["pca_rna"]
    umap_rna = ro.globalenv["umap_rna"]
    adata.obsm["X_pca"] = pca_rna
    adata.obsm["X_umap"] = umap_rna
