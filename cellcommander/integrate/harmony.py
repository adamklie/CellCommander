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
from typing import Iterable, Optional, Tuple, Union

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

logger = logging.getLogger("cellcommander")

def run_harmonyR(
    adata: AnnData,
    correction_key: Union[str, Iterable[str]],
    counts_key: Optional[str] = None,
    data_key: Optional[str] = None,
    scale_data_key: Optional[str] = None,
    sctransform_var_genes_column: Optional[str] = "sctransform_genes",
    n_comps: int = 50,
    random_state: int = 1234,
):
    # R imports
    Seurat = importr("Seurat")
    SeuratObject = importr("SeuratObject")
    Matrix = importr("Matrix")

    # Grab out of adata
    if counts_key is None:
        counts_data = adata.X.T.astype("float32")
    else:
        if counts_key in adata.layers.keys():
            counts_data = adata.layers[counts_key].T.astype("float32")
        elif counts_key in adata.obsm.keys():
            counts_data = adata.obsm[counts_key].values.T.astype("float32")
        else:
            raise ValueError("counts_key not found in adata layers or obsm.")
    if data_key is None:
        data = adata.X.T.astype("float32")
    else:
        if data_key in adata.layers.keys():
            data = adata.layers[data_key].T.astype("float32")
        elif data_key in adata.obsm.keys():
            data = adata.obsm[data_key].values.T.astype("float32")
        else:
            raise ValueError("data_key not found in adata layers or obsm.")
    if scale_data_key is None:
        scale_data = adata.X.T.astype("float32")
    else:
        if scale_data_key in adata.layers.keys():
            scale_data = adata.layers[scale_data_key].T.astype("float32")
        elif scale_data_key in adata.obsm.keys():
            scale_data = adata.obsm[scale_data_key].values.T.astype("float32")
        else:
            raise ValueError("scale_data_key not found in adata layers or obsm.")
    cell_names = adata.obs_names
    gene_names = adata.var_names
    sct_gene_names = adata.var[adata.var[sctransform_var_genes_column]].index.values

    # Toss into R's global environment
    ro.globalenv["counts_data"] = counts_data
    ro.globalenv["data"] = data
    ro.globalenv["scale_data"] = scale_data
    ro.globalenv["cell_names"] = cell_names
    ro.globalenv["gene_names"] = gene_names
    ro.globalenv["sct_gene_names"] = sct_gene_names
    ro.globalenv["random_state"] = random_state
    ro.globalenv["pca_n_comps"] = n_comps
    ro.globalenv["batch_key"] = correction_key

     # Run 
    ro.r('''set.seed(random_state)
            counts_mtx = Matrix(counts_data, sparse = TRUE)
            rownames(counts_mtx) = gene_names
            colnames(counts_mtx) = cell_names
            mtx = Matrix(data, sparse = TRUE)
            rownames(mtx) = gene_names
            colnames(mtx) = cell_names
            rownames(scale_data) = sct_gene_names
            colnames(scale_data) = cell_names
            sobj = CreateSeuratObject(assay="SCT", counts=counts_mtx, data=mtx, min.cells=0, min.features=0)
            sobj@assays$SCT@scale.data <- scale_data
            VariableFeatures(sobj[["SCT"]]) <- sct_gene_names
            sobj = RunPCA(sobj, verbose = FALSE)
            sobj <- RunHarmony(
                object = sobj,
                group.by.vars = batch_key,
                assay.use = 'SCT',
                max.iter.harmony=25,
                project.dim = FALSE
            )
            sobj <- RunUMAP(sobj, reduction.use = "harmony", dims.use = 1:pca_n_comps)
            pca_embedding = Embeddings(sobj, reduction = "pca")
            harmony_embedding = Embeddings(sobj, reduction = "harmony")
            umap_embedding = Embeddings(sobj, reduction = "umap")
            ''')
    
    # Grab the embedding to use
    pca_embedding = ro.globalenv["pca_embedding"]
    harmony_embedding = ro.globalenv["harmony_embedding"]
    umap_embedding = ro.globalenv["umap_embedding"]
    adata.obsm["X_pca"] = pca_embedding
    adata.obsm["X_harmonyR_pca"] = harmony_embedding
    adata.obsm["X_umap"] = umap_embedding
