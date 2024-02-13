import logging
import os

import matplotlib
import scanpy as sc
import anndata2ri
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr
from scipy.sparse import coo_matrix
from anndata import AnnData
from mudata import MuData
anndata2ri.activate()
anndata2ri.scipy2ri.activate()
matplotlib.use("Agg")
from cellcommander.joint_integrate import consts
from typing import Iterable, Optional, Tuple

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')\

logger = logging.getLogger("cellcommander")


def run_seurat_wnn(
    mdata: MuData,
    rna_mod_key: str = "rna",
    atac_mod_key: str = "atac",
    rna_obsm_key: str = "X_pca",
    atac_obsm_key: str = "X_lsi",
    rna_dim_start: int = 1,
    rna_dim_end: int = 50,
    atac_dim_start: int = 2,
    atac_dim_end: int = 50,
    cluster_key=None,
    cluster_resolution=1,
    random_state=1234,
):
    logger.info("Importing necessary R libraries for WNN")
    Seurat = importr("Seurat")
    SeuratObject = importr("SeuratObject")
    Matrix = importr("Matrix")
    Signac = importr("Signac")

    # Grab modalities from MuData
    rna = mdata.mod[rna_mod_key]
    atac = mdata.mod[atac_mod_key]

    # Prepare rna data for Seurat
    logger.info("Preparing RNA data for Seurat")
    rna_data_mat = rna.X.T.astype("float32")
    rna_cell_names = rna.obs_names
    rna_gene_names = rna.var_names
    rna_dim = rna.obsm[rna_obsm_key]
    rna_dim_names = [f"rna_dim_{i}" for i in range(1, 1 + rna.obsm[rna_obsm_key].shape[1])]
    ro.globalenv["rna_data_mat"] = rna_data_mat
    ro.globalenv["rna_cell_names"] = rna_cell_names
    ro.globalenv["rna_gene_names"] = rna_gene_names
    ro.globalenv["rna_dim"] = rna_dim
    ro.globalenv["rna_dim_names"] = rna_dim_names

    # Prepare atac data for Seurat
    logger.info("Preparing ATAC data for Seurat")
    atac_data_mat = atac.X.T.astype("float32")
    atac_cell_names = atac.obs_names
    atac_feature_names = atac.var_names
    atac_dim = atac.obsm[atac_obsm_key]
    atac_dim_names = [f"atac_dim_{i}" for i in range(1, 1 + atac.obsm[atac_obsm_key].shape[1])]
    ro.globalenv["atac_data_mat"] = atac_data_mat
    ro.globalenv["atac_cell_names"] = atac_cell_names
    ro.globalenv["atac_feature_names"] = atac_feature_names
    ro.globalenv["atac_dim"] = atac_dim
    ro.globalenv["atac_dim_names"] = atac_dim_names

    # Others
    ro.globalenv["rna_dim_start"] = rna_dim_start
    ro.globalenv["rna_dim_end"] = rna_dim_end if rna_dim_end < rna_dim.shape[1] else rna_dim.shape[1]
    ro.globalenv["atac_dim_start"] = atac_dim_start
    ro.globalenv["atac_dim_end"] = atac_dim_end if atac_dim_end < atac_dim.shape[1] else atac_dim.shape[1]
    ro.globalenv["cluster_resolution"] = cluster_resolution
    ro.globalenv["seed"] = random_state

    # Run WNN
    logger.info("Running WNN")
    ro.r('''
        # Set seed
        set.seed(seed)
        # Create RNA Seurat object
        rna_mtx = Matrix(rna_data_mat, sparse = TRUE)
        rownames(rna_mtx) = rna_gene_names
        colnames(rna_mtx) = rna_cell_names
        rna_sobj = CreateSeuratObject(counts = rna_mtx, assay = "RNA")
        rownames(rna_dim) = rna_cell_names
        colnames(rna_dim) = rna_dim_names
        rna_reduction <- CreateDimReducObject(embeddings = rna_dim, key = "rna_dim", assay = "RNA")
        rna_sobj[["rna_dim"]] <- rna_reduction

        # Create ATAC Seurat object
        atac_mtx = Matrix(atac_data_mat, sparse = TRUE)
        rownames(atac_mtx) = atac_feature_names
        colnames(atac_mtx) = atac_cell_names
        atac_sobj = CreateSeuratObject(counts = atac_mtx, assay = "ATAC")
        rownames(atac_dim) = atac_cell_names
        colnames(atac_dim) = atac_dim_names
        atac_reduction <- CreateDimReducObject(embeddings = atac_dim, key = "atac_dim", assay = "ATAC")
        atac_sobj[["atac_dim"]] <- atac_reduction
        
        # Create multiome object 
        multiome <- rna_sobj
        multiome[["ATAC"]] <- CreateAssayObject(counts = atac_sobj@assays$ATAC@counts)
        multiome@reductions$atac_dim <- atac_sobj@reductions$atac_dim
            
        # Run WNN
        multiome <- FindMultiModalNeighbors(multiome, reduction.list=list('rna_dim', 'atac_dim'), dims.list=list(c(rna_dim_start:rna_dim_end), c(atac_dim_start:atac_dim_end)), verbose=TRUE)
        wnn <- as.data.frame(summary(multiome@graphs$wknn))
        RNA_weight = multiome$RNA.weight
        ATAC_weight = multiome$ATAC.weight
        ''')
    
    # Save WNN results
    wnn = ro.globalenv["wnn"]
    wnn["i"] = wnn["i"] - 1
    wnn["j"] = wnn["j"] - 1
    mdata.obsp["wnn_connectivities"] = coo_matrix(
        (wnn["x"], (wnn["i"], wnn["j"]))
    )
    mdata.obs["wnn_RNA_weight"] = ro.globalenv["RNA_weight"]
    mdata.obs["wnn_ATAC_weight"] = ro.globalenv["ATAC_weight"]
    
    # Optionally add Seurat clusters and dimensionality reduction with UMAP
    if cluster_key:
        logger.info("Running Seurat UMAP and clustering using Leiden algorithm")
        ro.r('''
            multiome <- RunUMAP(multiome, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
            multiome <- FindClusters(multiome, graph.name='wsnn', algorithm=4, resolution = cluster_resolution, verbose=FALSE)
            umap_wnn = Embeddings(multiome, reduction = "umap.wnn")
            clusters = multiome$seurat_clusters 
            ''')
        mdata.obsm["X_umap_wnn"] = ro.globalenv["umap_wnn"]
        mdata.obs[cluster_key] = ro.globalenv["clusters"]
