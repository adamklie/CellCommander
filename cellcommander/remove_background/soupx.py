import logging
import os

import pickle
import matplotlib
import numpy as np
import pandas as pd
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


def run_soupx(
    adata: AnnData,
    adata_raw: AnnData,
    soupx_markers: pd.DataFrame,
    layer: Optional[str] = "soupx_counts",
    clust_num_hvgs: int = 3000,
    clust_n_neighbors: int = 30,
    clust_n_components: int = 50,
    clust_resolution: float = 1,
    umap_min_distance: float = 0.3,
    random_state: int = 1234,
    outdir_path: str = ".",
):
    # Run initial clustering
    logger.info("Performing initial clustering for SoupX.")
    adata_pp = adata.copy()
    
    # Normalize and log transform
    logger.info("Normalizing and log transforming with ScanPy defaults.")
    sc.pp.filter_genes(adata_pp, min_cells=20)
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    
    # Highly variable genes
    logger.info(f"Selecting {clust_num_hvgs} highly variable genes with ScanPy default.")
    sc.pp.highly_variable_genes(adata_pp, n_top_genes=clust_num_hvgs)
    adata_pp = adata_pp[:, adata_pp.var.highly_variable]
    
    # Dimensionality reduction
    logger.info(f"Running ScanPy dimensionality reduction with PCA "
                f"and ScanPy clustering with {clust_n_neighbors} neighbors, {clust_n_components} components, "
                f"and resolution {clust_resolution}.")
    sc.pp.pca(adata_pp, n_comps=clust_n_components, random_state=random_state)
    sc.pp.neighbors(adata_pp, n_neighbors=clust_n_neighbors, n_pcs=clust_n_components, random_state=random_state)
    sc.tl.leiden(adata_pp, key_added=f"pre_soupx_leiden_{clust_resolution}", resolution=clust_resolution, random_state=random_state)

    # UMAP
    logger.info(f"Running ScanPy UMAP with min_dist {umap_min_distance}.")
    sc.tl.umap(adata_pp, min_dist=umap_min_distance, random_state=random_state)
    with plt.rc_context():
        sc.pl.umap(adata_pp, color=[f"pre_soupx_leiden_{clust_resolution}"])
        plt.savefig(os.path.join(outdir_path, f"pre_soupx_leiden_{clust_resolution}_umap.png"))
        plt.close()

    # Add variables to original AnnData
    adata.obs[f"pre_soupx_leiden_{clust_resolution}"] = adata_pp.obs[f"pre_soupx_leiden_{clust_resolution}"]

    # Import libraries
    logger.info("Importing SoupX library.")
    SoupX = importr("SoupX")
    ro.r(f"set.seed({random_state})")
    
    # Preprocess variables for SoupX
    logger.info("Preprocessing variables for SoupX.")
    soupx_groups = adata_pp.obs[f"pre_soupx_leiden_{clust_resolution}"]
    metadata = pd.DataFrame(adata_pp.obsm["X_umap"])
    metadata.columns = ["RD1", "RD2"]
    metadata["Cluster"] = soupx_groups.values
    metadata["Annotation"] = soupx_groups.values
    metadata.index = adata_pp.obs.index
    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T.astype("float32")
    data_tod = adata_raw.X.T.astype("float32")
    raw_cells = adata_raw.obs_names
    raw_genes = adata_raw.var_names

    # Memory management
    del adata_pp, adata_raw

    # Pass in soupx_markers as R vector of characters
    keys = soupx_markers["Group"].values
    values = soupx_markers["Gene"].values
    soupx_markers_dict = dict(zip(keys, values))
    soupx_markers_r = ro.ListVector(soupx_markers_dict)

    # Set up R objects
    logger.info("Setting up R objects.")
    ro.globalenv["soupx_markers"] = soupx_markers_r
    ro.globalenv["data"] = data
    ro.globalenv["data_tod"] = data_tod
    ro.globalenv["genes"] = genes
    ro.globalenv["cells"] = cells
    ro.globalenv["raw_genes"] = raw_genes
    ro.globalenv["raw_cells"] = raw_cells
    ro.globalenv["metadata"] = metadata
    ro.globalenv["outdir_path"] = outdir_path

    # Run SoupX!
    logger.info(f"Running SoupX and saving RDS file to {os.path.join(outdir_path, 'soupx.rds')}.")
    ro.r('''
        nonExpressedGeneList = soupx_markers
        data <- as(data, "dgCMatrix")
        data_tod <- as(data_tod, "dgCMatrix")
        rownames(data) <- genes
        colnames(data) <- cells
        rownames(data_tod) <- raw_genes
        colnames(data_tod) <- raw_cells
        soupc = SoupChannel(tod=data_tod, toc=data)
        soupc = setDR(soupc, metadata[colnames(soupc$toc), c("RD1", "RD2")])
        soupc = setClusters(soupc, setNames(metadata$Cluster, rownames(metadata)))
        useToEst = estimateNonExpressingCells(soupc, nonExpressedGeneList = nonExpressedGeneList)
        soupc = calculateContaminationFraction(soupc, nonExpressedGeneList, useToEst=useToEst)
        contam_frac = 100*exp(coef(soupc$fit))[[1]]
        out = adjustCounts(soupc, roundToInt = TRUE)
        saveRDS(soupc, file = file.path(outdir_path, "soupx.rds"))
    ''')
    
    # Grab the output
    logger.info("Cleaning up the output and generating some summary statistics.")
    contam_frac_py = ro.globalenv["contam_frac"]
    out_py = ro.globalenv["out"]

    # Manage memory
    del data
    del data_tod

    # Save output
    soupx_out_dict = {}
    soupx_out_dict["gene_markers_used"] = soupx_markers["Gene"].values
    soupx_out_dict["cluster_groups_for_markers"] = soupx_markers["Group"].values
    soupx_out_dict["soup_contamination_fraction"] = np.round(contam_frac_py[0], 2)

    # Add in the counts adjustments to the gene metadata
    cnt_soggy = adata.X.sum(axis=0).A.squeeze()
    cnt_strained = out_py.T.sum(axis=0).A.squeeze()
    num_zeroed_cells = (cnt_soggy - cnt_strained)
    adata.var["num_soupx_zeroed_cells"] = num_zeroed_cells
    
    # Total counts removed
    soupx_out_dict["soup_total_counts_removed"] = np.round(cnt_soggy.sum() - cnt_strained.sum(), 2)

    # Top 10 genes with the most zeroed cells
    soupx_out_dict["top_10_genes"] = adata.var.sort_values("num_soupx_zeroed_cells", ascending=False).head(10).index.tolist()

    # Add in the new counts
    logger.info(f"Adding SoupX counts to AnnData, will be in '.X' and '{layer}' layers.")
    adata.layers["counts"] = adata.X
    adata.layers[layer] = out_py.T
    adata.X = adata.layers[layer]

    # Reprocess after SoupX
    logger.info("Reprocessing after SoupX to generate clustering and UMAP with ScanPy. Useful to compare to pre-SoupX.")
    adata_pp = adata.copy()
    adata_pp.X = adata.layers[layer]
    sc.pp.filter_genes(adata_pp, min_cells=20)
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.highly_variable_genes(adata_pp, n_top_genes=clust_num_hvgs)
    adata_pp = adata_pp[:, adata_pp.var.highly_variable]
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp, n_neighbors=clust_n_neighbors, n_pcs=clust_n_components, random_state=random_state)
    sc.tl.leiden(adata_pp, key_added=f"post_soupx_leiden_{clust_resolution}", resolution=clust_resolution, random_state=random_state)
    sc.tl.umap(adata_pp, min_dist=umap_min_distance, random_state=random_state)
    with plt.rc_context():
        sc.pl.umap(adata_pp, color=[f"post_soupx_leiden_{clust_resolution}"])
        plt.savefig(os.path.join(outdir_path, f"post_soupx_leiden_{clust_resolution}_umap.png"))
        plt.close()

    # Add variables to original AnnData
    adata.obs[f"post_soupx_leiden_{clust_resolution}"] = adata_pp.obs[f"post_soupx_leiden_{clust_resolution}"]

    # Save the dict
    logger.info(f"Dumping SoupX stats to {os.path.join(outdir_path, 'soupx_stats.pickle')}.")
    pickle.dump(soupx_out_dict, open(os.path.join(outdir_path, f"soupx_stats.pickle"), "wb"))
