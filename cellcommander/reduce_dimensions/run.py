"""Single run of reduce_dimensions, given input arguments."""

import argparse
import logging
import os
import sys
import traceback
from datetime import datetime
from typing import Dict, Optional, Tuple, Union

import numpy as np
import matplotlib
import muon as mu
import pandas as pd
import psutil
import scanpy as sc
import snapatac2 as snap
from mudata import MuData
from anndata import AnnData

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')
import seaborn as sns

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)
from cellcommander.utils import describe_anndata
from cellcommander.reduce_dimensions import consts
from cellcommander.reduce_dimensions.lsi import run_muon_lsi
from cellcommander.reduce_dimensions.scanpy import run_scanpy_default
from cellcommander.reduce_dimensions.seurat import run_seurat_default
from cellcommander.reduce_dimensions.spectral import run_snapatac2_spectral
logger = logging.getLogger("cellcommander")


def run_reduce_dimensions(args: argparse.Namespace):
    """The full script for the command line tool to perform reduce_dimensions and filtering.

    Args:
        args: Inputs from the command line, already parsed using argparse.

    Note: Returns nothing, but writes output to a file(s) specified from command line.

    """

    # Log the start time.
    logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Running reduce-dimensions command")

    try:
            
        # Read in single h5 file
        logger.info(f"Reading h5 file from {args.input_file}")
        adata = sc.read_h5ad(args.input_file)
        adata.var_names_make_unique()
        describe_anndata(adata)
        
        # Run methods
        if "lsi" == args.method:
            logger.info("Using Muon to run LSI on the data")
            # If variable features are passed in subset to those
            if args.variable_features_key is not None:
                logger.info(f"Subsetting to variable features in {args.variable_features_key}")
                adata_pp = adata[:, adata.var[args.variable_features_key]].copy()
            else:
                logger.info(f"No variable features passed in, using all features")
                adata_pp = adata.copy()

            # If layer is passed in, use that as the data
            if args.layer is not None:
                logger.info(f"Setting .X to come from {args.layer} for Muon since it doesn't allow to pass in a layer")
                adata_pp.X = adata_pp.layers[args.layer].copy()
            
            # If obsm is passed in, use that as the data
            elif args.obsm_key is not None:
                raise NotImplementedError("'obsm' functionality not implemented yet, use --layer argument instead")
            
            # Run the method
            run_muon_lsi(adata_pp, scale_embeddings=args.scale_data, n_comps=args.n_components)

            # Remove components if specified
            if args.components_to_remove is not None:
                logger.info(f"Removing components {args.components_to_remove} from LSI, {adata_pp.obsm['X_lsi'].shape[1]} components total to start")
                adata_pp.obsm['X_lsi'] = np.delete(adata_pp.obsm['X_lsi'], args.components_to_remove, axis=1)
                adata_pp.varm["LSI"] = np.delete(adata_pp.varm["LSI"], args.components_to_remove, axis=1)
                adata_pp.uns["lsi"]["stdev"] = np.delete(adata_pp.uns["lsi"]["stdev"], args.components_to_remove, axis=0)
                logger.info(f"Removed components {args.components_to_remove} from LSI, {adata_pp.obsm['X_lsi'].shape[1]} components total remaining")
                n_components = args.n_components - len(args.components_to_remove)
            else:
                n_components = args.n_components

            # Final dimensionality reduction
            adata_pp.obsm["X_reduced"] = adata_pp.obsm["X_lsi"].copy()

        elif "spectral" == args.method:
            # Handle variable features
            if args.variable_features_key is None:
                logger.info(f"No variable features passed in, using all features")
            else:
                logger.info(f"Passing in variable features to SnapATAC2 from {args.variable_features_key}")
            adata_pp = adata.copy()
            run_snapatac2_spectral(
                adata=adata_pp,
                features_key=args.variable_features_key,
                n_comps=args.n_components,
                random_state=args.random_state,
            )

            # Remove components if specified
            if args.components_to_remove is not None:
                logger.info(f"Removing components {args.components_to_remove} from spectral, {adata_pp.obsm['X_spectral'].shape[1]} components total to start")
                adata_pp.obsm['X_spectral'] = np.delete(adata_pp.obsm['X_spectral'], args.components_to_remove, axis=1)
                logger.info(f"Removed components {args.components_to_remove} from spectral, {adata_pp.obsm['X_spectral'].shape[1]} components total remaining")
                n_components = args.n_components - len(args.components_to_remove)
            else:
                n_components = args.n_components

            # Run UMAP
            logger.info(f"Running UMAP on spectral embeddings using SnapATAC2 implementation")
            snap.tl.umap(adata_pp, use_rep="X_spectral", random_state=args.random_state)

            # Final dimensionality reduction
            adata_pp.obsm["X_reduced"] = adata_pp.obsm["X_spectral"].copy()

        elif "scanpy_default" == args.method:
            logger.info("Using ScanPy to run PCA on the data")
            # If variable features are passed in subset to those
            if args.variable_features_key is not None:
                logger.info(f"Subsetting to variable features in {args.variable_features_key}")
                adata_pp = adata[:, adata.var[args.variable_features_key]].copy()
            else:
                logger.info(f"No variable features passed in, using all features")
                adata_pp = adata.copy()

            # If layer is passed in, use that as the data
            if args.layer is not None:
                logger.info(f"Setting .X to come from {args.layer} and running PCA on that layer (ScanPy default)")
                adata_pp.X = adata_pp.layers[args.layer].copy()
                run_scanpy_default(adata_pp, layer=args.layer, n_comps=args.n_components, random_state=args.random_state)

            elif args.obsm_key is not None:
                logger.info(f"Running PCA on obsm key {args.obsm_key} (ScanPy default).")
                run_scanpy_default(adata_pp, obsm_key=args.obsm_key, n_comps=args.n_components, random_state=args.random_state)
            else:
                raise NotImplementedError("Must specify either layer or obsm_key to run ScanPy default PCA.")
            
            # Remove components if specified
            if args.components_to_remove is not None:
                logger.info(f"Removing components {args.components_to_remove} from PCA, {adata_pp.obsm['X_pca'].shape[1]} components total to start")
                adata_pp.obsm['X_pca'] = np.delete(adata_pp.obsm['X_pca'], args.components_to_remove, axis=1)
                adata_pp.varm["PCs"] = np.delete(adata_pp.varm["PCs"], args.components_to_remove, axis=1)
                logger.info(f"Removed components {args.components_to_remove} from PCA, {adata_pp.obsm['X_pca'].shape[1]} components total remaining")
                n_components = args.n_components - len(args.components_to_remove)
            else:
                n_components = args.n_components

            # Final dimensionality reduction
            adata_pp.obsm["X_reduced"] = adata_pp.obsm["X_pca"].copy()

        elif "seurat_default" == args.method:
            logger.info("Using Seurat to run PCA and UMAP on the data")
            # Handle variable features
            if args.variable_features_key is None:
                logger.info(f"No variable features passed in, using all features")
            else:
                logger.info(f"Passing in variable features to Seurat from {args.variable_features_key}")
            adata_pp = adata.copy()

            # If layer is passed in, use that as the data
            if args.layer is not None:
                logger.info(f"Using layer {args.layer} as scaled data input, make sure you want to do this")
                run_seurat_default(
                    adata_pp, 
                    scale_data_key=args.layer,
                    var_genes_column=args.variable_features_key,
                    n_comps=args.n_components, 
                    random_state=args.random_state
                )
            elif args.obsm_key is not None:
                logger.info(f"Using obsm key {args.obsm_key} as scaled data input")
                run_seurat_default(
                    adata_pp, 
                    scale_data_key=args.obsm_key,
                    var_genes_column=args.variable_features_key,
                    n_comps=args.n_components, 
                    random_state=args.random_state
                )

            # Remove components if specified
            if args.components_to_remove is not None:
                logger.info(f"Removing components {args.components_to_remove} from PCA, {adata_pp.obsm['X_pca'].shape[1]} components total to start")
                adata_pp.obsm['X_pca'] = np.delete(adata_pp.obsm['X_pca'], args.components_to_remove, axis=1)
                adata_pp.varm["PCs"] = np.delete(adata_pp.varm["PCs"], args.components_to_remove, axis=1)
                logger.info(f"Removed components {args.components_to_remove} from PCA, {adata_pp.obsm['X_pca'].shape[1]} components total remaining")
                n_components = args.n_components - len(args.components_to_remove)
            else:
                n_components = args.n_components
            adata_pp.obsm["X_reduced"] = adata_pp.obsm["X_pca"].copy()
            
        # Get neighborhood graph
        logger.info(f"Creating kNN graph using {args.n_neighbors} neighbors and {n_components} components of the reduced data")
        sc.pp.neighbors(adata_pp, use_rep="X_reduced", n_neighbors=args.n_neighbors, n_pcs=n_components, random_state=args.random_state)

        # Run leiden clustering
        logger.info(f"Running leiden clustering with resolution {args.clust_resolution}")
        key_added = f"leiden_{args.clust_resolution}"
        sc.tl.leiden(adata_pp, resolution=args.clust_resolution, random_state=args.random_state, key_added=key_added)

        # Plot UMAP with leiden clusters
        if "X_umap" not in adata_pp.obsm.keys():
            logger.info(f"No UMAP coordinates found, running UMAP with {args.umap_min_distance} min_dist using ScanPy")
            sc.tl.umap(adata_pp, min_dist=args.umap_min_distance, random_state=args.random_state)
        with plt.rc_context({"figure.figsize": (6, 6)}):
            sc.pl.umap(adata_pp, color=[key_added], show=False, legend_loc="on data")
            plt.savefig(os.path.join(args.output_dir, f"{key_added}_umap.png"))
            plt.close()

        # Save results in old adata object
        adata.obsm[f"X_{args.method}"] = adata_pp.obsm["X_reduced"].copy()
        adata.obsm[f"X_{args.method}_umap"] = adata_pp.obsm["X_umap"].copy()
        adata.obs[key_added] = adata_pp.obs[key_added].copy()

        # Save a TSV
        adata.obs.to_csv(os.path.join(args.output_dir, "cell_metadata.tsv"), sep="\t")

        # Save the dim reductions as well
        pca_df = pd.DataFrame(adata.obsm[f"X_{args.method}"], index=adata.obs_names, columns=[f"{args.method}_{i}" for i in range(adata.obsm[f"X_{args.method}"].shape[1])])
        pca_df.to_csv(os.path.join(args.output_dir, f"{args.method}_dim_reductions.tsv"), sep="\t")
        umap_df = pd.DataFrame(adata.obsm[f"X_{args.method}_umap"], index=adata.obs_names, columns=[f"{args.method}_umap_{i}" for i in range(adata.obsm[f"X_{args.method}_umap"].shape[1])])
        umap_df.to_csv(os.path.join(args.output_dir, f"{args.method}_umap.tsv"), sep="\t")

        # Save the adata
        logger.info(
            f"Saving adata with dimensionality reduction and iniitial clustering to {os.path.join(args.output_dir, f'{args.output_prefix}.h5ad')}"
        )
        adata.write(os.path.join(args.output_dir, f"{args.output_prefix}.h5ad"))

        # Log the end time
        logger.info("Completed reduce-dimensions")
        logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))

    # The exception allows user to end inference prematurely with CTRL-C.
    except KeyboardInterrupt:
        # If partial output has been saved, delete it.
        full_file = args.output_dir + args.output_prefix + ".h5ad"

        # Name of the filtered (cells only) file.
        if os.path.exists(full_file):
            os.remove(full_file)

        logger.info("Keyboard interrupt.  Terminated without saving\n")
