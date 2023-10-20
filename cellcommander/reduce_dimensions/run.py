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

        # Make a copy of the AnnData object to run everything on
        

        # If variable features are passed in subset to those
        if args.variable_features_key is not None:
            logger.info(f"Subsetting to variable features in {args.variable_features_key}")
            adata_pp = adata[:, adata.var[args.variable_features_key]].copy()
        else:
            logger.info(f"No variable features passed in, using all features")
            adata_pp = adata.copy()
        
        # Run methods
        if "lsi" == args.method:
            if args.layer is not None:
                logger.info(f"Setting .X to come from {args.layer}")
                adata_pp.X = adata_pp.layers[args.layer]
            elif args.obsm is not None:
                raise NotImplementedError("obsm functionality not implemented yet, use --layer argument instead")
            run_muon_lsi(adata_pp, scale_embeddings=args.scale_data, n_comps=args.n_components)
            if args.components_to_remove is not None:
                logger.info(f"Removing components {args.components_to_remove} from LSI, {adata_pp.obsm['X_lsi'].shape[1]} components total to start")
                adata_pp.obsm['X_lsi'] = np.delete(adata_pp.obsm['X_lsi'], args.components_to_remove, axis=1)
                adata_pp.varm["LSI"] = np.delete(adata_pp.varm["LSI"], args.components_to_remove, axis=1)
                adata_pp.uns["lsi"]["stdev"] = np.delete(adata_pp.uns["lsi"]["stdev"], args.components_to_remove, axis=0)
                logger.info(f"Removed components {args.components_to_remove} from LSI, {adata_pp.obsm['X_lsi'].shape[1]} components total remaining")
                n_components = args.n_components - len(args.components_to_remove)
            else:
                n_components = args.n_components
            adata_pp.obsm["X_reduced"] = adata_pp.obsm["X_lsi"]

        # Run UMAP
        logger.info(f"Creating kNN graph using {args.n_neighbors} neighbors and {n_components} components of the reduced data")
        sc.pp.neighbors(adata_pp, use_rep="X_reduced", n_neighbors=args.n_neighbors, n_pcs=n_components, random_state=args.random_state)
        sc.tl.umap(adata_pp, min_dist=args.umap_min_distance, random_state=args.random_state)
        key_added = f"initial_leiden_{args.initial_clust_resolution}"
        sc.tl.leiden(adata_pp, resolution=args.initial_clust_resolution, random_state=args.random_state, key_added=key_added)

        # Plot UMAP with leiden clusters
        with plt.rc_context({"figure.figsize": (6, 6)}):
            sc.pl.umap(adata_pp, color=[key_added], show=False)
            plt.savefig(os.path.join(args.output_dir, f"{key_added}_umap.png"))
            plt.close()

        # Save results in old adata object
        adata.obsm[f"X_{args.method}"] = adata_pp.obsm["X_reduced"]
        adata.obsm[f"X_{args.method}_umap"] = adata_pp.obsm["X_umap"]
        adata.obs[key_added] = adata_pp.obs[key_added]

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
