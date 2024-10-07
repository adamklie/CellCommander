"""Single run of integrate, given input arguments."""

import argparse
import logging
import os
import sys
import yaml
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
import anndata as ad

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
from cellcommander.integrate import consts
from cellcommander.integrate.harmony import run_harmonyR, run_harmony

logger = logging.getLogger("cellcommander")


def run_integrate(args: argparse.Namespace):
    """The full script for the command line tool to perform integrate and filtering.

    Args:
        args: Inputs from the command line, already parsed using argparse.

    Note: Returns nothing, but writes output to a file(s) specified from command line.

    """

    # Log the start time.
    logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Running integrate command")

    try:
            
        # Read in single h5 file
        logger.info(f"Reading h5 file from {args.input_file}")
        adata = sc.read_h5ad(args.input_file)
        bcs = adata.obs.index
        adata.var_names_make_unique()
        describe_anndata(adata)

        # Merge batch_file info
        if args.batch_file:
            batch_info = pd.read_csv(args.batch_file, sep="\t", index_col=0)
            adata_obs = adata.obs.merge(batch_info, left_index=True, right_index=True, how="left")
            adata_obs.index = bcs
            adata.obs = adata_obs

        # Make a copy
        adata_pp = adata.copy()

        # Remove components prior to correction if specified
        if args.components_to_remove is not None:
            logger.info(f"Removing components {args.components_to_remove} from spectral, {adata_pp.obsm[args.obsm_key].shape[1]} components total to start")
            adata_pp.obsm[args.obsm_key] = np.delete(adata_pp.obsm[args.obsm_key], args.components_to_remove, axis=1)
            logger.info(f"Removed components {args.components_to_remove} from spectral, {adata_pp.obsm[args.obsm_key].shape[1]} components total remaining")
            n_components = args.n_components - len(args.components_to_remove)
        else:
            n_components = args.n_components

        # Integrate
        if args.method == "harmonyR":
            logger.info("Running harmonyR correction.")
            run_harmonyR(
                adata=adata_pp, 
                obsm_key=args.obsm_key,
                corrected_obsm_key=args.corrected_obsm_key,
                vars_to_correct=args.vars_to_correct,
                max_iter_harmony=args.max_iter_harmony,
                random_state=args.random_state
            )
        elif args.method == "harmony":
            logger.info("Running harmony correction.")
            run_harmony(
                adata=adata_pp, 
                obsm_key=args.obsm_key,
                corrected_obsm_key=args.corrected_obsm_key,
                vars_to_correct=args.vars_to_correct,
                max_iter_harmony=args.max_iter_harmony,
                random_state=args.random_state
            )

        # Run kNN
        logger.info(f"Creating kNN graph using {args.n_neighbors} neighbors and {n_components} components of the reduced data")
        sc.pp.neighbors(adata_pp, use_rep=args.corrected_obsm_key, n_neighbors=args.n_neighbors, random_state=args.random_state, n_pcs=args.n_components)
        
        # Run UMAP
        logger.info(f"No UMAP coordinates found, running UMAP with {args.umap_min_distance} min_dist using ScanPy")
        sc.tl.umap(adata_pp, min_dist=args.umap_min_distance, random_state=args.random_state)

        # Cluster unintegrated data
        logger.info(f"Running leiden clustering with resolution {args.clust_resolution}")
        key_added = f"{args.corrected_obsm_key}_leiden_{args.clust_resolution}"
        sc.tl.leiden(adata_pp, resolution=args.clust_resolution, random_state=args.random_state, key_added=key_added)

        # Plot the umap with
        ncol = 3
        nrow = 1
        figsize = 4
        wspace = 0.5
        fig, axs = plt.subplots(
            nrow, ncol, figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
        )
        plt.subplots_adjust(wspace=wspace)
        sc.pl.embedding(adata_pp, basis="X_seurat_default_umap", color='sample', s=10, ax=axs[0], show=False)
        sc.pl.embedding(adata_pp, basis="X_umap", color='sample', s=10, ax=axs[1], show=False)
        sc.pl.embedding(adata_pp, basis="X_umap", color=f"{args.corrected_obsm_key}_leiden_{args.clust_resolution}", s=10, ax=axs[2], show=False)
        fig.tight_layout()
        plt.savefig(os.path.join(args.output_dir, f"{key_added}_umap.png"))
        plt.close()
        
        # Generate barplots of metadata in each cluster
        for var_to_correct in args.vars_to_correct:
            cross_tab = pd.crosstab(adata_pp.obs[f"{args.corrected_obsm_key}_leiden_{args.clust_resolution}"], adata_pp.obs[var_to_correct])
            cross_tab_norm = cross_tab.div(cross_tab.sum(axis=1), axis=0)

            # Make a pretty horizontal bar plot with x-axis showing sample and y-axis showing manual_cellid_annotation
            cross_tab_norm.plot.barh(stacked=True, figsize=(4, 10))

            # Put the legend outside the plot
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

            # Save the figure
            plt.savefig(os.path.join(args.output_dir, f"{key_added}_{var_to_correct}_barplot.png"), bbox_inches='tight')
            plt.close()

        # Save results in old adata object
        adata.obsm[args.corrected_obsm_key] = adata_pp.obsm[args.corrected_obsm_key]
        adata.obsm[args.corrected_obsm_key + "_umap"] = adata_pp.obsm["X_umap"]
        adata.obs[f"{args.corrected_obsm_key}_leiden_{args.clust_resolution}"] = adata_pp.obs[f"{args.corrected_obsm_key}_leiden_{args.clust_resolution}"]

        # Save the adata
        logger.info(
            f"Saving integrated adata to {os.path.join(args.output_dir, f'{args.output_prefix}.h5ad')}"
        )
        adata.write(os.path.join(args.output_dir, f"{args.output_prefix}.h5ad"))

        # Log the end time
        logger.info("Completed integrate")
        logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))

    # The exception allows user to end inference prematurely with CTRL-C.
    except KeyboardInterrupt:
        # If partial output has been saved, delete it.
        full_file = args.output_dir + args.output_prefix + ".h5ad"

        # Name of the filtered (cells only) file.
        if os.path.exists(full_file):
            os.remove(full_file)

        logger.info("Keyboard interrupt. Terminated without saving\n")

