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
from cellcommander.normalize.sctransform import run_sctransform
from cellcommander.reduce_dimensions.seurat import run_seurat_default
from cellcommander.integrate.harmony import run_harmonyR

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
        adata.var_names_make_unique()
        describe_anndata(adata)

        # Integrate
        if args.correction_method == "harmonyR":
            logger.info("Running harmonyR correction.")
            run_harmonyR(
                adata=adata,
                correction_key=args.batch_key,
                counts_key=None,
                data_key=None,
                scale_data_key=args.scale_data_key,
                n_comps=args.n_components,
                random_state=args.random_state,
            )
            use_rep = "X_pca_harmonyR"
            
        elif args.correction_method == "none":
            logger.info("No correction method specified. Using 'X_pca' in '.obsm' calculated on normalized data.")
            use_rep = args.obsm_key

        # Find neighbors on this
        sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=args.n_neighbors, random_state=args.random_state, n_pcs=args.n_components)
        
        # Run UMAP
        if "X_umap" not in adata.obsm.keys():
            logger.info(f"No UMAP coordinates found, running UMAP with {args.umap_min_distance} min_dist using ScanPy")
            sc.tl.umap(adata, min_dist=args.umap_min_distance, random_state=args.random_state)

        # Cluster unintegrated data
        sc.tl.leiden(adata, key_added=f"{args.normalization}_{args.correction_method}_leiden_{args.clust_resolution}", resolution=args.clust_resolution, random_state=args.random_state)

        # Plot the umap with
        plot_keys = [f"{args.normalization}_{args.correction_method}_leiden_{args.clust_resolution}"] + args.plot_keys + [args.names_key]
        logger.info(f"Plotting with keys: {plot_keys}")
        with plt.rc_context():
            sc.pl.umap(adata, color=plot_keys)
            plt.savefig(os.path.join(args.output_dir, f"{args.normalization}_{args.correction_method}_leiden_{args.clust_resolution}_umap.png"), bbox_inches='tight')
            plt.close()

        # Generate barplots of metadata in each cluster
        for key in plot_keys:
            cross_tab = pd.crosstab(adata.obs[f"{args.normalization}_{args.correction_method}_leiden_{args.clust_resolution}"], adata.obs[key])
            cross_tab_norm = cross_tab.div(cross_tab.sum(axis=1), axis=0)

            # Make a pretty horizontal bar plot with x-axis showing sample and y-axis showing manual_cellid_annotation
            cross_tab_norm.plot.barh(stacked=True, figsize=(20, 10))

            # Put the legend outside the plot
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

            # Save the figure
            plt.savefig(os.path.join(args.output_dir, f"{key}_barplot.png"), bbox_inches='tight')
            plt.close()

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
