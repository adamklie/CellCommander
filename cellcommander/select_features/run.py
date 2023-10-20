"""Single run of select_features, given input arguments."""

import argparse
import logging
import os
import sys
import traceback
from datetime import datetime
from typing import Dict, Optional, Tuple, Union

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
from cellcommander.select_features import consts
from cellcommander.select_features.scanpy import run_seurat, run_cell_ranger, run_seurat_v3
from cellcommander.select_features.deviance import run_deviance
from cellcommander.select_features.signac import run_signac
from cellcommander.select_features.utils import is_counts, is_mostly_counts

logger = logging.getLogger("cellcommander")


def run_select_features(args: argparse.Namespace):
    """The full script for the command line tool to perform select_features and filtering.

    Args:
        args: Inputs from the command line, already parsed using argparse.

    Note: Returns nothing, but writes output to a file(s) specified from command line.

    """

    # Log the start time.
    logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Running select-features command")

    try:
            
        # Read in single h5 file
        logger.info(f"Reading h5 file from {args.input_file}")
        adata = sc.read_h5ad(args.input_file)
        adata.var_names_make_unique()
        describe_anndata(adata)

        # Run methods
        if "seurat" in args.methods:
            run_seurat(
                adata=adata,
                layer=args.layer,
                min_mean=args.min_mean,
                max_mean=args.max_mean,
                min_disp=args.min_disp,
                max_disp=args.max_disp,
                n_bins=args.n_bins,
                key_added="highly_variable_seurat",
            )

        if "seurat_v3" in args.methods:
            run_seurat_v3(
                adata=adata,
                layer=args.layer,
                n_top_genes=args.n_top_genes,
                span=args.span,
                key_added="highly_variable_seurat_v3",
            )

        if "cell_ranger" in args.methods:
            run_cell_ranger(
                adata=adata,
                min_mean=args.min_mean,
                max_mean=args.max_mean,
                min_disp=args.min_disp,
                max_disp=args.max_disp,
                n_bins=args.n_bins,
                key_added="highly_variable_cell_ranger",
            )

        if "deviance" in args.methods:
            run_deviance(
                adata=adata,
                layer=args.layer,
                n_top_genes=args.n_top_genes,
                key_added="highly_variable_deviance",
            )

        if "signac" in args.methods:
            run_signac(
                adata=adata,
                layer=args.layer,
                key_added="highly_variable_signac",
            )

        # Plot means vs dispersions for each method
        if "means" not in adata.var.columns or "dispersions" not in adata.var.columns:
            mtx = adata.layers[args.layer].A
            if is_counts(mtx) or is_mostly_counts(mtx, percent=0.95):
                logger.info(f"Data in layer {args.layer} is likely counts, normalizing for mean and dispersion calculation")
                sc.pp.normalize_total(adata, target_sum=1e4, layer=args.layer)
            sc.pp.highly_variable_genes(adata, layer=args.layer)
        highly_var_cols = [c for c in adata.var.columns if "highly_variable" in c]
        _, axes = plt.subplots(
            nrows=1,
            ncols=len(highly_var_cols),
            figsize=(len(highly_var_cols) * 5, 5),
        )
        for i, col in enumerate(highly_var_cols):
            ax = axes[i]
            sns.scatterplot(data=adata.var, x="means", y="dispersions", hue=highly_var_cols[i], s=5, ax=ax)
            ax.set_title(col)
        plt.savefig(os.path.join(args.output_dir, f"means_vs_dispersion_scatterplots.png"))
        plt.close()

        # Save a tsv with all the highly variable annotations
        adata.var.to_csv(os.path.join(args.output_dir, f"feature_metadata.tsv"), sep="\t")

        # Save the adata
        logger.info(
            f"Saving adata with highly variable gene annotations in `.var` to {os.path.join(args.output_dir, f'{args.output_prefix}.h5ad')}"
        )
        adata.write(os.path.join(args.output_dir, f"{args.output_prefix}.h5ad"))

        # Log the end time
        logger.info("Completed select-features")
        logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))

    # The exception allows user to end inference prematurely with CTRL-C.
    except KeyboardInterrupt:
        # If partial output has been saved, delete it.
        full_file = args.output_dir + args.output_prefix + ".h5ad"

        # Name of the filtered (cells only) file.
        if os.path.exists(full_file):
            os.remove(full_file)

        logger.info("Keyboard interrupt.  Terminated without saving\n")
