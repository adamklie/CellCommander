"""Single run of remove_background, given input arguments."""

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
from cellcommander.remove_background import consts
from cellcommander.remove_background.soupx import run_soupx

logger = logging.getLogger("cellcommander")


def run_remove_background(args: argparse.Namespace):
    """The full script for the command line tool to perform remove_background and filtering.

    Args:
        args: Inputs from the command line, already parsed using argparse.

    Note: Returns nothing, but writes output to a file(s) specified from command line.

    """

    # Log the start time.
    logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Running remove-background command")

    try:
        # Read in single h5 file
        logger.info(f"Reading h5 file from {args.input_file}")
        adata = sc.read_h5ad(args.input_file)
        adata.var_names_make_unique()
        describe_anndata(adata)

        # Run menthods
        if "soupx" == args.method:
            # Read in raw data and markers
            logger.info(f"Reading raw h5 file from {args.input_file}")
            adata_raw = sc.read_10x_h5(args.raw_h5_path)
            adata_raw.var_names_make_unique()
            describe_anndata(adata_raw)
            soupx_markers = pd.read_csv(args.markers_path, sep="\t")

            if args.sample_name is not None:
                logger.info(f"Adding sample name {args.sample_name} to barcode for raw data")
                logger.info(f"Before: {adata_raw.obs.index[0]}")
                adata_raw.obs.index = args.sample_name + "#" + adata_raw.obs.index
                logger.info(f"After: {adata_raw.obs.index[0]}")

            # Run SoupX
            logger.info("Running backgorund removal with SoupX.")
            run_soupx(
                adata,
                adata_raw,
                soupx_markers,
                layer=args.layer,
                clust_num_hvgs=args.clust_num_hvgs,
                clust_n_neighbors=args.clust_n_neighbors,
                clust_n_components=args.clust_n_components,
                clust_resolution=args.clust_resolution,
                umap_min_distance=args.umap_min_distance,
                random_state=args.random_state,
                outdir_path=args.output_dir,
            )
            # Recalculate QC metrics
            logger.info("Recalculating QC metrics")
            barcode_qc, gene_qc = sc.pp.calculate_qc_metrics(
                adata,
                qc_vars=["mt", "ribo"],
                inplace=False,
                percent_top=[20],
                log1p=False,
            )
            barcode_qc.columns = [f"{col}_post_soupx" for col in barcode_qc.columns]
            gene_qc.columns = [f"{col}_post_soupx" for col in gene_qc.columns]
            adata.obs = adata.obs.merge(barcode_qc, left_index=True, right_index=True)
            adata.obs["log1p_total_counts_post_soupx"] = np.log1p(adata.obs["total_counts_post_soupx"])
            adata.var = adata.var.merge(gene_qc, left_index=True, right_index=True)
        else:
            raise ValueError(f"Method {args.method} not recognized.")

        # Save the adata
        logger.info(
            f"Saving adata with corrected counts to {os.path.join(args.output_dir, f'{args.output_prefix}.h5ad')}"
        )
        adata.write(os.path.join(args.output_dir, f"{args.output_prefix}.h5ad"))

        # Log the end time
        logger.info("Completed remove-background")
        logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))

    # The exception allows user to end inference prematurely with CTRL-C.
    except KeyboardInterrupt:
        # If partial output has been saved, delete it.
        full_file = args.output_dir + args.output_prefix + ".h5ad"

        # Name of the filtered (cells only) file.
        if os.path.exists(full_file):
            os.remove(full_file)

        logger.info("Keyboard interrupt.  Terminated without saving\n")
