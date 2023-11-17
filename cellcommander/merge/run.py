"""Single run of merge, given input arguments."""

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
from cellcommander.merge import consts
from cellcommander.normalize.sctransform import run_sctransform
from cellcommander.reduce_dimensions.seurat import run_seurat_default

logger = logging.getLogger("cellcommander")


def run_merge(args: argparse.Namespace):
    """The full script for the command line tool to perform merge and filtering.

    Args:
        args: Inputs from the command line, already parsed using argparse.

    Note: Returns nothing, but writes output to a file(s) specified from command line.

    """

    # Log the start time.
    logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Running merge command")

    try:
            
        # Read in the AnnData(s)
        if args.names is None:
            args.names = [str(i) for i in range(len(args.input_files))]
            logger.info("No sample names provided. Using indices as names.")
        if args.make_bcs_unique:
            logger.info("Making barcodes unique by prepending sample names.")
        if len(args.input_files) > 1:
            logger.info("Reading in multiple AnnData objects.")
            adata_list = []
            samples = []
            for i, input_file in enumerate(args.input_files):
                sample = args.names[i]
                logger.info(f"Reading in AnnData for sample {sample}")
                adata = sc.read_h5ad(input_file)
                describe_anndata(adata)
                if args.make_bcs_unique:
                    logger.info(f"Adding sample name {sample} to barcode")
                    adata.obs.index = sample + "#" + adata.obs.index
                if args.layer is None:
                    logger.info("No layer provided. Assuming counts are in adata.X.")
                    adata_cp = sc.AnnData(adata.X.copy(), obs=adata.obs.copy(), var=adata.var.copy())
                else:
                    logger.info(f"Using layer {args.layer} as counts.")
                    adata_cp = sc.AnnData(adata.layers[args.layer].copy(), obs=adata.obs.copy(), var=adata.var.copy())
                adata_list.append(adata_cp)
                samples.append(sample)

            # Concat samples
            logger.info("Concatenating AnnData objects for samples")
            adata_concat = ad.concat(adata_list, label=args.names_key, keys=samples)
            del adata_list, adata, adata_cp

        else:
            logger.info("Reading in single AnnData object.")
            input_file = args.input_files[0]
            adata_concat = sc.read_h5ad(input_file)
        
        # Save counts in adata_concat.raw
        adata_concat.raw = adata_concat.copy()

        # Save counts in adata_concat.layers["counts"]
        if args.layer is not None:
            adata_concat.layers[args.layer] = adata_concat.X.copy()

        # Save the adata
        logger.info(
            f"Saving merged adata to {os.path.join(args.output_dir, f'{args.output_prefix}.h5ad')}"
        )
        adata_concat.write(os.path.join(args.output_dir, f"{args.output_prefix}.h5ad"))

        # Log the end time
        logger.info("Completed merge")
        logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))

    # The exception allows user to end inference prematurely with CTRL-C.
    except KeyboardInterrupt:
        # If partial output has been saved, delete it.
        full_file = args.output_dir + args.output_prefix + ".h5ad"

        # Name of the filtered (cells only) file.
        if os.path.exists(full_file):
            os.remove(full_file)

        logger.info("Keyboard interrupt. Terminated without saving\n")
