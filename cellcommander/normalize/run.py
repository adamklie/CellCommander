"""Single run of normalize, given input arguments."""

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

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)
from cellcommander.utils import describe_anndata
from cellcommander.normalize import consts
from cellcommander.normalize.tfidf import tfidf_recipe
from cellcommander.normalize.sctransform import sctransform_recipe
from cellcommander.normalize.log1p import log1p_recipe


logger = logging.getLogger("cellcommander")


def run_normalize(args: argparse.Namespace):
    """The full script for the command line tool to perform normalize and filtering.

    Args:
        args: Inputs from the command line, already parsed using argparse.

    Note: Returns nothing, but writes output to a file(s) specified from command line.

    """

    # Log the start time.
    logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Running normalize command")

    try:
            
        # Read in single h5 file
        logger.info(f"Reading h5 file from {args.input_file}")
        adata = sc.read_h5ad(args.input_file)
        adata.var_names_make_unique()
        describe_anndata(adata)

        # Run methods
        if "sctransform" in args.methods:
            logger.info(f"Using SCTransform to normalize data.")
            sctransform_recipe(
                adata=adata, 
                outdir_path=args.output_dir, 
                filter_genes=args.filter_features,
                save_normalized_mtx=args.save_normalized_mtx)

        if "tfidf" in args.methods:
            logger.info(f"Using TF-IDF to normalize data.")
            tfidf_recipe(adata, args.output_dir, scale_factor=args.tfidf_scale_factor, save_normalized_mtx=args.save_normalized_mtx)

        if "log1p" in args.methods:
            logger.info(f"Using log1p normalization.")
            log1p_recipe(
                adata=adata, 
                outdir_path=args.output_dir, 
                save_normalized_mtx=args.save_normalized_mtx)
            
        # Save the adata
        logger.info(
            f"Saving adata with normalizations in obsm and layers to {os.path.join(args.output_dir, f'{args.output_prefix}.h5ad')}"
        )
        adata.write(os.path.join(args.output_dir, f"{args.output_prefix}.h5ad"))

        # Log the end time
        logger.info("Completed normalize")
        logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))

    # The exception allows user to end inference prematurely with CTRL-C.
    except KeyboardInterrupt:
        # If partial output has been saved, delete it.
        full_file = args.output_dir + args.output_prefix + ".h5ad"

        # Name of the filtered (cells only) file.
        if os.path.exists(full_file):
            os.remove(full_file)

        logger.info("Keyboard interrupt.  Terminated without saving\n")
