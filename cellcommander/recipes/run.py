"""Single run of recipes, given input arguments."""

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
from cellcommander.recipes import consts
from cellcommander.recipes.snapatac2 import single_sample_recipe

logger = logging.getLogger("cellcommander")


def run_recipes(args: argparse.Namespace):
    """The full script for the command line tool to perform recipes and filtering.

    Args:
        args: Inputs from the command line, already parsed using argparse.

    Note: Returns nothing, but writes output to a file(s) specified from command line.

    """

    # Log the start time.
    logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Running recipes command")

    try:
            
        # Read in single h5 file
        if len(args.input_files) > 1:
            raise NotImplementedError("Only single h5 file input supported for now")
        else:
            input_file = args.input_files[0]

        # Read in params if passed in
        if args.params_file is not None:
            logger.info(f"Reading params file from {args.params_file}")
            with open(args.params_file) as f:
                params = yaml.load(f, Loader=yaml.FullLoader)
        else:
            logger.info(f"No params file passed in, using default params")
            params = consts.DEFAULT_PARAMS[args.method][args.mode]
        
        # Dump params to output dir in a yaml file
        logger.info(f"Dumping params to {os.path.join(args.output_dir, 'params.yaml')}")
        with open(os.path.join(args.output_dir, "params.yaml"), "w") as f:
            yaml.dump(params, f)

        # Check if user passed in RNA metadata
        if args.metadata_path is not None:
            logger.info(f"Reading metadata from {args.metadata_path}, should have bcs that match in first column")
            if args.metadata_path.endswith(".csv"):
                metadata = pd.read_csv(args.metadata_path, index_col=0)
            elif args.metadata_path.endswith(".tsv"):
                metadata = pd.read_csv(args.metadata_path, sep="\t", index_col=0)
        else:
            metadata = None

        if args.method == "snapatac2":
            if args.mode == "single-sample":
                single_sample_recipe(
                    frag_file=input_file,
                    outdir_path=args.output_dir,
                    sample_name=args.sample_name,
                    bin_size=params["feature_selection"]["bin_size"],
                    num_features=params["feature_selection"]["num_features"],
                    min_load_num_fragments=params["io"]["min_load_num_fragments"],
                    min_tsse=params["qc"]["min_tsse"],
                    min_num_fragments=params["qc"]["min_num_fragments"],
                    max_num_fragments=params["qc"]["max_num_fragments"],
                    sorted_by_barcode=params["io"]["sorted_by_barcode"],
                    chunk_size=params["io"]["chunk_size"],
                    clustering_resolution=params["analysis"]["clustering_resolution"],
                    gene_activity=params["io"]["gene_activity"],
                    metadata=metadata,
                    additional_doublets=params["qc"]["additional_doublets"],
                    save_intermediate=params["io"]["save_intermediate"],
                )
            else:
                raise NotImplementedError(f"Mode {args.mode} not implemented for tool {args.method}")
        else:
            raise NotImplementedError(f"Tool {args.method} not implemented")

        # Log the end time
        logger.info("Completed recipes")
        logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))

    # The exception allows user to end inference prematurely with CTRL-C.
    except KeyboardInterrupt:
        # If partial output has been saved, delete it.
        full_file = args.output_dir + args.output_prefix + ".h5ad"

        # Name of the filtered (cells only) file.
        if os.path.exists(full_file):
            os.remove(full_file)

        logger.info("Keyboard interrupt. Terminated without saving\n")
