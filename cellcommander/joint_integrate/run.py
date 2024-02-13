"""Single run of joint_integrate, given input arguments."""

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
from cellcommander.utils import describe_mudata
from cellcommander.joint_integrate import consts
from cellcommander.joint_integrate.wnn import run_seurat_wnn


logger = logging.getLogger("cellcommander")


def run_joint_integrate(args: argparse.Namespace):
    """The full script for the command line tool to perform joint_integrate and filtering.

    Args:
        args: Inputs from the command line, already parsed using argparse.

    Note: Returns nothing, but writes output to a file(s) specified from command line.

    """

    # Log the start time.
    logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Running joint-integrate command")

    try:
            
        # Read in rna and atac data
        logger.info(f"Reading RNA AnnData from {args.rna_h5ad_path}")
        rna = sc.read_h5ad(args.rna_h5ad_path)
        logger.info(f"Reading ATAC AnnData from {args.atac_h5ad_path}")
        atac = sc.read_h5ad(args.atac_h5ad_path)

        # Find the common cells
        common_idx = list(set(rna.obs_names).intersection(set(atac.obs_names)))
        logger.info(f"Found {len(common_idx)} common cells between RNA and ATAC data, subsetting to these cells in a MuData object")
        rna = rna[common_idx].copy()
        atac = atac[common_idx].copy()
        mdata = mu.MuData({"rna": rna, "atac": atac})
        describe_mudata(mdata)

        # Run methods
        if "wnn" in args.methods:
            logger.info("Running Seurat WNN")
            run_seurat_wnn(
                mdata=mdata, 
                rna_mod_key="rna",
                atac_mod_key="atac",
                rna_obsm_key=args.rna_dim_reduction,
                atac_obsm_key=args.atac_dim_reduction,
                cluster_key=args.cluster_key, 
                cluster_resolution=args.clust_resolution,
                random_state=args.random_state
            )
            if args.cluster_key:
                with plt.rc_context():
                    mu.pl.embedding(
                        mdata, color=[args.cluster_key, "wnn_RNA_weight", "wnn_ATAC_weight"], ncols=2, basis="X_umap_wnn", frameon=False
                    )
                    plt.savefig(os.path.join(args.output_dir, f"wnn_umap.png"))
                    plt.close()
            del mdata.obsp["wnn_connectivities"]

        # Save the mdata
        logger.info(
            f"Saving adata with cell identity annotations in `.obs` to {os.path.join(args.output_dir, f'{args.output_prefix}.h5ad')}"
        )
        mdata.write(os.path.join(args.output_dir, f"{args.output_prefix}.h5mu"))

        # Log the end time
        logger.info("Completed joint-integrate")
        logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))

    # The exception allows user to end inference prematurely with CTRL-C.
    except KeyboardInterrupt:
        # If partial output has been saved, delete it.
        full_file = args.output_dir + args.output_prefix + ".h5mu"

        # Name of the filtered (cells only) file.
        if os.path.exists(full_file):
            os.remove(full_file)

        logger.info("Keyboard interrupt.  Terminated without saving\n")
