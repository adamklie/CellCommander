"""Single run of detect_doublets, given input arguments."""

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
from cellcommander.detect_doublets import consts
from cellcommander.detect_doublets.scrublet import scrublet_recipe
from cellcommander.detect_doublets.scDblFinder import scDblFinder_recipe
from cellcommander.detect_doublets.amulet import amulet_recipe

logger = logging.getLogger("cellcommander")


def run_detect_doublets(args: argparse.Namespace):
    """The full script for the command line tool to perform detect_doublets and filtering.

    Args:
        args: Inputs from the command line, already parsed using argparse.

    Note: Returns nothing, but writes output to a file(s) specified from command line.

    """

    # Log the start time.
    logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Running detect_doublets command")

    try:
        # Read in single h5 file
        logger.info(f"Reading h5 file from {args.input_file}")
        if args.input_file.endswith(".h5ad"):
            adata = sc.read_h5ad(args.input_file)
        elif args.input_file.endswith(".h5"):
            adata = sc.read_10x_h5(args.input_file)
        else:
            adata = sc.read(args.input_file)

        # Make sure var_names are unique
        adata.var_names_make_unique()
        describe_anndata(adata)

        # Add fragments file to adata.uns
        if args.fragments_file is not None:
            logger.info(f"Adding fragments file to adata.uns['files']['fragments']")
            adata.uns["files"] = {"fragments": args.fragments_file}

        # Get a prelim clustering of the data for plotting
        logger.info("Generating preliminary clustering for plotting UMAP with scores on it.")
        adata_pp = adata.copy()
        sc.pp.filter_cells(adata_pp, min_genes=200)
        sc.pp.filter_genes(adata_pp, min_cells=20)
        sc.pp.normalize_total(adata_pp)
        sc.pp.log1p(adata_pp)
        sc.pp.highly_variable_genes(adata_pp, n_top_genes=args.clust_num_hvgs)
        adata_pp = adata_pp[:, adata_pp.var.highly_variable]
        sc.pp.pca(adata_pp, n_comps=args.clust_n_components, random_state=args.random_state)
        sc.pp.neighbors(adata_pp, n_neighbors=args.clust_n_neighbors, n_pcs=args.clust_n_components, random_state=args.random_state)
        sc.tl.umap(adata_pp, min_dist=args.umap_min_distance, random_state=args.random_state)
        sc.tl.leiden(adata_pp, resolution=args.clust_resolution, random_state=args.random_state, key_added=f"pre_doublet_filter_leiden_{args.clust_resolution}")

        # Get the methods to run
        if args.method != "consensus":
            methods = [args.method]
        else:
            methods = args.consensus_methods
            logger.info(f"Running {args.consensus_methods} for doublet detection and taking the {args.consensus_strategy} of the results.")

        # Run methods
        if "scrublet" in methods:
            scrublet_recipe(
                adata=adata, 
                outdir_path=args.output_dir, 
                random_state=args.random_state
            )

        if "scDblFinder" in methods:
            scDblFinder_recipe(
                adata=adata, 
                outdir_path=args.output_dir, 
                atac_params=args.scDblFinder_atac_params,
                random_state=args.random_state
            )
            
        if "amulet" in methods:
            amulet_recipe(
                adata=adata, 
                outdir_path=args.output_dir, 
            )
        
        if "cellranger" in methods:
            if "excluded_reason_cellranger" in adata.obs.columns:
                logger.info("Adding cellranger predicted doublets from adata.obs['excluded_reason_cellranger']")
                adata.obs["cellranger_predicted_doublet"] = (adata.obs["excluded_reason_cellranger"] == "1")
            else:
                logger.info("No cellranger predicted doublets column found, setting all to False")
                adata.obs["cellranger_predicted_doublet"] = False

        # Get doublets based on the consensus strategy
        if args.method == "consensus":
            doublet_cols = [c for c in adata.obs.columns if "_predicted_doublet" in c]
            if args.consensus_strategy == "union":
                adata.obs["doublet_filter"] = adata.obs[doublet_cols].any(axis=1)
            elif args.consensus_strategy == "intersection":
                adata.obs["doublet_filter"] = adata.obs[doublet_cols].all(axis=1)
            elif args.consensus_strategy == "majority":
                adata.obs["doublet_filter"] = adata.obs[doublet_cols].sum(axis=1) > (len(doublet_cols) / 2)
        else:
            adata.obs["doublet_filter"] = adata.obs[f"{methods[0]}_predicted_doublet"]
        
        # Filter out any user defined lists of barcodes
        barcodes_to_filter = []
        if args.barcode_exclusion_list_paths is not None:
            for path in args.barcode_exclusion_list_paths:
                logger.info(f"Excluding barcodes from {path}")
                barcodes_to_filter.append(pd.read_csv(path, header=None)[0].tolist())
            adata.obs["additional_barcode_exclusion_list"] = adata.obs.index.isin(barcodes_to_filter)
                        
        # Plot doublet scores on UMAP
        doublet_score_columns = [c for c in adata.obs.columns if "_score" in c]
        adata_pp.obs[doublet_score_columns + ["doublet_filter"]] = adata.obs[doublet_score_columns + ["doublet_filter"]]
        adata_pp.obs["doublet_filter"] = adata.obs["doublet_filter"].astype("category")
        adata.obs[f"pre_doublet_filter_leiden_{args.clust_resolution}"] = adata_pp.obs[f"pre_doublet_filter_leiden_{args.clust_resolution}"]
        plot_cols = [f"pre_doublet_filter_leiden_{args.clust_resolution}"] + doublet_score_columns + ["doublet_filter"]
        if args.barcode_exclusion_list_paths is not None:
            adata_pp.obs["additional_barcode_exclusion_list"] = adata.obs["additional_barcode_exclusion_list"].astype("category")
            plot_cols.append("additional_barcode_exclusion_list")
        logger.info("Plotting doublet scores on UMAP")
        with plt.rc_context():
            sc.pl.umap(
                adata_pp,
                color=plot_cols,
                vmin=0,
                vmax="p99", 
                sort_order=False, 
                frameon=False,
                cmap="Reds", 
            )
            plt.savefig(os.path.join(args.output_dir, "doublet_scores_umap.png"))
            plt.close()
        
        # Filter doublets if specified
        if not args.no_filter:
            logger.info(f"Filtering barcodes determined to be doublets")
            doublet_bcs = adata.obs[adata.obs["doublet_filter"] == True].index.tolist()
            adata = adata[~adata.obs["doublet_filter"], :].copy()
            logger.info(f"Filtered out {len(doublet_bcs)} doublets")
            logger.info(f"Saving doublet barcodes to {os.path.join(args.output_dir, 'doublet_barcodes.txt')}")
            with open(os.path.join(args.output_dir, "doublet_barcodes.txt"), "w") as f:
                for bc in doublet_bcs:
                    f.write(bc + "\n")

        # Save the adata
        logger.info(
            f"Saving adata to {os.path.join(args.output_dir, f'{args.output_prefix}.h5ad')}"
        )
        adata.write(os.path.join(args.output_dir, f"{args.output_prefix}.h5ad"))

        # Log the end time
        logger.info("Completed detect_doublets")
        logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))

    # The exception allows user to end inference prematurely with CTRL-C.
    except KeyboardInterrupt:
        # If partial output has been saved, delete it.
        full_file = args.output_dir + args.output_prefix + ".h5ad"

        # Name of the filtered (cells only) file.
        if os.path.exists(full_file):
            os.remove(full_file)

        logger.info("Keyboard interrupt.  Terminated without saving\n")
