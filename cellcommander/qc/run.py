"""Single run of qc, given input arguments."""

import argparse
import logging
import os
import sys
import traceback
from datetime import datetime
from typing import Dict, Optional, Tuple, Union

import matplotlib
import pandas as pd
import psutil
import scanpy as sc

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)
from cellcommander.qc.rna import rna_qc_triplet_plot, rna_qc, rna_outlier_filter

logger = logging.getLogger("cellcommander")


def run_qc(args: argparse.Namespace):
    """The full script for the command line tool to perform qc and filtering.

    Args:
        args: Inputs from the command line, already parsed using argparse.

    Note: Returns nothing, but writes output to a file(s) specified from command line.

    """

    # Log the start time.
    logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Running qc command")

    try:
        # Read in single h5 file
        logger.info(f"Reading h5 file from {args.input_file}")
        adata = sc.read_10x_h5(args.input_file)
        adata.var_names_make_unique()

        # Read metadata file if provided
        if args.metadata_file is not None:
            logger.info(f"Reading metadata from {args.metadata_file}")
            if args.metadata_source is None:
                args.metadata_source = os.path.basename(args.metadata_file).split(".")[
                    0
                ]
            metadata = pd.read_csv(args.metadata_file, index_col=0)
            metadata.columns = [
                c + "_" + args.metadata_source for c in metadata.columns
            ]
            adata.obs = adata.obs.merge(metadata, left_index=True, right_index=True)

        # Run QC
        if args.mode == "rna":
            # Log mode
            logger.info("Running in rna mode")
            
            # Calculate and store qc metrics in adata.obs
            logger.info("Calculating QC metrics and storing in adata.obs")
            rna_qc(
                adata=adata,
                n_top_genes=args.n_top_features,
                total_counts_nmads=args.total_counts_nmads,
                n_genes_by_counts_nmads=args.n_features_by_counts_nmads,
                pct_counts_in_top_genes_nmads=args.pct_counts_in_top_features_nmads,
                pct_counts_mt_nmads=args.pct_counts_mt_nmads,
                pct_counts_mt_threshold=args.pct_counts_mt_threshold,
                n_genes_high_threshold=args.n_features_high_threshold,
                n_genes_low_threshold=args.n_features_low_threshold,
            )

            # Pre filtering plots of QC metrics
            logger.info("Generating pre-filtering QC plots")
            rna_qc_triplet_plot(
                adata, args.output_dir, "pre", total_counts_bins=100
            )

            if not args.no_filter:
                # Log strategy and based on outlier definition in adata.obs, filter cells
                logger.info(f"Filtering using strategy: `{args.filtering_strategy}`")
                adata, filtered_bc = rna_outlier_filter(
                    adata,
                    outlier_cols=["outlier", "mt_outlier"],
                )

                # Save the barcodes of filtered cells to output directory
                logger.info(
                    f"Saving filtered barcodes to {os.path.join(args.output_dir, 'filtered_barcodes.txt')}"
                )
                filtered_bc_path = os.path.join(args.output_dir, "filtered_barcodes.txt")
                filtered_bc.to_series().to_csv(
                    filtered_bc_path, sep="\t", index=False, header=False
                )

                # Plot QC metrics after filtering
                logger.info("Generating post-filtering QC plots")
                rna_qc_triplet_plot(
                    adata, args.output_dir, "post", total_counts_bins=100
                )
            
        elif args.mode == "atac":
            raise NotImplementedError("ATAC-seq QC not yet implemented yet")

        # Optionally filter genes if min_cells is provided as arg
        if args.min_cells is not None:
            logger.info(f"Filtering genes based on cell count. Number of genes before filtering: {adata.n_vars}")
            sc.pp.filter_genes(adata, min_cells=args.min_cells)
            logger.info(f"Number of genes after filtering: {adata.n_vars}")

         # Save the filtered adata
        logger.info(
            f"Saving filtered adata to {os.path.join(args.output_dir, f'{args.output_prefix}.h5ad')}"
        )
        adata.write(os.path.join(args.output_dir, f"{args.output_prefix}.h5ad"))

        # Log the end time
        logger.info("Completed qc")
        logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))

    # The exception allows user to end inference prematurely with CTRL-C.
    except KeyboardInterrupt:
        # If partial output has been saved, delete it.
        full_file = args.output_dir + args.output_prefix + ".h5ad"

        # Name of the filtered (cells only) file.
        if os.path.exists(full_file):
            os.remove(full_file)

        logger.info("Keyboard interrupt.  Terminated without saving\n")
