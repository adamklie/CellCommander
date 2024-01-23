"""Single run of qc, given input arguments."""

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
from muon import atac as ac
from mudata import MuData

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)
from cellcommander.utils import describe_anndata, describe_mudata
from cellcommander.qc import consts
from cellcommander.qc.atac import atac_outlier_filter, atac_qc, atac_qc_triplet_plot
from cellcommander.qc.rna import rna_outlier_filter, rna_qc, rna_qc_triplet_plot
from cellcommander.qc.utils import sample_level_metrics

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
        
        if args.multimodal_input:
            # Read in multimodal h5 file
            logger.info("Reading multimodal input using Muon. Expecting to return a MuData object.")
            data = mu.read_10x_h5(args.input_file)
            data.var_names_make_unique()
            describe_mudata(data)

        else:
            if args.mode == "rna":
                # Read in single modality h5 file
                logger.info("Reading single modality RNA input using Scanpy. Expecting to return an AnnData object.")
                data = sc.read_10x_h5(args.input_file)
                data.var_names_make_unique()
                describe_anndata(data)

            elif args.mode == "atac":
                # Read in single modality h5 file
                logger.info("Reading single modality ATAC input using Muon. Expecting to return an AnnData object.")
                data = ac.read_10x_h5(args.input_file)
                data.var_names_make_unique()
                describe_anndata(data)
        
        # If rna
        if args.mode == "rna":

            # Grab the RNA only
            if isinstance(data, MuData) and "rna" in data.mod:
                adata = data.mod["rna"]
                
                # Keep only features that are genes, just in case
                adata[:, adata.var["feature_types"] == "Gene Expression"]
            else:
                adata = data

            # Sample level metrics
            metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]

        
        # If atac
        elif args.mode == "atac":
        
            # Grab the ATAC only
            if isinstance(data, MuData) and "atac" in data.mod:
                adata = data.mod["atac"]
                
                # Keep only features that are Peaks, just in case
                adata[:, adata.var["feature_types"] == "Peaks"]
            else:
                adata = data 

            # Make sure we have a fragments file. TODO: Move to cli checks
            if "files" not in adata.uns:
                logger.info("No fragment file found in adata.uns['files']. Checking arguments.")
                if args.fragments_file is not None:
                    logger.info(f"Using provided fragments file from {args.fragments_file}")
                    adata.uns["files"] = {}
                    adata.uns["files"]["fragments"] = args.fragments_file
                else:
                    logger.info("No fragments file provided. Checking for fragments.tsv.gz or atac_fragments.tsv.gz file in input directory.")
                    if os.path.exists(os.path.join(os.path.dirname(args.input_file), "fragments.tsv.gz")):
                        logger.info("Found fragments.tsv.gz file in input directory.")
                        adata.uns["files"] = {}
                        adata.uns["files"]["fragments"] = os.path.join(os.path.dirname(args.input_file), "fragments.tsv.gz")
                    elif os.path.exists(os.path.join(os.path.dirname(args.input_file), "atac_fragments.tsv.gz")):
                        logger.info("Found atac_fragments.tsv.gz file in input directory.")
                        adata.uns["files"] = {}
                        adata.uns["files"]["fragments"] = os.path.join(os.path.dirname(args.input_file), "atac_fragments.tsv.gz")
                    else:
                        logger.info("No fragments file found. Cannot calculate metrics.")
                        sys.exit(1)
                
            # Sample level metrics
            metrics = [
                "frip",
                "pct_counts_mt",
                "nucleosome_signal",
                "tss_score",
                "total_counts",
            ]

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
            adata.obs = adata.obs.merge(metadata, left_index=True, right_index=True, how="left")
    
        # Add sample name to barcode with # in between
        if args.sample_name is not None:
            logger.info(f"Adding sample name {args.sample_name} to barcode")
            logger.info(f"Before: {adata.obs.index[0]}")
            adata.obs.index = args.sample_name + "#" + adata.obs.index
            logger.info(f"After: {adata.obs.index[0]}")
            
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
                n_genes_by_counts_nmads=args.n_features_nmads,
                pct_counts_in_top_genes_nmads=args.pct_counts_in_top_features_nmads,
                pct_counts_mt_nmads=args.pct_counts_mt_nmads,
                pct_counts_mt_hi=args.pct_counts_mt_hi,
                n_genes_hi=args.n_features_hi,
                n_genes_low=args.n_features_low,
            )

            # Pre filtering plots of QC metrics
            logger.info("Generating pre-filtering QC plots")
            rna_qc_triplet_plot(adata, args.output_dir, "pre", total_counts_bins=100)

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
                filtered_bc_path = os.path.join(
                    args.output_dir, "filtered_barcodes.txt"
                )
                filtered_bc.to_series().to_csv(
                    filtered_bc_path, sep="\t", index=False, header=False
                )

                # Filter out any user defined lists of barcodes
                if args.barcode_exclusion_list_paths is not None:
                    for path in args.barcode_exclusion_list_paths:
                        logger.info(f"Excluding barcodes from {path}")
                        barcodes_to_filter = pd.read_csv(path, header=None)[0].tolist()
                        adata = adata[~adata.obs.index.isin(barcodes_to_filter), :]
                        logger.info(f"Number of cells after filtering: {adata.n_obs}")

                # Plot QC metrics after filtering
                logger.info("Generating post-filtering QC plots")
                rna_qc_triplet_plot(
                    adata, args.output_dir, "post", total_counts_bins=100
                )

        elif args.mode == "atac":
            # Log mode
            logger.info("Running in atac mode")

            # If metadata source is cellranger, we can add some additional metrics
            if args.metadata_source == "cellranger":
                if args.multimodal_input:
                    adata.obs["frip"] = adata.obs["atac_peak_region_fragments_cellranger"]/adata.obs["atac_fragments_cellranger"]
                    adata.obs["pct_counts_mt"] = (adata.obs["atac_mitochondrial_reads_cellranger"]/adata.obs["atac_raw_reads_cellranger"])

            # Calculate and store qc metrics in adata.obs
            logger.info("Calculating QC metrics and storing in adata.obs")
            if args.atac_qc_tool == "muon":
                atac_qc(
                    adata=adata,
                    n_top_features=args.n_top_features,
                    pct_counts_in_top_features_nmads=args.pct_counts_in_top_features_nmads,
                    n_for_ns_calc=args.n_for_ns_calc,
                    ns_nmads=args.ns_nmads,
                    ns_hi=args.ns_hi,
                    n_tss=args.n_tss,
                    tss_annot=None,
                    tss_nmads=args.tss_nmads,
                    tss_hi=args.tss_hi,
                    tss_low=args.tss_low,
                    total_counts_nmads=args.total_counts_nmads,
                    total_counts_low=args.total_counts_low,
                    total_counts_hi=args.total_counts_hi,
                    n_features_nmads=args.n_features_nmads,
                    n_features_low=args.n_features_low,
                    random_state=args.random_state,
                )

            # Pre filtering plots of QC metrics
            logger.info("Generating pre-filtering QC plots")
            atac_qc_triplet_plot(
                adata=adata,
                outdir_path=args.output_dir,
                output_prefix="pre",
                total_counts_bins=100,
                total_counts_low=args.total_counts_low,
                total_counts_hi=args.total_counts_hi,
                tss_low=args.tss_low,
                tss_hi=args.tss_hi,
            )

            if not args.no_filter:
                # Log strategy and based on outlier definition in adata.obs, filter cells
                logger.info(f"Filtering using strategy: `{args.filtering_strategy}`")
                adata, filtered_bc = atac_outlier_filter(
                    adata,
                    outlier_cols=["outlier", "mt_outlier"],
                )

                # Save the barcodes of filtered cells to output directory
                logger.info(
                    f"Saving filtered barcodes to {os.path.join(args.output_dir, 'filtered_barcodes.txt')}"
                )
                filtered_bc_path = os.path.join(
                    args.output_dir, "filtered_barcodes.txt"
                )
                filtered_bc.to_series().to_csv(
                    filtered_bc_path, sep="\t", index=False, header=False
                )

                # Plot QC metrics after filtering
                logger.info("Generating post-filtering QC plots")
                atac_qc_triplet_plot(
                    adata, args.output_dir, "post", total_counts_bins=100
                )

        # Optionally filter genes if min_cells is provided as arg
        if args.min_cells_per_feature is not None:
            logger.info(
                f"Filtering genes based on cell count. Number of genes before filtering: {adata.n_vars}"
            )
            sc.pp.filter_genes(adata, min_cells=args.min_cells_per_feature)
            logger.info(f"Number of genes after filtering: {adata.n_vars}")

        # Add counts to layers
        adata.layers["counts"] = adata.X.copy()
        
        # Save the filtered adata
        logger.info(
            f"Saving filtered adata to {os.path.join(args.output_dir, f'{args.output_prefix}.h5ad')}"
        )
        adata.write(os.path.join(args.output_dir, f"{args.output_prefix}.h5ad"))

        # Sample level metrics dict
        sample_metrics = sample_level_metrics(adata, metrics)
        sample_metrics_df = pd.DataFrame.from_dict(sample_metrics, orient="index")
        sample_metrics_df.columns = ["median"]
        sample_metrics_df.to_csv(
            os.path.join(args.output_dir, "sample_metrics.tsv"), sep="\t"
        )

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
