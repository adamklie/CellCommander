"""Single run of annotate, given input arguments."""

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
from cellcommander.annotate import consts
from cellcommander.annotate.utils import check_marker_genes, get_user_input


logger = logging.getLogger("cellcommander")


def run_annotate(args: argparse.Namespace):
    """The full script for the command line tool to perform annotate and filtering.

    Args:
        args: Inputs from the command line, already parsed using argparse.

    Note: Returns nothing, but writes output to a file(s) specified from command line.

    """

    # Log the start time.
    logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Running annotate command")

    try:
            
        # Read in single h5 file
        logger.info(f"Reading h5 file from {args.input_file}")
        adata = sc.read_h5ad(args.input_file)
        adata.var_names_make_unique()
        describe_anndata(adata)
        
        # Make a copy of the adata so we don't overwrite the original
        adata_plot = adata.copy()

        # Make .X the passed in layer if provided
        if args.layer is not None:
            adata_plot.X = adata_plot.layers[args.layer]
            logger.info(f"Using {args.layer} layer as .X")
        
        # Get the dim reductions
        adata_plot.obsm["X_umap"] = adata_plot.obsm[args.dim_reduction]
        logger.info(f"Using {args.dim_reduction} as .obsm['X_umap']")

        # Get the clustering
        if args.cluster_key is None:
            logger.info("No cluster key provided, running leiden clustering")
            adata_plot.obsm["X_pca"] = adata_plot.obsm[args.clust_dim_reduction]
            logger.info(f"Using {args.clust_dim_reduction} as .obsm['X_pca']")
            sc.pp.neighbors(adata_plot, use_rep="X_pca", n_neighbors=args.clust_n_neighbors, random_state=args.random_state, n_pcs=args.clust_n_components)
            cluster_key = f"annotate_leiden_{args.clust_resolution}"
            sc.tl.leiden(adata_plot, resolution=args.clust_resolution, random_state=args.random_state, key_added=cluster_key)

        else:
            logger.info(f"Using {args.cluster_key} as cluster key")
            cluster_key = args.cluster_key

        # Plot the clustering as it's own UMAP
        with plt.rc_context({"figure.figsize": (5, 5)}):
            sc.pl.umap(adata_plot, color=cluster_key, legend_loc="on data", show=False)
            plt.savefig(os.path.join(args.output_dir, f"{cluster_key}_umap.png"))
            plt.close()
        
        # Plot markers if file provided
        if args.marker_gene_list is not None:
            logger.info(f"Plotting markers from {args.marker_gene_list}")
            marker_genes_df = pd.read_csv(args.marker_gene_list)
            marker_genes_dict = marker_genes_df.groupby("cell_id")["gene"].apply(list).to_dict()
            marker_genes_in_data = check_marker_genes(adata_plot, marker_genes_dict)
            for ct in marker_genes_in_data.keys():
                logger.info(f"Plotting markers for {ct} to {os.path.join(args.output_dir, 'marker_plots', f'{ct}_markers.png')}")
                if not os.path.exists(os.path.join(args.output_dir, "marker_plots")):
                        os.makedirs(os.path.join(args.output_dir, "marker_plots"))
                logger.info(f"Relavent {ct} markers: {marker_genes_in_data[ct]}")
                with plt.rc_context():
                    sc.pl.umap(
                        adata_plot,
                        color=marker_genes_in_data[ct],
                        vmin=0,
                        vmax="p99.5",  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
                        sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
                        frameon=False,
                        cmap="Reds",  # or choose another color map e.g. from here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
                    )
                    plt.savefig(os.path.join(args.output_dir, "marker_plots", f"{ct}_markers_umap.png"))
                    plt.close()

            # Plot a dotplot of the markers
            with plt.rc_context():
                logger.info(f"Plotting dotplot of markers to {os.path.join(args.output_dir, 'marker_gene_dotplot.png')}")
                sc.pl.dotplot(
                    adata_plot,
                    groupby=cluster_key,
                    var_names=marker_genes_in_data,
                    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
                )
                plt.savefig(os.path.join(args.output_dir, "marker_gene_dotplot.png"))
                plt.close()

            # Get quantified marker genes and plot a dotplot
            sc.tl.rank_genes_groups(
                adata_plot, groupby=cluster_key, method="wilcoxon", key_added=f"dea_{cluster_key}"
            )
            # Filter these results to genes that are highly specific to each cluster
            sc.tl.filter_rank_genes_groups(
                adata_plot,
                min_in_group_fraction=0.2,
                max_out_group_fraction=0.2,
                key=f"dea_{cluster_key}",
                key_added=f"dea_{cluster_key}_filtered",
            )
            with plt.rc_context():
                # Dot plot those markers
                logger.info(f"Plotting dotplot of DEA markers to {os.path.join(args.output_dir, f'dea_{cluster_key}_filtered_dotplot.png')}")
                sc.pl.rank_genes_groups_dotplot(
                    adata_plot,
                    groupby=cluster_key,
                    standard_scale="var",
                    n_genes=5,
                    key=f"dea_{cluster_key}_filtered",
                )
                plt.savefig(os.path.join(args.output_dir, f"dea_{cluster_key}_filtered_dotplot.png"))
                plt.close()

        if "manual" in args.methods:
            logger.info("Running manual annotation, will prompt user to annotate clusters")
            numerically_sorted_clusters = sorted(adata_plot.obs[cluster_key].unique(), key=lambda x: int(x))
            if args.marker_gene_list is not None:
                cell_types_to_choose_from = marker_genes_df.iloc[:, 1].unique().tolist() + ["other"]
            else:
                cell_types_to_choose_from = None
            cluster_to_celltype = get_user_input(numerically_sorted_clusters, cell_types_to_choose_from)

            # Add the annotation to obs
            adata_plot.obs[args.annotation_key] = adata_plot.obs[cluster_key].map(cluster_to_celltype)

            with plt.rc_context({"figure.figsize": (5, 5)}):
                sc.pl.umap(adata_plot, color=[args.annotation_key, cluster_key], legend_loc="on data", show=False)
                plt.savefig(os.path.join(args.output_dir, f"{args.annotation_key}umap.png"))
                plt.close()

        # Save a tsv
        adata_plot.obs.to_csv(os.path.join(args.output_dir, f"{args.output_prefix}_metadata.tsv"), sep="\t")

        # Save the adata
        #del adata.uns[f"dea_{cluster_key}_filtered"]
        adata.obs[cluster_key] = adata_plot.obs[cluster_key]
        adata.obs[args.annotation_key] = adata_plot.obs[args.annotation_key]
        logger.info(
            f"Saving adata with cell identity annotations in `.obs` to {os.path.join(args.output_dir, f'{args.output_prefix}.h5ad')}"
        )
        adata.write(os.path.join(args.output_dir, f"{args.output_prefix}.h5ad"))

        # Log the end time
        logger.info("Completed annotate")
        logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))

    # The exception allows user to end inference prematurely with CTRL-C.
    except KeyboardInterrupt:
        # If partial output has been saved, delete it.
        full_file = args.output_dir + args.output_prefix + ".h5ad"

        # Name of the filtered (cells only) file.
        if os.path.exists(full_file):
            os.remove(full_file)

        logger.info("Keyboard interrupt.  Terminated without saving\n")
