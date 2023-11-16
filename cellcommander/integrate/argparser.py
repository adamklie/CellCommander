import argparse

from cellcommander.integrate import consts


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for integrate.

    Args:
        subparsers: Parser object before addition of arguments specific to integrate.

    Returns:
        parser: Parser object with additional parameters.

    """

    subparser = subparsers.add_parser(
        "integrate",
        description="Integrates multiple objects together.",
        help="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparser.add_argument(
        "--input_h5ad_path",
        nargs=None,
        type=str,
        dest="input_file",
        required=True,
        help="Data file on which to run tool. "
        "The following input formats are supported: "
        "AnnData (.h5ad). "
    )
    subparser.add_argument(
        "--outdir_path",
        nargs=None,
        type=str,
        dest="output_dir",
        required=True,
        help="Output directory location. " 
        "If it does not exist, it will be created.",
    )
    subparser.add_argument(
        "--output_prefix",
        nargs=None,
        type=str,
        dest="output_prefix",
        required=False,
        default="integrate",
        help="Prefix for output files. " "If not provided, the prefix will be 'reduce-dimensions'.",
    )
    subparser.add_argument(
        "--batch-key",
        nargs=None,
        type=str,
        dest="batch_key",
        required=False,
        default="sample",
        help="Name of the key in the AnnData.obs to use for integration. ",
    )
    subparser.add_argument(
        "--counts-key",
        nargs=None,
        type=str,
        dest="counts_key",
        required=False,
        default=None,
        help="Key in the AnnData.layers that contains counts. "
        "This is simply required if you use an R method that needs "
        "to build a Seurat object."
        "If not provided, will use AnnData.X which usually works fine"
        "when operating on the dimensionality reduction."
    )
    subparser.add_argument(
        "--data-key",
        nargs=None,
        type=str,
        dest="data_key",
        required=False,
        default=None,
        help="Key in the AnnData.layers that contains data. "
        "This is simply required if you use an R method that needs "
        "to build a Seurat object."
        "If not provided, will use AnnData.X which usually works fine "
        "when operating on the dimensionality reduction."
    )
    subparser.add_argument(
        "--scale-data-key",
        nargs=None,
        type=str,
        dest="scale_data_key",
        required=False,
        default=None,
        help="Key in the AnnData.layers that contains scaled data. "
        "This is simply required if you use an R method that needs "
        "to build a Seurat object."
        "If not provided, will use AnnData.X which usually works fine "
        "when operating on the dimensionality reduction."
    )
    subparser.add_argument(
        "--obsm_key",
        nargs=None,
        type=str,
        dest="obsm_key",
        required=False,
        default=None,
        help="Key in the AnnData.obsm to use as dimensionality reduction for correction."
    )
    subparser.add_argument(
        "--correction-method",
        nargs=None,
        type=str,
        dest="correction_method",
        choices=["none", "harmonyR"],
        required=False,
        default="none",
        help="Integration methods to use for analysis. "
    )
    subparser.add_argument(
        "--n-components",
        type=int,
        default=consts.DEFAULT_N_COMPONENTS,
        dest="n_components",
        help="Number of components to use for kNN graph calculation on corrected dimensionality reduction.",
    )
    subparser.add_argument(
        '--n-neighbors',
        type=int,
        default=consts.DEFAULT_N_NEIGHBORS,
        dest="n_neighbors",
        help="Number of neighbors to use for kNN graph calculation on corrected dimensionality reduction.",
    )
    subparser.add_argument(
        '--clust-resolution',
        type=float,
        default=consts.DEFAULT_CLUST_RESOLUTION,
        dest="clust_resolution",
        help="Resolution for an clustering on corrected dimensionality reduction.",
    )
    subparser.add_argument(
        "--umap-min-distance",
        type=float,
        default=consts.DEFAULT_UMAP_MIN_DISTANCE,
        dest="umap_min_distance",
        help="Minimum distance for UMAP.",
    )
    subparser.add_argument(
        "--plot-keys",
        nargs="+",
        type=str,
        dest="plot_keys",
        required=False,
        default=None,
        help="Additional keys to use for plots.",
    )
    subparser.add_argument(
        "--random-state",
        nargs=None,
        type=int,
        dest="random_state",
        required=False,
        default=consts.RANDOM_STATE,
        help="Random state to use for random number generators.",
    )
    subparser.add_argument(
        "--cpu-threads",
        type=int,
        default=None,
        dest="n_threads",
        help="Number of threads to use when pytorch is run "
        "on CPU. Defaults to the number of logical cores.",
    )
    subparser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        help="Including the flag --debug will log "
        "extra messages useful for debugging.",
    )
    return subparsers
