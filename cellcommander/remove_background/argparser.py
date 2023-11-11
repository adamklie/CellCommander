import argparse

from cellcommander.remove_background import consts


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for remove-background.

    Args:
        subparsers: Parser object before addition of arguments specific to remove-background.

    Returns:
        parser: Parser object with additional parameters.

    """

    subparser = subparsers.add_parser(
        "remove-background",
        description="Remove background from scRNA-seq data.",
        help="Run ambient and background RNA correction with SoupX.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparser.add_argument(
        "--input_h5ad_path",
        nargs=None,
        type=str,
        dest="input_file",
        required=True,
        help="Data file on which to run background removal. "
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
        default="remove_background",
        help="Prefix for output files. " 
        "If not provided, the prefix will be 'remove_background'.",
    )
    subparser.add_argument(
        "--method",
        nargs=None,
        type=str,
        dest="method",
        choices=["soupx", "cellbender"],
        required=False,
        default="soupx",
        help="Method to use for background removal. "
        "Currently, only soupx is supported. ",
    )
    subparser.add_argument(
        "--raw-h5-path",
        nargs=None,
        type=str,
        dest="raw_h5_path",
        required=False,
        default=None,
        help="Path to a raw data file. "
        "Needed for SoupX correction. ",
    )
    subparser.add_argument(
        "--markers_path",
        nargs=None,
        type=str,
        dest="markers_path",
        required=False,
        default=None,
        help="Path to a file containing a list of marker genes. "
        "If provided, the list of marker genes will be used to "
        "by soupx to remove background. "
    )
    subparser.add_argument(
        "--layer",
        nargs=None,
        type=str,
        dest="layer",
        required=False,
        default="corrected_counts",
        help="Layer in the AnnData to add removed counts to. "
        "If not provided, the counts will be added to the "
        "corrected_counts layer. ",
    )
    subparser.add_argument(
        "--initial-clust-num-hvgs",
        nargs=None,
        type=str,
        dest="initial_clust_num_hvgs",
        required=False,
        default=consts.DEFAULT_INITIAL_NUM_HVGS,
        help="Number of highly variable genes to use for "
        "initial clustering. If not provided, 3000 "
        "highly variable genes will be used. ",
    )
    subparser.add_argument(
        '--initial-clust-n-neighbors',
        type=int,
        default=consts.DEFAULT_INITIAL_N_NEIGHBORS,
        dest="initial_clust_n_neighbors",
        help="Number of neighbors for kNN graph calculation.",
    )
    subparser.add_argument(
        "--initial-clust-n-components",
        type=int,
        default=consts.DEFAULT_INITIAL_N_COMPONENTS,
        dest="initial_clust_n_components",
        help="Number of components to use for kNN graph calculation.",
    )
    subparser.add_argument(
        '--initial-clust-resolution',
        type=float,
        default=consts.DEFAULT_INITIAL_CLUST_RESOLUTION,
        dest="initial_clust_resolution",
        help="Resolution for an initial clustering.",
    )
    subparser.add_argument(
        "--umap-min-distance",
        type=float,
        default=consts.DEFAULT_UMAP_MIN_DISTANCE,
        dest="umap_min_distance",
        help="Minimum distance for UMAP.",
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
