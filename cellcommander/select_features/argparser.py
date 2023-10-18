import argparse

from cellcommander.select_features import consts


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for select-features.

    Args:
        subparsers: Parser object before addition of arguments specific to select-features.

    Returns:
        parser: Parser object with additional parameters.

    """

    subparser = subparsers.add_parser(
        "select-features",
        description="Selects features from single-cell matrices with a variety of methods.",
        help="Run feature selection on single-cell data. ",
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
        help="Output directory location. " "If it does not exist, it will be created.",
    )
    subparser.add_argument(
        "--output_prefix",
        nargs=None,
        type=str,
        dest="output_prefix",
        required=False,
        default="select-features",
        help="Prefix for output files. " "If not provided, the prefix will be 'select-features'.",
    )
    subparser.add_argument(
        "--methods",
        nargs="+",
        type=str,
        dest="method",
        choices=["seurat", "deviance"],
        required=False,
        default="seurat",
        help="Method to use for normalization. "
        "Currently supports 'seurat' for scRNA-seq data "
    )
    subparser.add_argument(
        '--initial-clust-n-neighbors',
        type=int,
        default=consts.DEFAULT_INITIAL_CLUST_N_NEIGHBORS,
        dest="initial_clust_n_neighbors",
        help="Number of neighbors for initial clustering.",
    )
    subparser.add_argument(
        '--initial_clust_resolution',
        type=float,
        default=consts.DEFAULT_INITIAL_CLUST_RESOLUTION,
        dest="initial_clust_resolution",
        help="Resolution for initial clustering.",
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
