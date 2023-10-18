import argparse

from cellcommander.normalize import consts


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for qc.

    Args:
        subparsers: Parser object before addition of arguments specific to qc.

    Returns:
        parser: Parser object with additional parameters.

    """

    subparser = subparsers.add_parser(
        "normalize",
        description="Normalize single-cell data with a variety of methods.",
        help="Runs normalization methods. ",
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
        default="normalize",
        help="Prefix for output files. " "If not provided, the prefix will be 'normalize'.",
    )
    subparser.add_argument(
        "--methods",
        nargs="+",
        type=str,
        dest="method",
        choices=["log1p", "scran", "pearson", "depth", "sctransform", "lsi"],
        required=False,
        default="log1p",
        help="Method to use for normalization. "
        "Currently supports 'log1p', 'scran', 'pearson', 'depth', 'sctransform' for scRNA-seq data "
        "and 'lsi' for scATAC-seq data. ",
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
