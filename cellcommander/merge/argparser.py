import argparse

from cellcommander.merge import consts


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for merge.

    Args:
        subparsers: Parser object before addition of arguments specific to merge.

    Returns:
        parser: Parser object with additional parameters.

    """

    subparser = subparsers.add_parser(
        "merge",
        description="Merges multiple objects together.",
        help="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparser.add_argument(
        "--input_paths",
        nargs="+",
        type=str,
        dest="input_files",
        required=True,
        help="Data files on which load and merge.",
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
        default="merge",
        help="Prefix for output files. " "If not provided, the prefix will be 'reduce-dimensions'.",
    )
    subparser.add_argument(
        "--names",
        nargs="+",
        type=str,
        dest="names",
        required=False,
        default=None,
        help="Name of each file to use for integration. "
        "If not provided, the sample names will be "
        "extracted from the input files. "
    )
    subparser.add_argument(
        "--names-key",
        nargs=None,
        type=str,
        dest="names_key",
        required=False,
        default="sample",
        help="Name of the key in the AnnData.obs to use for integration. ",
    )
    subparser.add_argument(
        "--layer",
        nargs=None,
        type=str,
        dest="layer",
        required=False,
        default=None,
        help="Layer in the AnnData to use for dimensionality reduction. "
        "If not provided, the '.X' matrix will be used. "
        "Note that choice of layer should be consistent with the choice of method. "
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
