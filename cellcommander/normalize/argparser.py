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
        "--save-normalized-mtx",
        dest="save_normalized_mtx",
        action="store_true",
        help="Including the flag --save-normalized-matrix will save the normalized matrix to a file, "
        "in addition to a barcodes.tsv and features.tsv file.",
    )
    subparser.add_argument(
        "--methods",
        nargs="+",
        type=str,
        dest="methods",
        choices=["log1p", "sctransform", "tfidf"],
        required=False,
        default="log1p",
        help="Method to use for normalization. "
        "Currently supports 'log1p', 'sctransform' for scRNA-seq data "
        "and 'lsi' for scATAC-seq data. ",
    )
    subparser.add_argument(
        "--filter_features",
        action="store_true",
        dest="filter_features",
        required=False,
        default=False,
        help="Whether to filter features before normalization. "
        "This currently only applies to the 'sctransform' method. "
        "and filters out genes expressed in less than 20 cells"
    )
    subparser.add_argument(
        "--vars_to_regress",
        nargs="+",
        type=str,
        dest="vars_to_regress",
        required=False,
        default=None,
        help="Variables to regress out of the data. "
        "This currently only applies to the 'sctransform' method."
    )
    subparser.add_argument(
        "--log_idf",
        nargs=None,
        type=bool,
        dest="log_idf",
        required=False,
        default=True,
        help="Whether to log the inverse document frequency values for the tfidf method. "
        "To run the Signac version of tfidf "
        "set this to False, '--log-tf' to False, and '--log-tf-idf to True.",
    )
    subparser.add_argument(
        "--log-tf",
        nargs=None,
        type=bool,
        dest="log_tf",
        required=False,
        default=True,
        help="Whether to log the term frequency values for the tfidf method. "
        "To run the Signac version of tfidf "
        "set this to False, '--log-idf' to False, and '--log-tf-idf to True.",
    )
    subparser.add_argument(
        "--log-tf-idf",
        nargs=None,
        type=bool,
        dest="log_tfidf",
        required=False,
        default=False,
        help="Whether to log the tfidf values for the tfidf method. "
        "To run the Signac version of tfidf "
        "If True, 'log-tf' and 'log-idf' are automatically set to False.",
    )
    subparser.add_argument(
        "--tfidf-scale-factor",
        nargs=None,
        type=float,
        dest="tfidf_scale_factor",
        required=False,
        default=1e4,
        help="Scale factor for tfidf normalization.",
    )
    subparser.add_argument(
        '--clust-n-neighbors',
        type=int,
        default=consts.DEFAULT_CLUST_N_NEIGHBORS,
        dest="clust_n_neighbors",
        help="Number of neighbors for clustering.",
    )
    subparser.add_argument(
        '--clust_resolution',
        type=float,
        default=consts.DEFAULT_CLUST_RESOLUTION,
        dest="clust_resolution",
        help="Resolution for clustering.",
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
