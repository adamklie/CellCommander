import argparse

from cellcommander.detect_doublets import consts


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for detect_doublets.

    Args:
        subparsers: Parser object before addition of arguments specific to qc.

    Returns:
        parser: Parser object with additional parameters.

    """

    subparser = subparsers.add_parser(
        "detect-doublets",
        description="Single sample doublet detection",
        help="Runs doublet detection using one or more methods. "
        "Currently supports running Scrublet and scDblFinder "
        "for either scRNA-seq data scATAC-seq data, and "
        "AMULET for scATAC-seq data. ",
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
        default="detect_doublets",
        help="Prefix for output files. " "If not provided, the prefix will be 'detect_doublets'.",
    )
    subparser.add_argument(
        "--method",
        nargs=None,
        type=str,
        dest="method",
        choices=["scrublet", "scDblFinder", "amulet", "cellranger", "consensus"],
        required=False,
        default="scrublet",
        help="Mode of the data. "
        "Currently supports 'scrublet', 'scDblFinder', 'cellranger', 'amulet' or 'consenesus`. "
        "If 'amulet' is selected, the data must be scATAC-seq. "
        "If cellranger is selected, the the column 'excluded_reason_cellranger' must be present in adata.obs. "
        "If 'consensus' is selected, all methods will be run doublets will be removed "
        "based on the --consensus-method argument. ",
    )
    subparser.add_argument(
        "--consensus-methods",
        nargs='+',
        type=str,
        dest="consensus_methods",
        required=False,
        default=None,
        help="Methods to include for consensus doublet detection. "
        "'method' argument must be set to 'consensus'."
        "Currently supports 'scrublet', 'scDblFinder', 'amulet' and 'cellranger'. "
    )
    subparser.add_argument(
        "--consensus-strategy",
        nargs=None,
        type=str,
        dest="consensus_strategy",
        choices=["union", "intersection", "majority"],
        required=False,
        default="union",
        help="Method for consensus doublet detection. "
        "Currently supports 'union' and 'intersection'. "
        "If 'union' is selected, doublets will be removed if they are predicted to be doublets "
        "by any of the methods. "
        "If 'intersection' is selected, doublets will be removed if they are predicted to be doublets "
        "by all of the methods. "
        "If 'majority' is selected, doublets will be removed if they are predicted to be doublets "
        "by a majority of the methods. ",
    )
    subparser.add_argument(
        "--scDblFinder-atac-params",
        action="store_true",
        dest="scDblFinder_atac_params",
        help="If included, the scDblFinder parameters will be set to the defaults for scATAC-seq data. "
        "See https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scATAC.html#the-clamulet-method "
        "for more details. ",
    )
    subparser.add_argument(
        "--no-filter",
        action="store_true",
        dest="no_filter",
        required=False,
        default=False,
        help="If included, Only the doublet metrics will be calculated "
        "and the saved AnnData object will not be filtered. "
        "By default, the AnnData object will be filtered. ",
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
