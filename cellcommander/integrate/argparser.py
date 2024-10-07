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
        "--batch_file",
        nargs=None,
        type=str,
        dest="batch_file",
        required=False,
        default=None,
        help="File containing batch information for each barcode. "
        "Expected format is a tab-separated file with the first column as the barcode and the remaining columns as batch information.",
    )
    subparser.add_argument(
        "--vars_to_correct",
        nargs="+",
        type=str,
        dest="vars_to_correct",
        required=False,
        default=None,
        help="Variables to correct for."
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
        "--components-to-remove",
        nargs="+",
        type=int,
        dest="components_to_remove",
        required=False,
        default=None,
        help="Number of components to remove from the dimensionality reduction. "
        "This is often utilized for scATAC-seq data, where the first component is "
        "can be correlated with the number of fragments. "
    )
    subparser.add_argument(
        "--method",
        nargs=None,
        type=str,
        dest="method",
        choices=["harmony", "harmonyR"],
        required=False,
        default="none",
        help="Integration methods to use for analysis. "
    )
    subparser.add_argument(
        "--corrected_obsm_key",
        nargs=None,
        type=str,
        dest="corrected_obsm_key",
        required=False,
        default=None,
        help="Key in the AnnData.obsm to store corrected dimensionality reduction."
    )
    subparser.add_argument(
        "--max-iter-harmony",
        type=int,
        default=consts.DEFAULT_MAX_ITER_HARMONY,
        dest="max_iter_harmony",
        help="Maximum number of iterations to run harmony correction.",
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
