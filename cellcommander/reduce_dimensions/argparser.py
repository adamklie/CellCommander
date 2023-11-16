import argparse

from cellcommander.reduce_dimensions import consts


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for reduce-dimensions.

    Args:
        subparsers: Parser object before addition of arguments specific to reduce-dimensions.

    Returns:
        parser: Parser object with additional parameters.

    """

    subparser = subparsers.add_parser(
        "reduce-dimensions",
        description="Performs dimensionality reduction on single-cell data with a variety of methods.",
        help="This command runs an initial dimensionality reduction of choice on a chosen layer of the passed in object. "
            " By default, this command runs a UMAP on that dimensionality reduction, performs an inital clustering with "
            " the Leiden algorithm, and plots the results. It also saves a new object with all this information stored in it.",
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
        default="reduce_dimensions",
        help="Prefix for output files. " "If not provided, the prefix will be 'reduce-dimensions'.",
    )
    subparser.add_argument(
        "--method",
        nargs=None,
        type=str,
        dest="method",
        choices=["scanpy_default", "seurat_default", "seurat_sctransform", "signac_default", "lsi"],
        required=False,
        default="scanpy_default",
        help="Method to use for dimensionality reduction. "
        "By default, uses Scanpy to perform PCA and UMAP. ",
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
        "--obsm_key",
        nargs=None,
        type=str,
        default=None,
        dest="obsm_key",
        required=False,
        help="Key in the .obsm attribute to use for dimensionality reduction. "
        "If specified, the --layer argument will be ignored. "
        "Only applicable to AnnData objects. ",
    )
    subparser.add_argument(
        "--variable-features-key",
        nargs=None,
        type=str,
        dest="variable_features_key",
        required=False,
        default=None,
        help="Key in the variable metadata to use for variable features. "
        "Will subset a copy of the object before performing dimensionality reduction. ",
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
        "--scale-data",
        dest="scale_data",
        action="store_true",
        help="Including the flag --scale-data will scale the data "
        "before performing dimensionality reduction.",
    )
    subparser.add_argument(
        "--scale-max",
        nargs=None,
        type=float,
        dest="scale_max",
        required=False,
        default=consts.DEFAULT_SCALE_MAX,
        help="Maximum value for scaling.",
    )
    subparser.add_argument(
        '--n-neighbors',
        type=int,
        default=consts.DEFAULT_N_NEIGHBORS,
        dest="n_neighbors",
        help="Number of neighbors for kNN graph calculation.",
    )
    subparser.add_argument(
        "--n-components",
        type=int,
        default=consts.DEFAULT_N_COMPONENTS,
        dest="n_components",
        help="Number of components to use for kNN graph calculation.",
    )
    subparser.add_argument(
        "--umap-min-distance",
        type=float,
        default=consts.DEFAULT_UMAP_MIN_DISTANCE,
        dest="umap_min_distance",
        help="Minimum distance for UMAP.",
    )
    subparser.add_argument(
        '--clust-resolution',
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
