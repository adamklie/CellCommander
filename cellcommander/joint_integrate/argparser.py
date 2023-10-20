import argparse

from cellcommander.joint_integrate import consts


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for joint-integrate.

    Args:
        subparsers: Parser object before addition of arguments specific to joint-integrate.

    Returns:
        parser: Parser object with additional parameters.

    """

    subparser = subparsers.add_parser(
        "joint-integrate",
        description="Runs joint integration between multiple data modalities",
        help="Currently configured to integrate RNA and ATAC objects. "
        "Future iterations will take in MuData objects directly and make less assumptions about the data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparser.add_argument(
        "--rna_h5ad_path",
        nargs=None,
        type=str,
        dest="rna_h5ad_path",
        required=True,
        help="RNA data file on which to run tool. "
        "The following input formats are supported: "
        "AnnData (.h5ad). "
    )
    subparser.add_argument(
        "--atac_h5ad_path",
        nargs=None,
        type=str,
        dest="atac_h5ad_path",
        required=True,
        help="ATAC data file on which to run tool. "
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
        default="joint_integrate",
        help="Prefix for output files. " "If not provided, the prefix will be 'joint_integrate'.",
    )
    subparser.add_argument(
        "--methods",
        nargs="+",
        type=str,
        dest="methods",
        required=False,
        default="wnn",
        help="Methods to use for cell identity annotation. "
        "Currently supports 'wnn' which uses Seurat's Weighted Nearest Neighbor algorithm.",
    )
    subparser.add_argument(
        "--rna-dim-reduction",
        nargs=None,
        type=str,
        dest="rna_dim_reduction",
        required=False,
        default="X_pca",
        help="Name of key containing dimensionality reduction to use for RNA modality. "
        "E.g. for ScanPy, this will often be 'X_pca' in the .obsm dictionary of the AnnData. "
    )
    subparser.add_argument(
        "--atac-dim-reduction",
        nargs=None,
        type=str,
        dest="atac_dim_reduction",
        required=False,
        default="X_lsi",
        help="Name of key containing dimensionality reduction to use for ATAC modality. "
        "E.g. for ScanPy, this will often be 'X_lsi' in the .obsm dictionary of the AnnData. "
    )
    subparser.add_argument(
        "--cluster-key",
        nargs=None,
        type=str,
        dest="cluster_key",
        required=False,
        default=None,
        help="If provided, will also cluster the joint data using the Leiden algorithm. "
        "Will then store the cluster labels in the .obs of the AnnData under this key.",
    )
    subparser.add_argument(
        "--clust-resolution",
        type=float,
        default=consts.DEFAULT_CLUST_RESOLUTION,
        dest="clust_resolution",
        help="Resolution for initial clustering if --cluster-key is not provided.",
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
