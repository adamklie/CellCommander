import argparse

from cellcommander.annotate import consts


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for annotate.

    Args:
        subparsers: Parser object before addition of arguments specific to annotate.

    Returns:
        parser: Parser object with additional parameters.

    """

    subparser = subparsers.add_parser(
        "annotate",
        description="Runs a cell annotation pipeline on an input single cell dataset",
        help="Takes in a single file and adds a column to the metadata with a description of cell "
            "identity.",
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
        default="annotate",
        help="Prefix for output files. " "If not provided, the prefix will be 'annotate'.",
    )
    subparser.add_argument(
        "--methods",
        nargs="+",
        type=str,
        dest="methods",
        choices=["manual"],
        required=False,
        default="manual",
        help="Method to use for cell identity annotation. "
        "Currently supports 'manual' annotation. "
        "Future versions will support 'CellTypis' and 'DecoupleR based annotation'.",
    )
    subparser.add_argument(
        "--annotation-key",
        nargs=None,
        type=str,
        dest="annotation_key",
        required=False,
        default="manual_cellid_annotation",
        help="Name of key in the metadata to use for annotation. "
        "If not provided, the default will be 'manual_cellid_annotation'. "
    )
    subparser.add_argument(
        "--layer",
        nargs=None,
        type=str,
        dest="layer",
        required=False,
        default=None,
        help="Layer in the AnnData to use for plotting dimensionality reduction. "
        "If not provided, the '.X' matrix will be used. "
    )
    subparser.add_argument(
        "--dim-reduction",
        nargs=None,
        type=str,
        dest="dim_reduction",
        required=False,
        default="X_umap",
        help="Name of key containing dimensionality reduction to use for plotting. "
        "E.g. for ScanPy, this will often be 'X_umap' in the .obsm dictionary of the AnnData. "
    )
    subparser.add_argument(
        "--marker-gene-list",
        nargs=None,
        type=str,
        dest="marker_gene_list",
        required=False,
        default=None,
        help="Path to marker gene list. "
        "At mininum this should be a single column file with gene names. "
        "If a second column is provided, it will be interpreted as the cell identitiy "
        "associated with the gene. Subsequent plots will be saved separately for each cell identity. "
    )
    subparser.add_argument(
        "--cluster-key",
        nargs=None,
        type=str,
        dest="cluster_key",
        required=False,
        default=None,
        help="Name of key in the metadata containing cluster assignments. "
        "If not provided, clustering will be performed using Leiden clustering "
        "with the parameters provided in --clust-n-neighbors and --clust-resolution. "
    )
    subparser.add_argument(
        '--clust-n-neighbors',
        type=int,
        default=consts.DEFAULT_CLUST_N_NEIGHBORS,
        dest="clust_n_neighbors",
        help="Number of neighbors for clustering if --cluster-key is not provided.",
    )
    subparser.add_argument(
        "--clust-dim-reduction",
        nargs=None,
        type=str,
        dest="clust_dim_reduction",
        required=False,
        default="X_pca",
        help="Name of key containing dimensionality reduction to use for clustering. "
        "E.g. for ScanPy, this will often be 'X_pca' in the .obsm dictionary of the AnnData. "
    )
    subparser.add_argument(
        "--clust-n-components",
        type=int,
        default=consts.DEFAULT_CLUST_N_COMPONENTS,
        dest="clust_n_components",
        help="Number of components for clustering if --cluster-key is not provided.",
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
