import argparse

from cellcommander.recipes import consts


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for recipes.

    Args:
        subparsers: Parser object before addition of arguments specific to recipes.

    Returns:
        parser: Parser object with additional parameters.

    """

    subparser = subparsers.add_parser(
        "recipes",
        description="Runs method workflows.",
        help="These are chained workflows that give less fine-grained control over the parameters of the methods and"
        "and don't allow for as much interchanging in steps. These are usually run if you want to get a quick "
        "overview of the method behaves on a dataset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparser.add_argument(
        "--input_paths",
        nargs="+",
        type=str,
        dest="input_files",
        required=True,
        help="Data files on which to run tool. "
        "This is largely dependent on the input data, tool and mode "
        "For example, for scATAC-seq data analysis methods, this would often be a fragment file. ",
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
        "--sample_name",
        nargs=None,
        type=str,
        dest="sample_name",
        required=False,
        default=None,
        help="Name of the sample. "
        "Will be prepended to barcode followed by '#' if provided. ",
    )
    subparser.add_argument(
        "--metadata_path",
        nargs=None,
        type=str,
        dest="metadata_path",
        required=False,
        default=None,
        help="Path to a metadata file. "
        "This is a tab-separated file that contains metadata for the cells. "
        "The first column should be the barcodes. Will subset the data to only include cells in the metadata."
        "If not provided, no metadata will be used.",
    )
    subparser.add_argument(
        "--metadata_source",
        nargs=None,
        type=str,
        dest="metadata_source",
        required=False,
        default="external",
        help="Name of the metadata source for the metadata file."
        "This is used as the prefix for the column names "
        "added to the adata.obs dataframe."
        "If not provided, the source will be 'external'.",
    )
    subparser.add_argument(
        "--method",
        nargs=None,
        type=str,
        dest="method",
        choices=["snapatac2"],
        required=True,
        help="Tool to use for analysis. "
        "Currently supported tools are: "
        "snapatac2",
    )
    subparser.add_argument(
        "--mode",
        nargs=None,
        type=str,
        dest="mode",
        choices=["single-sample", "multi-sample"],
        required=True,
        help="Mode to run the tool in. "
        "Currently supported modes are: "
        "single-sample, multi-sample",
    )
    subparser.add_argument(
        "--params_path",
        nargs=None,
        type=str,
        dest="params_file",
        required=False,
        default=None,
        help="Path to a params file. "
        "This is a YAML file that contains parameters for the tool. "
        "If not provided, default parameters will be used.",
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
