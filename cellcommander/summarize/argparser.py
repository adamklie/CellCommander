
import argparse

from cellcommander.summarize import consts


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for summarize.

    Args:
        subparsers: Parser object before addition of arguments specific to summarize.

    Returns:
        parser: Parser object with additional parameters.

    """

    subparser = subparsers.add_parser(
        "summarize",
        description="Summarizes different outputs from workflows.",
        help="Summarizes different outputs from workflows.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparser.add_argument(
        "--input_dirs",
        nargs="+",
        type=str,
        dest="input_dirs",
        required=True,
        help="Directories in which previous CellCommander command outputs were saved.",
    )
    subparser.add_argument(
        "--outdir_path",
        nargs=None,
        type=str,
        dest="output_dir",
        required=False,
        default=None,
        help="Output directory location. " 
        "If it does not exist, it will be created. "
        "If not provided, the output directory will be the same as the input directory.",
    )
    subparser.add_argument(
        "--summary",
        nargs=None,
        type=str,
        dest="summary",
        required=False,
        default="CellCommander summary run",
        help="Summary of the experiment. "
    )
    subparser.add_argument(
        "--ignore",
        nargs="+",
        type=str,
        dest="ignore",
        required=False,
        default=None,
        help="List of paths to ignore when summarizing. "
        "These must be relative paths to the input directory. "
        "For example, if you want to ignore the file "
        "input_dir/qc, "
        "then you would pass in the argument "
        "--ignore qc",
    )
    subparser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        help="Including the flag --debug will log "
        "extra messages useful for debugging.",
    )
    return subparsers
