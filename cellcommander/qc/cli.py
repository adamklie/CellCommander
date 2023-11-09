"""Command-line tool functionality for qc."""

import argparse
import logging
import os
import sys

import cellcommander
from cellcommander.base_cli import AbstractCLI, get_version
from cellcommander.qc.checkpoint import create_workflow_hashcode
from cellcommander.qc.run import run_qc


class CLI(AbstractCLI):
    """CLI implements AbstractCLI from the cellcommander package."""

    def __init__(self):
        self.name = "qc"
        self.args = None

    def get_name(self) -> str:
        return self.name

    @staticmethod
    def validate_args(args) -> argparse.Namespace:
        """Validate parsed arguments."""

        # Ensure that if there's a tilde for $HOME in the file path, it works.
        try:
            args.input_file = os.path.expanduser(args.input_file)
            args.output_dir = os.path.expanduser(args.output_dir)
            if args.metadata_file:
                args.metadata_file = os.path.expanduser(args.metadata_file)
        except TypeError:
            raise ValueError("Problem with provided input and output paths.")

        # Ensure write access to the save directory.
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)
        if args.output_dir:
            assert os.access(args.output_dir, os.W_OK), (
                f"Cannot write to specified output directory {args.output_dir}. "
                f"Ensure the directory exists and is write accessible."
            )
        # If filtering_strategy is "mad", make all thresholds None
        if args.filtering_strategy == "mad":
            args.n_features_low = None
            args.n_features_hi = None
            args.total_counts_low = None
            args.total_counts_hi = None
            args.ns_hi = None
            args.tss_low = None
            args.tss_hi = None

        # Do the opposite for nmads if the strategy is "threshold"
        elif args.filtering_strategy == "threshold":
            args.total_counts_nmads = None
            args.n_features_nmads = None
            args.pct_counts_in_top_features_nmads = None
            args.pct_counts_mt_nmads = None
            args.ns_nmads = None
            args.tss_nmads = None

        # Validate numerical arguments that should have logical minimums
        if args.total_counts_nmads is not None:
            assert (
                args.total_counts_nmads > 0
            ), "--total_counts_nmads must be a positive number."

        if args.n_features_nmads is not None:
            assert (
                args.n_features_nmads > 0
            ), "--n_features_nmads must be a positive number."

        if args.n_top_features is not None:
            assert (
                args.n_top_features > 0
            ), "--n_top_features must be an integer greater than 0."

        if args.pct_counts_in_top_features_nmads is not None:
            assert (
                args.pct_counts_in_top_features_nmads > 0
            ), "--pct_counts_in_top_features_nmads must be a positive number."

        if args.pct_counts_mt_nmads is not None:
            assert (
                args.pct_counts_mt_nmads > 0
            ), "--pct_counts_mt_nmads must be a positive number."

        # Validate that percentage values are within a logical range
        if args.pct_counts_mt_hi is not None:
            assert (
                0 <= args.pct_counts_mt_hi <= 100
            ), "--pct_counts_mt_hi must be a percentage between 0 and 100."

        # Validate that n_features_low and n_features_hi are logical
        if args.n_features_low is not None and args.n_features_hi is not None:
            assert (
                args.n_features_low < args.n_features_hi
            ), "--n_features_low should be less than --n_features_hi."
            assert (
                args.n_features_low > 0
            ), "--n_features_low must be a positive integer."
            assert args.n_features_hi > 0, "--n_features_hi must be a positive integer."

        # Validate that min_cells_per_feature is logical
        if args.min_cells_per_feature is not None:
            assert (
                args.min_cells_per_feature > 0
            ), "--min_cells_per_feature must be an integer greater than 0."

        # Make sure n_threads makes sense.
        if args.n_threads is not None:
            assert args.n_threads > 0, "--cpu-threads must be an integer >= 1"

        # Return the validated arguments.
        return args

    @staticmethod
    def run(args):
        """Run the main tool functionality on parsed arguments."""

        # Run the tool.
        return main(args)


def setup_and_logging(args):
    """Take command-line input, parse arguments, and run tests or tool."""

    # Send logging messages to stdout as well as a log file.
    output_dir = args.output_dir
    log_file = os.path.join(output_dir, "qc.log")
    logger = logging.getLogger("cellcommander")  # name of the logger
    logger.setLevel(logging.INFO if not args.debug else logging.DEBUG)
    formatter = logging.Formatter("cellcommander:qc: %(message)s")
    file_handler = logging.FileHandler(filename=log_file, mode="w", encoding="UTF-8")
    console_handler = logging.StreamHandler()
    file_handler.setFormatter(formatter)  # set the file format
    console_handler.setFormatter(formatter)  # use the same format for stdout
    logger.addHandler(file_handler)  # log to file
    logger.addHandler(console_handler)  # log to stdout

    # Log the command as typed by user.
    logger.info("Command:\n" + " ".join(["cellcommander", "qc"] + sys.argv[2:]))
    logger.info("cellcommander " + get_version())

    # Set up checkpointing by creating a unique workflow hash.
    hashcode = create_workflow_hashcode(
        module_path=os.path.dirname(cellcommander.__file__),
        args_to_remove=(
            [
                "output_dir",
                "debug",
                "cpu_threads",
            ]
        ),
        args=args,
    )
    args.checkpoint_filename = hashcode  # store this in args
    logger.info(f"(Workflow hash {hashcode})")
    return args, file_handler


def main(args):
    """Take command-line input, parse arguments, and run tests or tool."""

    args, file_handler = setup_and_logging(args)

    # Run the tool.
    run_qc(args)
    file_handler.close()

    return
