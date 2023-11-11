"""Command-line tool functionality for normalize."""

import argparse
import logging
import os
import sys

import cellcommander
from cellcommander.base_cli import AbstractCLI, get_version
from cellcommander.normalize.checkpoint import create_workflow_hashcode


class CLI(AbstractCLI):
    """CLI implements AbstractCLI from the cellcommander package."""

    def __init__(self):
        self.name = "normalize"
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

        # If either log-tf or log-idf is set to True, set log-tf-idf to False
        if args.log_tf or args.log_idf:
            if args.log_tfidf:
                args.log_tfidf = False
                logging.warning(
                    "If either log-tf and log-idf are set to True, "
                    "log-tf-idf cannot be set to True. Setting log-tf-idf to False."
                )
        # If log-tf-idf is set to True, set log-tf and log-idf to False
        if args.log_tfidf:
            args.log_tf = False
            args.log_idf = False
            logging.warning(
                "log-tfidf is set to True. Setting log-tf and log-idf to False."
            )
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
    log_file = os.path.join(output_dir, "normalize.log")
    logger = logging.getLogger("cellcommander")  # name of the logger
    logger.setLevel(logging.INFO if not args.debug else logging.DEBUG)
    formatter = logging.Formatter("cellcommander:normalize: %(message)s")
    file_handler = logging.FileHandler(filename=log_file, mode="w", encoding="UTF-8")
    console_handler = logging.StreamHandler()
    file_handler.setFormatter(formatter)  # set the file format
    console_handler.setFormatter(formatter)  # use the same format for stdout
    logger.addHandler(file_handler)  # log to file
    logger.addHandler(console_handler)  # log to stdout

    # Log the command as typed by user.
    logger.info("Command:\n" + " ".join(["cellcommander", "normalize"] + sys.argv[2:]))
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
    from cellcommander.normalize.run import run_normalize
    run_normalize(args)
    file_handler.close()

    return
