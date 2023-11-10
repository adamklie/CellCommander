"""Single run of summarize, given input arguments."""

import argparse
import logging
import os
import glob
import sys
import yaml
import traceback
from datetime import datetime
from typing import Dict, Optional, Tuple, Union

import numpy as np
import matplotlib
import muon as mu
import pandas as pd
import psutil
import scanpy as sc
from mudata import MuData
from anndata import AnnData

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')
import seaborn as sns

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)
import cellcommander
from cellcommander.utils import describe_anndata
from cellcommander.summarize import consts
from cellcommander.summarize import rna

logger = logging.getLogger("cellcommander")


def run_summarize(args: argparse.Namespace):
    """The full script for the command line tool to perform summarize and filtering.

    Args:
        args: Inputs from the command line, already parsed using argparse.

    Note: Returns nothing, but writes output to a file(s) specified from command line.

    """

    # Log the start time.
    logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    logger.info("Running summarize command")

    try:
        # Grab all subcommand directories
        subcommand_dirs = []
        for input_dir in args.input_dirs:
            logger.info("Looking for subcommands in {}".format(input_dir))
            subcommand_dirs += glob.glob(os.path.join(input_dir, "*"))
        subcommand_dirs = [x for x in subcommand_dirs if os.path.isdir(x)]

        # Remove ignored subcommands
        if args.ignore:
            logger.info("Ignoring subcommands: {}".format(args.ignore))
            subcommand_dirs = [x for x in subcommand_dirs if os.path.basename(x) not in args.ignore]

        # Subcommands kept
        subcommands = [os.path.basename(x) for x in subcommand_dirs]
        
        # Set up markdown text
        markdown_text = "# CellCommander summary\n"

        # Add date and version
        markdown_text += "Date: {}<br>".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        markdown_text += "CellCommander version: {}\n\n".format(cellcommander.__version__)

        # Add summary as an <aside> tag
        markdown_text += f"<aside class='summary'>{args.summary}</aside>\n\n"

        if "qc" in subcommands:
            # Get subcommand dir from index of the subcoomand
            subcommand_dir = subcommand_dirs[subcommands.index("qc")]

            # Add qc heading
            markdown_text += "# `qc`\n\n"

            # List all the files in the `qc` directory
            qc_files = glob.glob(os.path.join(subcommand_dir, "*"))

            # Open up the file with *.log
            log_file = [x for x in qc_files if x.endswith(".log")][0]
            log_data = open(log_file).readlines()
            
            # Write to markdown
            markdown_text = rna.write_qc_log_data(log_data, markdown_text)
            markdown_text = rna.add_qc_images(subcommand_dir, args.output_dir, markdown_text)

            # Open up the anndata file 
            anndata_file = [x for x in qc_files if x.endswith(".h5ad")][0]
            adata = sc.read_h5ad(anndata_file)
            adata_obs_head = adata.obs.head()
            markdown_text += "The first 5 rows of the `adata.obs` dataframe:\n\n"
            markdown_text += adata_obs_head.to_markdown() + "\n\n"

        if "remove_background" in subcommands:
            # Get subcommand dir from index of the subcoomand
            subcommand_dir = subcommand_dirs[subcommands.index("remove_background")]

            # Add remove_background heading
            markdown_text += "# `remove_background`\n\n"

            # List all the files in the `remove_background` directory
            remove_background_files = glob.glob(os.path.join(subcommand_dir, "*"))

            # Open up the file with *.log
            log_file = [x for x in remove_background_files if x.endswith(".log")][0]
            log_data = open(log_file).readlines()

            # Write to markdown
            markdown_text = rna.write_remove_background_log_data(log_data, markdown_text)
            markdown_text = rna.add_remove_background_images(subcommand_dir, args.output_dir, markdown_text)

        if "doublet_detection" in subcommands:
            # Get subcommand dir from index of the subcoomand
            subcommand_dir = subcommand_dirs[subcommands.index("doublet_detection")]
            
            # Add doublet_detection heading
            markdown_text += "# `doublet_detection`\n\n"

            # List all the files in the `doublet_detection` directory
            doublet_detection_files = glob.glob(os.path.join(subcommand_dir, "*"))

            # Open up the file with *.log
            log_file = [x for x in doublet_detection_files if x.endswith(".log")][0]
            log_data = open(log_file).readlines()

            # Write to markdown
            markdown_text = rna.write_doublet_detection_log_data(log_data, markdown_text)
            markdown_text = rna.add_doublet_detection_images(subcommand_dir, args.output_dir, markdown_text)

        if "normalization" in subcommands:
            # Get subcommand dir from index of the subcoomand
            subcommand_dir = subcommand_dirs[subcommands.index("normalization")]

            # Add normalization heading
            markdown_text += "# `normalization`\n\n"

            # List all the files in the `normalization` directory
            normalization_files = glob.glob(os.path.join(subcommand_dir, "*"))

            # Open up the file with *.log
            log_file = [x for x in normalization_files if x.endswith(".log")][0]
            log_data = open(log_file).readlines()

            # Write to markdown
            markdown_text = rna.write_normalization_log_data(log_data, markdown_text)
            markdown_text = rna.add_normalization_images(subcommand_dir, args.output_dir, markdown_text)

        if "feature_selection" in subcommands:
            # Get subcommand dir from index of the subcoomand
            subcommand_dir = subcommand_dirs[subcommands.index("feature_selection")]

            # Add feature_selection heading
            markdown_text += "# `feature_selection`\n\n"

            # List all the files in the `feature_selection` directory
            feature_selection_files = glob.glob(os.path.join(subcommand_dir, "*"))

            # Open up the file with *.log
            log_file = [x for x in feature_selection_files if x.endswith(".log")][0]
            log_data = open(log_file).readlines()

            # Write to markdown
            markdown_text = rna.write_feature_selection_log_data(log_data, markdown_text)
            markdown_text = rna.add_feature_selection_images(subcommand_dir, args.output_dir, markdown_text)

        if "dimensionality_reduction" in subcommands:
            # Get subcommand dir from index of the subcoomand
            subcommand_dir = subcommand_dirs[subcommands.index("dimensionality_reduction")]

            # Add dimensionality_reduction heading
            markdown_text += "# `dimensionality_reduction`\n\n"

            # List all the files in the `dimensionality_reduction` directory
            dimensionality_reduction_files = glob.glob(os.path.join(subcommand_dir, "*"))

            # Open up the file with *.log
            log_file = [x for x in dimensionality_reduction_files if x.endswith(".log")][0]
            log_data = open(log_file).readlines()

            # Write to markdown
            markdown_text = rna.write_dimensionality_reduction_log_data(log_data, markdown_text)
            markdown_text = rna.add_dimensionality_reduction_images(subcommand_dir, args.output_dir, markdown_text)

        if "annotate" in subcommands:
            # Get subcommand dir from index of the subcoomand
            subcommand_dir = subcommand_dirs[subcommands.index("annotate")]

            # Add annotate heading
            markdown_text += "# `annotate`\n\n"

            # List all the files in the `annotate` directory
            annotate_files = glob.glob(os.path.join(subcommand_dir, "*"))

            # Open up the file with *.log
            #log_file = [x for x in annotate_files if x.endswith(".log")][0]
            #log_data = open(log_file).readlines()
            
            # Write to markdown
            markdown_text = rna.write_annotate_log_data(markdown_text)
            markdown_text = rna.add_annotate_images(subcommand_dir, args.output_dir, markdown_text)

            # Open up the anndata file, and add the string representation of it to the markdown
            anndata_file = [x for x in annotate_files if x.endswith(".h5ad")][0]
            adata = sc.read_h5ad(anndata_file)
            markdown_text += adata.__repr__() + "\n\n"
            
        # Write markdown text to file
        with open(os.path.join(args.output_dir, "summary.md"), "w") as f:
            f.write(markdown_text)      

        # Log the end time
        logger.info("Completed summarize")
        logger.info(datetime.now().strftime("%Y-%m-%d %H:%M:%S\n"))

    # The exception allows user to end inference prematurely with CTRL-C.
    except KeyboardInterrupt:
        # If partial output has been saved, delete it.
        full_file = args.output_dir + args.output_prefix + ".h5ad"

        # Name of the filtered (cells only) file.
        if os.path.exists(full_file):
            os.remove(full_file)

        logger.info("Keyboard interrupt. Terminated without saving\n")
