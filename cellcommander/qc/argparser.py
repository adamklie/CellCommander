import argparse

from cellcommander.qc import consts


def add_subparser_args(subparsers: argparse) -> argparse:
    """Add tool-specific arguments for qc.

    Args:
        subparsers: Parser object before addition of arguments specific to qc.

    Returns:
        parser: Parser object with additional parameters.

    """

    subparser = subparsers.add_parser(
        "qc",
        description="Single sample quality control and filtering",
        help="Runs quality control and filtering on a single sample. "
        "Currently supports scRNA-seq data ('rna' mode) scATAC-seq data ('atac' mode).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparser.add_argument(
        "--input_h5_path",
        nargs=None,
        type=str,
        dest="input_file",
        required=True,
        help="Data file on which to run tool. "
        "The following input formats are supported: "
        "CellRanger v2 and v3 (.h5) or AnnData (.h5ad). "
        "This will work for 10x multiome data as well.",
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
        default="qc",
        help="Prefix for output files. " "If not provided, the prefix will be 'qc'.",
    )
    subparser.add_argument(
        "--metadata_path",
        nargs=None,
        type=str,
        dest="metadata_file",
        required=False,
        default=None,
        help="Path to metadata file "
        "containing information about cell barcodes. "
        "This is most commonly the per_barcode_metrics.csv "
        "file output by CellRanger, but can be any csv file "
        "where the first column is the cell barcode",
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
        "--barcode_exclusion_list_paths",
        nargs="+",
        type=str,
        dest="barcode_exclusion_list_paths",
        required=False,
        default=None,
        help="Paths to files containing barcodes (one per line) to exclude. "
        "If provided, all barcodes in these files will be "
        "excluded from the analysis (e.g. detected doublets). ",
    )
    subparser.add_argument(
        "--mode",
        nargs=None,
        type=str,
        dest="mode",
        required=False,
        choices=["rna", "atac"],
        default="rna",
        help="Mode of the data. "
        "Currently only supports 'rna' and 'atac' mode. "
        "If 'rna' mode is selected, the features for "
        "calculating qcs and thresholds are assumed to be "
        "gene counts."
        "If 'atac' mode is selected, the features for "
        "calculating qcs and thresholds are assumed to be "
        "fragment counts. ",
    )
    subparser.add_argument(
        "--no-filter",
        action="store_true",
        dest="no_filter",
        required=False,
        default=False,
        help="If included, Only the qc metrics will be calculated "
        "and the saved AnnData object will not be filtered. "
        "By default, the AnnData object will be filtered. ",
    )
    subparser.add_argument(
        "--filtering_strategy",
        nargs=None,
        type=str,
        dest="filtering_strategy",
        required=False,
        choices=["mad", "threshold"],
        default="mad",
        help="Filtering strategy to use. "
        "Currently supports either mean absolute deviation (mad) "
        "or user-defined thresholds. "
        "If mad is selected, the thresholds for filtering "
        "are calculated based on the median absolute deviation "
        "of each metric. "
        "If thresholds are selected, the thresholds for filtering "
        "must be provided by the user as arguments.",
    )
    subparser.add_argument(
        "--total_counts_nmads",
        nargs=None,
        type=float,
        dest="total_counts_nmads",
        required=False,
        default=consts.DEFAULT_TOTAL_COUNTS_NMADS,
        help="Number of median absolute deviations "
        "from the median total counts to use as the cutoff "
        "for filtering cells based on total counts. ",
    )
    subparser.add_argument(
        "--n_features_nmads",
        nargs=None,
        type=float,
        dest="n_features_nmads",
        required=False,
        default=consts.DEFAULT_N_FEATURES_NMADS,
        help="Number of median absolute deviations "
        "from the median number of features by counts "
        "to use as the cutoff for filtering cells based "
        "on number of features by counts.",
    )
    subparser.add_argument(
        "--n_top_features",
        nargs=None,
        type=int,
        dest="n_top_features",
        required=False,
        default=20,
        help="Number of top features to calculate the "
        "percentage of total counts in. "
        "Only used if filtering strategy is `mad`.",
    )
    subparser.add_argument(
        "--pct_counts_in_top_features_nmads",
        nargs=None,
        type=int,
        dest="pct_counts_in_top_features_nmads",
        required=False,
        default=consts.DEFAULT_PCT_COUNTS_IN_TOP_FEATURES_NMADS,
        help="All cells with that have a mad for percentage of total counts "
        "in top features greater than this threshold "
        "will be filtered out. "
        "Only used if filtering strategy is `mad`.",
    )
    subparser.add_argument(
        "--pct_counts_mt_nmads",
        nargs=None,
        type=int,
        dest="pct_counts_mt_nmads",
        required=False,
        default=consts.DEFAULT_PCT_COUNTS_MT_NMADS,
        help="All cells with that have a mad for percentage of counts "
        "in mitochondrial features greater than this threshold "
        "will be filtered out. "
        "Only used in 'rna' mode if filtering strategy is `mad`.",
    )
    subparser.add_argument(
        "--ns_nmads",
        nargs=None,
        type=int,
        dest="ns_nmads",
        required=False,
        default=consts.DEFAULT_NS_NMADS,
        help="All cells with that have a mad for nucleosome signal "
        "greater than this threshold "
        "will be filtered out. "
        "Only used in 'atac' mode if filtering strategy is `mad`.",
    )
    subparser.add_argument(
        "--tss_nmads",
        nargs=None,
        type=int,
        dest="tss_nmads",
        required=False,
        default=consts.DEFAULT_TSS_NMADS,
        help="All cells with that have a mad for TSS enrichment "
        "greater than this threshold "
        "will be filtered out. "
        "Only used in 'atac' mode if filtering strategy is `mad`.",
    )
    subparser.add_argument(
        "--n_features_low",
        nargs=None,
        type=int,
        dest="n_features_low",
        required=False,
        default=consts.DEFAULT_N_FEATURES_LOW,
        help="All cells with a number of features by counts "
        "less than this threshold will be filtered out. "
        "Only used if filtering strategy is `threshold`.",
    )
    subparser.add_argument(
        "--n_features_hi",
        nargs=None,
        type=int,
        dest="n_features_hi",
        required=False,
        default=consts.DEFAULT_N_FEATURES_HI,
        help="All cells with a number of features by counts "
        "greater than this threshold will be filtered out. "
        "Only used if filtering strategy is `threshold`.",
    )
    subparser.add_argument(
        "--pct_counts_mt_hi",
        nargs=None,
        type=float,
        dest="pct_counts_mt_hi",
        required=False,
        default=consts.DEFAULT_PCT_COUNTS_MT_HI,
        help="All cells with a percentage of counts in mitochondrial "
        "features greater than this threshold will be filtered out. "
        "Only used in 'rna' mode regardless of filtering strategy.",
    )
    subparser.add_argument(
        "--total_counts_low",
        nargs=None,
        type=int,
        dest="total_counts_low",
        required=False,
        default=consts.DEFAULT_TOTAL_COUNTS_LOW,
        help="All cells with a total counts "
        "less than this threshold will be filtered out. "
        "Only used in 'atac' mode if filtering strategy is `threshold`.",
    )
    subparser.add_argument(
        "--total_counts_hi",
        nargs=None,
        type=int,
        dest="total_counts_hi",
        required=False,
        default=consts.DEFAULT_TOTAL_COUNTS_HI,
        help="All cells with a total counts "
        "greater than this threshold will be filtered out. "
        "Only used in 'atac' mode if filtering strategy is `threshold`.",
    )
    subparser.add_argument(
        "--ns_hi",
        nargs=None,
        type=float,
        dest="ns_hi",
        required=False,
        default=consts.DEFAULT_NS_HI,
        help="All cells with a nucleosome signal greater than this "
        "threshold will be filtered out. "
        "Only used in 'atac' mode when filtering strategy is `threshold`.",
    )
    subparser.add_argument(
        "--tss_low",
        nargs=None,
        type=float,
        dest="tss_low",
        required=False,
        default=consts.DEFAULT_TSS_LOW,
        help="All cells with a TSS enrichment less than this "
        "threshold will be filtered out. "
        "Only used in 'atac' mode when filtering strategy is `threshold`.",
    )
    subparser.add_argument(
        "--tss_hi",
        nargs=None,
        type=float,
        dest="tss_hi",
        required=False,
        default=consts.DEFAULT_TSS_HI,
        help="All cells with a TSS enrichment greater than this "
        "threshold will be filtered out. "
        "Only used in 'atac' mode when filtering strategy is `threshold`.",
    )
    subparser.add_argument(
        "--atac_qc_tool",
        nargs=None,
        type=str,
        choices=["muon", "pycistopic"],
        dest="atac_qc_tool",
        required=False,
        default="muon",
        help="Tool to use for ATAC-seq QC. "
        "Currently supports either Muon or PyCisTopic. "
        "At a minimum, the tool will be used to calculate the per barcode TSS enrichment"
        "Other metrics may be calculated depending on the tool. "
        "For instance, if pycistopic is selected, it will also output sample level metrics "
        "in the output directory.",
    )
    subparser.add_argument(
        "--n_for_ns_calc",
        nargs=None,
        type=int,
        dest="n_for_ns_calc",
        required=False,
        default=consts.DEFAULT_N_FOR_NS_CALC,
        help="Number of cells to use for calculating the "
        "nucleosome signal. Used when 'atac_qc_tool' is 'muon'.",
    )
    subparser.add_argument(
        "--n_tss",
        nargs=None,
        type=int,
        dest="n_tss",
        required=False,
        default=consts.DEFAULT_N_TSS,
        help="Number of TSS's to use for calculating the "
        "TSS enrichment. Used when 'atac_qc_tool' is 'muon'.",
    )
    subparser.add_argument(
        "--min_cells_per_feature",
        nargs=None,
        type=int,
        dest="min_cells_per_feature",
        required=False,
        default=consts.DEFAULT_MIN_CELLS_PER_FEATURE,
        help="The number of cells a feature must be detected in "
        "to be included in the analysis "
        "Only use this if you would like to remove very rare features "
        "in the early stages of the analysis. ",
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
