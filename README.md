# CellCommander
Common single cell analysis tasks made runnable with the command line. This tool is still under active development, so the `dev` branch is the best place to get the latest updates/working code.

# Why?
Single cell analysis is often characterized by “kitchensink” modeling (e.g. throwing everything at your data and seeing what sticks). Many of the most popular data analysis steps and methods are captured in these online books:

 - [Single-cell best practices — Single-cell best practices](https://www.sc-best-practices.org/preamble.html)
 - [Orchestrating Single-Cell Analysis with Bioconductor](https://bioconductor.org/books/3.12/OSCA/)

Though frameworks like scverse, Bioconductor, Seurat and ScanPy have made so many of these methods easily accessible, these methods are still difficult to manipulate into pipelines that you can run for automated analysis of many samples.

**My goal was to build a tool that automates testing out pipelines end-to-end on a variety of data types.**

Below we enumerate the modules that are currently available in CellCommander that one can think of as the building blocks of a single cell analysis pipeline. Each module contains the following in common:
1. An input argument, usually a path to a cell x feature matrix or fragment file
2. An output directory argument, where the output of the module will be saved
3. Generates a log file that contains a history of the commands that were run

# `qc`
⭐ This is a command line module for performing initial QC and filtering of your cell x feature matrix. Use `cellcommander qc --help` to see all of the available options.

## `--modes ["rna", "atac"]`
Determines what metrics should be calculated for the data.
 - [] So far, I’m liking SnapATAC2 and ArchR’s way of calculating QCs — Muon may not be the best option here. `--atac_qc_tool**` where we theoretically could use either of these two to load in the files for initial QC
 
## `--filtering_strategy ["mad", "threshold"]`

- Determines how cell outliers are calculated
- So far, it seems like MAD usually keeps a lot of low quality cells

# `remove-background`

- Change initial clust
- Must be a single modality (RNA) adata
- Currently hardcoded to 4 markers and 4 groups
- Adds 2 obs columns (pre and post soupx clustering)
- Adds soupx counts and raw counts to layers
- Moves soupx counts to .X
- Plots pre and post removal UMAPs and clusterings
- Saves the soupx object as an RDS
- Saves statistics derived from the SoupX run in a pickle file

# `detect-doublets`

<aside>
⭐ This is a command line module for identifying and filtering out doublets from your cell x feature matrix.

</aside>

## `—-method ["scrublet", "scDblFinder", "amulet", "cellranger", "consensus"]`

- Determines what algorithm is used for doublet detection

## `normalize`

<aside>
⭐ This is a command line module for normalizing your cell x feature matrix.

</aside>

## `—-method ["log1p", "scran", "pearson", "depth", "sctransform", "tfidf"]`

- Determines how the matrix is normalized

# `select-features`

<aside>
⭐ This is a command line module for selecting columns of interest from an input cell x feature matrix.

</aside>

# `reduce-dimensions`

<aside>
⭐ This is a command line module for creating a dimensionality reduced version of an input cell x feature matrix.

</aside>

- Will always compute a UMAP and some initial clustering on this dimensionality reduction

## `—-method ["scanpy_default", "seurat_default", "seurat_sctransform", "signac_default", "lsi"]`

- Determines what dimensionality reduction method to run

# `annotate`

<aside>
⭐ This is a command line module for annotating the cells with cell identity and state either manually or using prior knowledge.

</aside>

# `integrate`

<aside>
⭐ This is a command line module for *horizontally* integrating multiple input cell x feature matrices.

</aside>

# `joint-integrate`

<aside>
⭐ This is a command line module for **********vertically********** integrating multiple input cell x feature matrices.

</aside>

# `recipes`

<aside>
⭐ This is a command line module for running pipelines.

</aside>

## `snapatac2`

### `multi-sample`

## `archr`

### `multi-sample`

`**initialize**`

1. Make `arrow` files from list of fragment file or bam paths
2. Make an `ArchR` project object from all of these `arrow` files
3. Write out metadata
4. Add doublet scores and filter out
5. Filter out doublets
    1. Based on ArchR
    2. Based on any passed in doublet list (path to text file)
6. Write out metadata
7. Save the ArchR project

`**preprocess**`

1. Read in `ArchR` project
2. Read in any metadata and inner join on index?
3. Dimensionality reduction with iterative LSI
4. Clustering with Seurat method
5. Optional Harmony correction
6. UMAP

`**visualize**`

## A note on interoperability
https://scverse-tutorials.readthedocs.io/en/latest/notebooks/scverse_data_interoperability.html
https://scverse-tutorials.readthedocs.io/en/latest/index.html
[10/29/2023 Brainstorm — Ideas](https://www.notion.so/10-29-2023-Brainstorm-Ideas-65e4f4f3e5374510829559e013f3d70d?pvs=21)
