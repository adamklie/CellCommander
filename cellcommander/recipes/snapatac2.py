
import os
import logging
from tqdm.auto import tqdm
import pandas as pd
import snapatac2 as snap
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np


logger = logging.getLogger("cellcommander")

def single_sample_recipe(
    frag_file: str,
    outdir_path: str,
    sample_name: str = None,
    bin_size: int = 500,
    num_features: int = 50000,
    min_load_num_fragments: int = 500,
    min_tsse: int = 4,
    min_num_fragments: int = 1000,
    max_num_fragments: int = None,
    sorted_by_barcode: bool = True,
    chunk_size: int = 2000,
    clustering_resolution: float = 1.0,
    gene_activity: bool = True,
    metadata: pd.DataFrame = None,
    additional_doublets: str = None,
    save_intermediate: bool = False,
):

    # Log snapATAC version
    logger.info(f"Running standard single sample workflow with snapATAC version {snap.__version__}")

    # Load in from fragment file into memory
    logger.info("Loading in data using `import_data` function without file backing")
    adata = snap.pp.import_data(
        fragment_file=frag_file,
        chrom_sizes=snap.genome.hg38,
        min_num_fragments=min_load_num_fragments,
        sorted_by_barcode=sorted_by_barcode,
        chunk_size=chunk_size,
        n_jobs=-1,
    )

    # Add sample name to barcode with # in between
    if sample_name is not None:
        logger.info(f"Adding sample name {sample_name} to barcode")
        logger.info(f"Before: {adata.obs.index[0]}")
        adata.obs.index = sample_name + "#" + adata.obs.index
        logger.info(f"After: {adata.obs.index[0]}")
        logger.info("If passing in metadata, make sure to add sample name to barcode column as well")

    # Plot fragment size distribution
    logger.info("Plotting fragment size distribution")
    snap.pl.frag_size_distr(adata, interactive=False, out_file=os.path.join(outdir_path, "frag_size_distr.png"))

    # Plot TSSe distribution vs number of fragments
    logger.info("Plotting TSSe distribution vs number of fragments")
    snap.metrics.tsse(adata, snap.genome.hg38)
    snap.pl.tsse(adata, interactive=False, out_file=os.path.join(outdir_path, "nfrag_vs_tsse.png"))

    # Save the processed data
    if save_intermediate:
        logger.info("Saving the qc data prior to filtering")
        adata.write(os.path.join(outdir_path, f"qc.h5ad"))

    # Filter out low quality cells
    logger.info(f"Filtering out low quality cells with tsse<{min_tsse} and min_num_fragments<{min_num_fragments}, max_num_fragments>{max_num_fragments}")
    adata.obs["log_n_fragment"] = np.log10(adata.obs["n_fragment"] + 1)
    snap.pp.filter_cells(adata, min_tsse=min_tsse, min_counts=min_num_fragments, max_counts=max_num_fragments)

    # Report number of cells after filtering
    logger.info(f"Number of cells after filtering: {adata.shape[0]}")

    # Add a 5kb tile matrix
    logger.info(f"Adding a {bin_size}bp tile matrix")
    snap.pp.add_tile_matrix(adata, bin_size=bin_size)

    # Select the top accessible features
    logger.info(f"Selecting the top {num_features} accessible features")
    snap.pp.select_features(adata, n_features=num_features)

    # Run scrublet
    logger.info("Running scrublet")
    snap.pp.scrublet(adata)

    # Filter out doublets
    logger.info("Filtering out doublets")
    snap.pp.filter_doublets(adata)

    # Filter out additional doublets if passed in
    if additional_doublets is not None:
        logger.info(f"Filtering out additional doublets from {additional_doublets}")
        additional_doublets = pd.read_csv(additional_doublets, header=None, index_col=0).index
        adata = adata[~adata.obs.index.isin(additional_doublets)]
    
    # Report number of cells after filtering
    logger.info(f"Number of cells after filtering doublets: {adata.shape[0]}")
    
    # Add in metadata if passed in
    logger.info("Subsetting data to metadata if passed in")
    if metadata is not None:
        num_intersecting_cells = len(set(metadata.index).intersection(set(adata.obs.index)))
        logger.info(f"Number of cells found in metadata: {num_intersecting_cells}")

        # Subset the object to only include cells in the metadata
        adata = adata[adata.obs.index.isin(metadata.index)]
        
        # Add in the metadata
        adata.obs = adata.obs.merge(metadata, left_index=True, right_index=True, suffixes=("", "_rna"))

        # Check
        logger.info(f"Number of cells after subsetting to metadata: {adata.shape[0]}")

    # Run the spectral embedding
    logger.info("Running spectral embedding")
    snap.tl.spectral(adata)

    # Plot first spectral embedding against log_n_fragment    
    with plt.rc_context({"figure.figsize": (5, 5)}):
        sc.pl.embedding(basis="X_spectral", adata=adata, color="log_n_fragment", show=False)
        plt.savefig(os.path.join(outdir_path, "spectral_embedding.png"))
        plt.close()

    # Run UMAP
    logger.info("Running UMAP")
    snap.tl.umap(adata, use_dims=list(range(1, adata.obsm["X_spectral"].shape[1])))

    # Find nearest neighbor graph
    logger.info("Finding nearest neighbor graph")
    snap.pp.knn(adata, use_rep="X_spectral", use_dims=list(range(1, adata.obsm["X_spectral"].shape[1])))

    # Cluster data
    logger.info(f"Clustering data using Leiden algorithm with resolution {clustering_resolution}")
    snap.tl.leiden(adata, resolution=clustering_resolution, key_added=f"leiden_{clustering_resolution}")

    # Plot the UMAP with clusters
    logger.info("Plotting UMAP with clusters")
    with plt.rc_context({"figure.figsize": (5, 5)}):
        sc.pl.umap(adata, color=["log_n_fragment", "tsse", f"leiden_{clustering_resolution}"], show=False)
        plt.savefig(os.path.join(outdir_path, "umap.png"))
        plt.close()

    # Save updated data
    logger.info("Saving clustered data")
    adata.write(os.path.join(outdir_path, f"clustered.h5ad"))
    adata.obs.to_csv(os.path.join(outdir_path, f"cell_metadata.tsv"), sep="\t")

    # Create a gene matrix
    if gene_activity:
        # Creating gene matrix
        logger.info("Creating gene activity matrix")
        gene_matrix = snap.pp.make_gene_matrix(adata=adata, gene_anno=snap.genome.hg38)

        # Clean up the gene matrix
        logger.info("Filtering and normalizing the gene activity matrix")
        sc.pp.filter_genes(gene_matrix, min_cells=3)
        sc.pp.normalize_total(gene_matrix)
        sc.pp.log1p(gene_matrix)

        # Run MAGIC
        logger.info("Running MAGIC on the gene activity matrix for imputation")
        sc.external.pp.magic(gene_matrix, solver="approximate")

        # Transfer the UMAP from the original data to the gene matrix
        gene_matrix.obsm["X_umap"] = adata.obsm["X_umap"]

        # Save the gene matrix
        logger.info("Saving gene activity matrix")
        gene_matrix.write(os.path.join(outdir_path, f"gene_matrix.h5ad"))
