
import os
import time
import logging
import argparse
from tqdm.auto import tqdm
import snapatac2 as snap
import scanpy as sc

import matplotlib.pyplot as plt  # This needs to be after matplotlib.use('Agg')

logger = logging.getLogger("cellcommander")

def single_sample_recipe(
    frag_file: str,
    outdir_path: str,
    bin_size: int = 500,
    num_features: int = 50000,
    min_load_tsse: int = 1,
    min_load_num_fragments: int = 500,
    min_num_fragments: int = 1000,
    max_num_fragments: int = None,
    sorted_by_barcode: bool = True,
    chunk_size: int = 2000,
    low_memory: bool = True,
    clustering_resolution: float = 1.0,
    gene_activity: bool = True,
    save_intermediate: bool = False,
):

    # Log snapATAC version
    logger.info(f"Running standard single sample workflow with snapATAC version {snap.__version__}")

    # Load in from fragment file into memory
    logger.info("Loading in data using `import_data` function without file backing")
    adata = snap.pp.import_data(
        fragment_file=frag_file,
        genome=snap.genome.hg38,
        min_tsse=min_load_tsse,
        min_num_fragments=min_load_num_fragments,
        sorted_by_barcode=sorted_by_barcode,
        chunk_size=chunk_size,
        low_memory=low_memory,
    )

    # Plot TSSe distribution vs number of fragments
    logger.info("Plotting TSSe distribution vs number of fragments")
    snap.pl.tsse(adata, interactive=False, out_file=os.path.join(outdir_path, "nfrag_vs_tsse.png"))

    # Filter out low quality cells
    logger.info(f"Filtering out low quality cells with tsse<{min_tsse} and min_num_fragments<{min_num_fragments}, max_num_fragments>{max_num_fragments}")
    snap.pp.filter_cells(adata, min_tsse=min_tsse, min_counts=min_num_fragments, max_counts=max_num_fragments)

    # Add a 5kb tile matrix
    logger.info(f"Adding a {bin_size}bp tile matrix")
    snap.pp.add_tile_matrix(adata, bin_size=bin_size)

    # Select the top accessible features
    logger.info(f"Selecting the top {num_features} accessible features")
    snap.pp.select_features(adata, n_features=num_features)

    # Run scrublet
    logger.info("Running scrublet")
    snap.pp.scrublet(adata)

    # Save the processed data
    if save_intermediate:
        logger.info("Saving the qc data prior to doublet filtering")
        adata.write(os.path.join(outdir_path, f"qc.h5ad"))

    # Filter out doublets
    logger.info("Filtering out doublets")
    snap.pp.filter_doublets(adata)
    
    # Run the spectral embedding
    logger.info("Running spectral embedding")
    snap.tl.spectral(adata)

    # Run UMAP
    logger.info("Running UMAP")
    snap.tl.umap(adata)

    # Find nearest neighbor graph
    logger.info("Finding nearest neighbor graph")
    snap.pp.knn(adata, use_rep="X_spectral")

    # Cluster data
    logger.info(f"Clustering data using Leiden algorithm with resolution {clustering_resolution}")
    snap.tl.leiden(adata, resolution=clustering_resolution, key_added=f"leiden_{clustering_resolution}")

    # Plot the UMAP with clusters
    logger.info("Plotting UMAP with clusters")
    snap.pl.umap(adata, color=f"leiden", interactive=False, out_file=os.path.join(outdir_path, f"umap_leiden_{clustering_resolution}.png"))

    # Save updated data
    logger.info("Saving clustered data")
    adata.write(os.path.join(outdir_path, f"clustered.h5ad"))

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
