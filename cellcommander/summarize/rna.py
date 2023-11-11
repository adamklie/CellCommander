import os
import glob


def write_qc_log_data(log_data, markdown_text):
    markdown_text += " - Filtered initial cells with < 20 genes\n"
    for line in log_data:
        if "Low nfeatures threshold" in line:
            to_add = line.split("INFO - ")[1].strip() + "\n"
            markdown_text += " - " + to_add
        if "High nfeatures threshold" in line:
            to_add = line.split("INFO - ")[1].strip() + "\n"
            markdown_text += " - " + to_add
        if "Percent counts in MT threshold" in line:
            to_add = line.split("INFO - ")[1].strip() + "\n"
            markdown_text += " - " + to_add
        if "Number of cells after filtering of low quality cells:" in line:
            to_add = line.split("Number of cells after filtering of low quality cells:")[1].strip()
            markdown_text += " - Retained " + to_add + " cells\n"
    markdown_text += " - No gene filtering\n\n"
    return markdown_text


def add_qc_images(input_dir, outdir_path, markdown_text):
    # Add image to markdown, make sure it's centered
    if os.path.exists(os.path.join(input_dir, "total_counts_vs_n_genes_after_qc.png")):
        os.system(f"cp {os.path.join(input_dir, 'total_counts_vs_n_genes_after_qc.png')} {outdir_path}")
        image_caption = "Plot of total counts vs number of genes after QC colored by mitochondrial count percentage. These should correspond " \
                        "to the thresholds listed"
        html_text = f"<p align=\"center\"><img src=\"total_counts_vs_n_genes_after_qc.png\" width=\"500\"><br>{image_caption}</p>\n\n"
        markdown_text += html_text
    return markdown_text


def write_remove_background_log_data(log_data, markdown_text):
    # Add markers that were used
    for line in log_data:
        if "Reading in SoupX markers from" in line:
            to_add = line.split("INFO - Reading in SoupX markers from")[1].strip()
            markdown_text += " - SoupX markers: `" + to_add + "`\n\n"
    return markdown_text


def add_remove_background_images(input_dir, outdir_path, markdown_text):
    # Add image to markdown, make sure it's centered
    if os.path.exists(os.path.join(input_dir, "initial_soupx_groups_umap.png")):
        os.system(f"cp {os.path.join(input_dir, 'initial_soupx_groups_umap.png')} {outdir_path}")
        image_caption = "UMAP of data fed into SoupX with initial clustering of cells. " \
            "This clustering is used by SoupX to inform background removal. See the method page for more details. " \
            "The UMAP plot is generated after shifted logarithm normalization, subsetting down to HVGs and running PCA with ScanPy. " \
            "Leiden clustering performed with 30 neighbors, 50 components and a resolution of 0.5.\n\n"
        html_text = f"<p align=\"center\"><img src=\"initial_soupx_groups_umap.png\" width=\"500\"><br>{image_caption}</p>\n\n"
        markdown_text += html_text
    return markdown_text


def write_doublet_detection_log_data(log_data, markdown_text):
    # Add markers that were used
    methods_used = []
    used_cellranger = True
    cells_retained = None
    for line in log_data:
        if "Running" in line and "for doublet detection" in line:
            method = line.split("Running ")[1].split(" for doublet detection")[0]
            methods_used.append(method)
    markdown_text += " - Called doublets with " + ", ".join(methods_used) + "\n"
    if used_cellranger:
        markdown_text += " - Filtered out CellRanger called doublets\n"
    if cells_retained:
        markdown_text += " - " + cells_retained + " retained\n\n"
    else:
        markdown_text += " - No doublet filtering reported\n\n"
    return markdown_text


def add_doublet_detection_images(input_dir, outdir_path, markdown_text):
    # Add image to markdown, make sure it's centered
    if os.path.exists(os.path.join(input_dir, "doublet_leiden_scores_umap.png")):
        os.system(f"cp {os.path.join(input_dir, 'doublet_leiden_scores_umap.png')} {outdir_path}")
        image_caption = "UMAPs with doublet scores from all methods used overlayed. These UMAP plots is generated after shifted logarithm normalization, " \
            "subsetting down to HVGs and running PCA with ScanPy. Leiden clustering performed with 30 neighbors, 50 components and a resolution of 0.5.\n\n"
        html_text = f"<p align=\"center\"><img src=\"doublet_leiden_scores_umap.png\" width=\"1500\"><br>{image_caption}</p>\n\n"
        markdown_text += html_text
    return markdown_text


def write_normalization_log_data(log_data, markdown_text):
    # Add markers that were used
    methods_used = []
    genes_filtered = False
    for line in log_data:
        if "Running" in line and "normalization" in line:
            method = line.split("Running ")[1].split(" normalization")[0]
            methods_used.append(method)
        if "Number of genes after cell filter:" in line:
            genes_filtered = True
            to_add = line.split("Number of genes after cell filter:")[1].strip()
    markdown_text += " - Normalized with " + ", ".join(methods_used) + "\n"
    if genes_filtered:
        markdown_text += " - Retained " + to_add + " genes\n\n"
    else:
        markdown_text += " - No gene filtering\n\n"
    return markdown_text


def add_normalization_images(input_dir, outdir_path, markdown_text):
    # Add all plots that end with _distribution.png in a side by side grid, make sure it's centered
    normalization_files = glob.glob(os.path.join(input_dir, "*"))
    normalization_files = [x for x in normalization_files if x.endswith("_distribution.png")]
    normalization_files = [x for x in normalization_files if "unnormalized_total_counts_distribution.png" not in x]
    normalization_files = [x for x in normalization_files if "depth_normalization_distribution.png" not in x]
    html_text = ""
    for file in normalization_files:
        os.system(f"cp {file} {outdir_path}")
        html_text += f"<img src=\"{os.path.basename(file)}\" width=\"600\">"
    html_text += "\n\n"
    markdown_text += html_text
    return markdown_text


def write_feature_selection_log_data(log_data, markdown_text):
    # Add markers that were used
    methods_used = []
    for line in log_data:
        if "Running" in line and "feature selection" in line:
            method = line.split("Running ")[1].split(" feature selection")[0]
            methods_used.append(method)
    markdown_text += " - Selected features with " + ", ".join(methods_used) + "method\n\n"
    return markdown_text


def add_feature_selection_images(input_dir, outdir_path, markdown_text):
    # Add all plots that end with .png in a side by side grid, make sure it's centered
    feature_selection_files = glob.glob(os.path.join(input_dir, "*"))
    feature_selection_files = [x for x in feature_selection_files if x.endswith(".png")]
    feature_selection_files = [x for x in feature_selection_files if "depth_normalization" not in x]
    html_text = ""
    for file in feature_selection_files:
        os.system(f"cp {file} {outdir_path}")
        html_text += f"<img src=\"{os.path.basename(file)}\" width=\"750\">"
    html_text += "\n\n"
    markdown_text += html_text
    return markdown_text


def write_dimensionality_reduction_log_data(log_data, markdown_text):
    # Add markers that were used
    markdown_text += " - Ran PCA on SCTransform normalized data (no feature selection needed), then UMAP on top of this. The PCA and UMAP used here were run in Seurat.\n"
    markdown_text += " - Performed an initial clustering with Leiden with 30 neighbors, 50 components (PCs) and resolution of 0.5. Implementation from Scanpy.\n\n"
    return markdown_text


def add_dimensionality_reduction_images(input_dir, outdir_path, markdown_text):
    # Add all plots that end with pdf, make sure it's centered
    dimensionality_reduction_files = glob.glob(os.path.join(input_dir, "*"))
    dimensionality_reduction_files = [x for x in dimensionality_reduction_files if x.endswith("pdf")]
    html_text = ""
    for file in dimensionality_reduction_files:
        os.system(f"cp {file} {outdir_path}")
        image_caption = "UMAP of SCTransform normalized data.  Leiden clustering performed with 30 neighbors, 50 components and a resolution of 0.5."
        html_text += f"<p align=\"center\"><img src=\"{os.path.basename(file)}\" width=\"500\"><br>{image_caption}</p>\n\n"
    html_text += "\n\n"
    markdown_text += html_text
    return markdown_text


def write_annotate_log_data(markdown_text):
    # Add markers that were used
    markdown_text += " - Use the same UMAP but plot some new clusters on it at higher resolution\n"
    markdown_text += "- In conjunction with the accompanying dot plot of known marker genes: annotation/18Oct23/SC.islet.marker_genes.csv\n\n"
    return markdown_text


def add_annotate_images(input_dir, outdir_path, markdown_text):
    # If annotate_clustering_umap.png in annotate directory, copy it over
    html_text = ""
    if os.path.exists(os.path.join(input_dir, "annotate_clustering_umap.png")):
        os.system(f"cp {os.path.join(input_dir, 'annotate_clustering_umap.png')} {outdir_path}")
        image_caption = "UMAP of SCTransform normalized data.  Leiden clustering performed with 30 neighbors, 50 components and a resolution of 0.5."
        html_text += f"<p align=\"center\"><img src=\"annotate_clustering_umap.png\" width=\"500\"><br>{image_caption}</p>\n\n"
    # if marker_gene_dotplot.png in annotate directory, copy it over
    if os.path.exists(os.path.join(input_dir, "marker_gene_dotplot.png")):
        os.system(f"cp {os.path.join(input_dir, 'marker_gene_dotplot.png')} {outdir_path}")
        image_caption = "Dot plot of known sc islet marker genes"
        html_text += f"<p align=\"center\"><img src=\"marker_gene_dotplot.png\" width=\"500\"><br>{image_caption}</p>\n\n"
    markdown_text += html_text

    # if manual_cellid_annotation_umap.png in annotate directory, copy it over, center it
    html_text = ""
    if os.path.exists(os.path.join(input_dir, "manual_cellid_annotation_umap.png")):
        os.system(f"cp {os.path.join(input_dir, 'manual_cellid_annotation_umap.png')} {outdir_path}")
        image_caption = "User input annotations (left) for each of the clusters defined on the right."
        html_text += f"<p align=\"center\"><img src=\"manual_cellid_annotation_umap.png\" width=\"1500\"><br>{image_caption}</p>\n\n"
    markdown_text += html_text
    return markdown_text
