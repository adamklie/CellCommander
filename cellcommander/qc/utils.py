"""Utility functions for qc module."""

from typing import Iterable, Optional, Tuple

import numpy as np
from anndata import AnnData
from scipy.stats import median_abs_deviation


# Function to define outliers based on median absolute deviations
def is_outlier(adata, metric: str, nmads: int) -> np.ndarray:
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


# Function to calculate median sample level metrics
def sample_level_metrics(adata: AnnData, metrics: Iterable[str]):
    # Sample level metrics
    sample_metrics = {}
    for metric in metrics:
        if metric in adata.obs.columns:
            sample_metrics[f"median_{metric}"] = adata.obs[metric].median()
    return sample_metrics
