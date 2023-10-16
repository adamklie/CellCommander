"""Utility functions for qc module."""

from typing import Iterable, Optional, Tuple

import numpy as np
from scipy.stats import median_abs_deviation


# Function to define outliers based on median absolute deviations
def is_outlier(adata, metric: str, nmads: int) -> np.ndarray:
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier
