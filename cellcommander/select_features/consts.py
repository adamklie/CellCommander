"""Constant numbers used in remove-background."""

import numpy as np

# Seed for random number generators.
RANDOM_STATE = 1234

# Variable feature selection constants
DEFAULT_MIN_MEAN = 0.125
DEFAULT_MAX_MEAN = 3
DEFAULT_MIN_DISP = 0.5
DEFAULT_MAX_DISP = np.inf
DEFAULT_N_BINS = 20
DEFAULT_N_TOP_GENES = 3000
DEFAULT_SPAN = 0.3

# Constants across all modes
DEFAULT_INITIAL_CLUST_N_NEIGHBORS = 30
DEFAULT_INITIAL_CLUST_RESOLUTION = 0.5
