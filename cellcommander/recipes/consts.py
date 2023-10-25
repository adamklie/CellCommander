"""Constant numbers used in recipes."""

import numpy as np

RANDOM_STATE = 1234

# SnapATAC2 2.3.1 defaults
DEFAULT_PARAMS = {
    "snapatac2": {
        "single-sample": {
            "analysis": {
                "clustering_resolution": 1
            },
            "feature_selection": {
                "bin_size": 500,
                "num_features": 50000
            },
            "io": {
                "chunk_size": 2000,
                "low_memory": True,
                "min_load_num_fragments": 200,
                "min_load_tsse": 1,
                "sorted_by_barcode": True,
                "save_intermediate": False,
                "gene_activity": True
            },
            "qc": {
                "max_num_fragments": None,
                "min_num_fragments": 1000,
                "min_tsse": 5
            }
        }
    }
}