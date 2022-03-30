import numpy as np


def fibers_schema():
    return {
        "fibers_box_extent": {
            "_default": np.array([4000.0, 2000.0, 2000.0]),
            "_updater": "set",
            "_emit": True,
        },
        "fibers": {
            "*": {
                "type_name": {
                    "_default": "",
                    "_updater": "set",
                    "_emit": True,
                },
                "points": {  # list of shape (3) numpy arrays
                    "_default": [],
                    "_updater": "set",
                    "_emit": True,
                },
            }
        },
    }
