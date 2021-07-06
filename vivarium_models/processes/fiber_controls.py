import numpy as np

from vivarium.core.process import Deriver
from simularium_models_util.actin.actin_analyzer import ActinAnalyzer

class FiberControls(Deriver):
    defaults = {}

    def __init__(self, parameters):
        super().__init__(parameters)

    def ports_schema(self):
        return {
            "topologies": {
                "*": {
                    "type": {"_default": ""},
                    "particles": {
                        "*": {
                            "type": {
                                "_default": "",
                            },
                            "position": {
                                "_default": np.zeros(3),
                            },
                            "neighbors": {
                                "_default": [],
                            },
                        },
                    },
                }
            },
            "fibers": {
            }
        }

    def next_update(self, timestep, states):
        import ipdb; ipdb.set_trace()

        

        for topology_key, topology in states['topologies'].items():
            pass
