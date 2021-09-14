import numpy as np

from vivarium.core.process import Deriver
from vivarium.core.engine import Engine, pf

from simularium_models_util.actin import ActinGenerator, ActinTestData, FiberData
from ..util import create_monomer_update


class FiberToMonomer(Deriver):
    defaults = {}

    def __init__(self, parameters=None):
        super().__init__(parameters)

    def ports_schema(self):
        return {
            "fibers": {
                "*": {
                    "type_name": {
                        "_default": "",
                        "_updater": "set",
                        "_emit": True,
                    },
                    "points": {
                        "_default": [],
                        "_updater": "set",
                        "_emit": True,
                    },
                }
            },
            "monomers": {
                "topologies": {
                    "*": {
                        "type_name": {
                            "_default": "",
                            "_updater": "set",
                            "_emit": True,
                        },
                        "particle_ids": {
                            "_default": [],
                            "_updater": "set",
                            "_emit": True,
                        },
                    }
                },
                "particles": {
                    "*": {
                        "type_name": {
                            "_default": "",
                            "_updater": "set",
                            "_emit": True,
                        },
                        "position": {
                            "_default": np.zeros(3),
                            "_updater": "set",
                            "_emit": True,
                        },
                        "neighbor_ids": {
                            "_default": [],
                            "_updater": "set",
                            "_emit": True,
                        },
                    }
                },
            },
        }

    def next_update(self, timestep, states):
        fiber_data = states["fibers"]
        previous_monomers = states["monomers"]

        fibers = [
            FiberData(fiber_id, fiber_data[fiber_id]["points"])
            for fiber_id in fiber_data
        ]
        fiber_monomers = ActinGenerator.get_monomers(fibers, 0)
        for particle_id in fiber_monomers["particles"]:
            fiber_monomers["particles"][particle_id] = dict(
                fiber_monomers["particles"][particle_id]
            )
        return create_monomer_update(previous_monomers, fiber_monomers)


def get_initial_fiber_data():
    fibers = ActinTestData.linear_actin_fiber()
    fibers_dict = {}
    for fiber in fibers:
        fibers_dict[fiber.fiber_id] = dict(fiber)
    return fibers_dict


def test_fiber_to_monomer():
    fiber_data = get_initial_fiber_data()
    fiber_to_monomer = FiberToMonomer()

    engine = Engine(
        {
            "processes": {"fiber_to_monomer": fiber_to_monomer},
            "topology": {
                "fiber_to_monomer": {
                    "fibers": ("fibers",),
                    "monomers": ("monomers",),
                }
            },
            "initial_state": {"fibers": fiber_data, "monomers": {}},
        }
    )

    engine.update(1.0)

    output = engine.emitter.get_data()
    print(pf(output))


if __name__ == "__main__":
    test_fiber_to_monomer()
