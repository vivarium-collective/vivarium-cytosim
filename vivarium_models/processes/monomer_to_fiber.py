import numpy as np

from vivarium.core.process import Deriver
from vivarium.core.engine import Engine, pf
from simularium_models_util.actin import ActinUtil, ActinTestData

from ..util import agents_update


class MonomerToFiber(Deriver):
    defaults = {}

    def __init__(self, parameters=None):
        super().__init__(parameters)

    def ports_schema(self):
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
                    "points": {
                        "_default": [],
                        "_updater": "set",
                        "_emit": True,
                    },
                }
            },
            "monomers": {
                "box_center": {
                    "_default": np.array([3000.0, 1000.0, 1000.0]),
                    "_updater": "set",
                    "_emit": True,
                },
                "box_size": {
                    "_default": 500.0,
                    "_updater": "set",
                    "_emit": True,
                },
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
        print("in monomer to fiber deriver next update")

        monomers = states["monomers"]
        monomer_box_center = monomers["box_center"]
        monomer_box_size = monomers["box_size"]
        previous_fibers = states["fibers"]

        print(f"box_size = {monomer_box_size}")
        monomer_fibers = MonomerToFiber.generate_fibers_from_monomers(
            monomers, monomer_box_center, monomer_box_size
        )

        fiber_update = agents_update(previous_fibers, monomer_fibers)

        return {"fibers": fiber_update}

    @staticmethod
    def get_next_actin_id(prev_actin_id, neighbor_ids):
        """
        Get the ID of the neighbor other than
        the neighbor with the prev_actin_id
        """
        if prev_actin_id < 0:
            return neighbor_ids[0]
        for neighbor_id in neighbor_ids:
            if neighbor_id != prev_actin_id:
                return neighbor_id
        return -1

    @staticmethod
    def get_actin_monomer_positions(
        start_actin_id, particles, box_center, prev_actin_id=-1, result=None
    ):
        """
        Get monomer positions for an actin chain starting
        at the pointed end's start_actin_id
        """
        if result is None:
            result = []
        start_actin = particles[start_actin_id]
        result.append(start_actin["position"] + box_center)
        next_actin_id = MonomerToFiber.get_next_actin_id(
            prev_actin_id, start_actin["neighbor_ids"]
        )
        if next_actin_id < 0:
            return result
        return MonomerToFiber.get_actin_monomer_positions(
            next_actin_id, particles, box_center, start_actin_id, result
        )

    @staticmethod
    def get_fiber_end_point(end_direction, positions, box_size):
        """
        Since we can only calculate axis position
        given a neighbor in each direction,
        use the positions of the first and second actin from the end
        to extrapolate the end axis position
        """
        axis_pos1 = ActinUtil.get_actin_axis_position(
            positions[0:3] if end_direction <= 0 else positions[-3:], box_size
        )
        axis_pos2 = ActinUtil.get_actin_axis_position(
            positions[1:4] if end_direction <= 0 else positions[-4:-1], box_size
        )
        return axis_pos1 - 1.5 * (axis_pos2 - axis_pos1)

    @staticmethod
    def get_fiber(topology_id, pointed_actin_id, monomers, box_center, box_size):
        """
        Get data for a fiber from a chain of particles
        """
        positions = MonomerToFiber.get_actin_monomer_positions(
            pointed_actin_id, monomers["particles"], box_center
        )
        return {
            "type_name": monomers["topologies"][topology_id]["type_name"],
            "points": [
                MonomerToFiber.get_fiber_end_point(-1, positions, box_size),
                MonomerToFiber.get_fiber_end_point(1, positions, box_size),
            ],
        }

    @staticmethod
    def generate_fibers_from_monomers(monomers, box_center, box_size=500.0):
        """
        Transform monomer data into fiber data
        """
        pointed_actin_ids = {}
        for topology_id in monomers["topologies"]:
            # assume first element is the pointed end
            pointed_actin_ids[topology_id] = monomers["topologies"][topology_id][
                "particle_ids"
            ][0]
        result = {}
        for topology_id in pointed_actin_ids:
            result[topology_id] = MonomerToFiber.get_fiber(
                topology_id,
                pointed_actin_ids[topology_id],
                monomers,
                box_center,
                box_size,
            )
        return result


def test_fiber_to_monomer():
    monomer_data = ActinTestData.linear_actin_monomers()
    monomer_to_fiber = MonomerToFiber()

    engine = Engine(
        {
            "processes": {"monomer_to_fiber": monomer_to_fiber},
            "topology": {
                "monomer_to_fiber": {
                    "fibers": ("fibers",),
                    "monomers": ("monomers",),
                }
            },
            "initial_state": {"fibers": {}, "monomers": monomer_data},
        }
    )

    engine.update(1.0)

    output = engine.emitter.get_data()
    print(pf(output))


if __name__ == "__main__":
    test_fiber_to_monomer()
