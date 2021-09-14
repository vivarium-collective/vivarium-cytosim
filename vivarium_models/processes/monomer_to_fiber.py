import numpy as np

from vivarium.core.process import Deriver
from vivarium.core.engine import Engine, pf
from simularium_models_util.actin import ActinUtil, ActinTestData


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


def get_actin_monomer_positions(
    start_actin_id, particles, box_size, prev_actin_id=-1, result=None
):
    """
    Get monomer positions for an actin chain starting
    at the pointed end's start_actin_id
    """
    if result is None:
        result = []
    start_actin = particles[start_actin_id]
    result.append(start_actin["position"])
    next_actin_id = MonomerToFiber.get_next_actin_id(
        prev_actin_id, start_actin["neighbor_ids"]
    )
    if next_actin_id < 0:
        return result
    return MonomerToFiber.get_actin_monomer_positions(
        next_actin_id, particles, box_size, start_actin_id, result
    )


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


def get_fiber(topology_id, pointed_actin_id, monomers, box_size):
    """
    Get data for a fiber from a chain of particles
    """
    positions = MonomerToFiber.get_actin_monomer_positions(
        pointed_actin_id, monomers["particles"], box_size
    )
    return {
        "type_name": monomers["topologies"][topology_id]["type_name"],
        "points": [
            MonomerToFiber.get_fiber_end_point(-1, positions, box_size),
            MonomerToFiber.get_fiber_end_point(1, positions, box_size),
        ],
    }


def generate_fibers_from_monomers(monomers, box_size=150.0):
    """
    Transform monomer data into fiber data
    """
    pointed_actin_ids = {}
    for topology_id in monomers["topologies"]:
        for particle_id in monomers["topologies"][topology_id]["particle_ids"]:
            if "pointed" in monomers["particles"][particle_id]["type_name"]:
                pointed_actin_ids[topology_id] = particle_id
                break
    result = {}
    for topology_id in pointed_actin_ids:
        result[topology_id] = MonomerToFiber.get_fiber(
            topology_id, pointed_actin_ids[topology_id], monomers, box_size
        )
    return result


class MonomerToFiber(Deriver):
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
        monomers = states["monomers"]
        previous_fibers = states["fibers"]

        monomer_fibers = generate_fibers_from_monomers(monomers)

        fiber_update = {}
        if len(previous_fibers) > 0:
            fiber_update["_delete"] = list(previous_fibers.keys())

        add_fibers = []
        for fiber_id in monomer_fibers:
            add_fibers.append(
                {
                    "key": fiber_id,
                    "state": {
                        "type_name": monomer_fibers[fiber_id]["type_name"],
                        "points": monomer_fibers[fiber_id]["points"],
                    },
                }
            )
        fiber_update["_add"] = add_fibers

        return {"fibers": fiber_update}


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
