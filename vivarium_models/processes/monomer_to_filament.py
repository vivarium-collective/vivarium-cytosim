import numpy as np

from vivarium.core.process import Deriver, Process
from vivarium.core.engine import Engine, pf
from simularium_models_util.actin import ActinUtil, ActinTestData

class MonomerToFilament(Deriver):
    defaults = {}

    def __init__(self, parameters=None):
        super().__init__(parameters)

    def ports_schema(self):
        return {
            "filaments": {
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
                    }}},
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
                        }}},
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
                        }}}}}

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
    def get_actin_monomer_positions(start_actin_id, particles, box_size, prev_actin_id=-1, result=None):
        """
        Get monomer positions for an actin chain starting 
        at the pointed end's start_actin_id
        """
        if result is None:
            result = []
        start_actin = particles[start_actin_id]
        result.append(start_actin["position"])
        next_actin_id = MonomerToFilament.get_next_actin_id(prev_actin_id, start_actin["neighbor_ids"])
        if next_actin_id < 0:
            return result
        return MonomerToFilament.get_actin_monomer_positions(next_actin_id, particles, box_size, start_actin_id, result)

    @staticmethod
    def get_filament_end_point(end_direction, positions, box_size):
        """
        Since we can only calculate axis position 
        given a neighbor in each direction,
        use the positions of the first and second actin from the end 
        to extrapolate the end axis position
        """
        axis_pos1 = ActinUtil.get_actin_axis_position(positions[0:3] if end_direction <= 0 else positions[-3:], box_size)
        axis_pos2 = ActinUtil.get_actin_axis_position(positions[1:4] if end_direction <= 0 else positions[-4:-1], box_size)
        return axis_pos1 - 1.5 * (axis_pos2 - axis_pos1)

    @staticmethod
    def get_filament(topology_id, pointed_actin_id, monomers, box_size):
        positions = MonomerToFilament.get_actin_monomer_positions(pointed_actin_id, monomers["particles"], box_size)
        return {
            "type_name": monomers["topologies"][topology_id]["type_name"],
            "points": [
                MonomerToFilament.get_filament_end_point(-1, positions, box_size),
                MonomerToFilament.get_filament_end_point(1, positions, box_size),
            ]
        }

    @staticmethod
    def generate_filaments_from_monomers(monomers, box_size=150.):
        pointed_actin_ids = {}
        for topology_id in monomers["topologies"]:
            for particle_id in monomers["topologies"][topology_id]["particle_ids"]:
                if "pointed" in monomers["particles"][particle_id]["type_name"]:
                    pointed_actin_ids[topology_id] = particle_id
                    break
        result = {}
        for topology_id in pointed_actin_ids:
            result[topology_id] = MonomerToFilament.get_filament(topology_id, pointed_actin_ids[topology_id], monomers, box_size)
        return result

    def next_update(self, timestep, states):
        monomers = states["monomers"]
        previous_filaments = states["filaments"]

        monomer_filaments = MonomerToFilament.generate_filaments_from_monomers(monomers)

        filament_update = {}
        if len(previous_filaments) > 0:
            filament_update["_delete"] = list(previous_filaments.keys())

        add_filaments = []
        for filament_id in monomer_filaments:
            add_filaments.append({
                "key": filament_id,
                "state": {
                    "type_name": monomer_filaments[filament_id]["type_name"],
                    "points": monomer_filaments[filament_id]["points"],
                },
            })
        filament_update["_add"] = add_filaments

        return {
            "filaments": filament_update}


def test_filament_to_monomer():
    monomer_data = ActinTestData.linear_actin_monomers()
    monomer_to_filament = MonomerToFilament()

    engine = Engine({
        "processes": {
            "monomer_to_filament": monomer_to_filament},
        "topology": {
            "monomer_to_filament": {
                "filaments": ("filaments",),
                "monomers": ("monomers",)}},
        "initial_state": {
            "filaments": {},
            "monomers": monomer_data}})

    engine.update(1.0)

    output = engine.emitter.get_data()
    print(pf(output))


if __name__ == "__main__":
    test_filament_to_monomer()
