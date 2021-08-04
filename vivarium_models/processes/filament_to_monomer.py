import copy
import numpy as np

from vivarium.core.process import Deriver, Process
from vivarium.core.engine import Engine, pf
from simularium_models_util.actin import ActinGenerator, ActinTestData, FiberData


class FilamentToMonomer(Deriver):
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

    def next_update(self, timestep, states):
        filament_data = states["filaments"]
        previous_monomers = states["monomers"]

        fibers = [FiberData(fiber_id, filament_data[fiber_id]["points"]) for fiber_id in filament_data]
        fiber_monomers = ActinGenerator.get_monomers(fibers, 0)
        for particle_id in fiber_monomers["particles"]:
            fiber_monomers["particles"][particle_id] = dict(fiber_monomers["particles"][particle_id])

        topology_update = {}
        particle_update = {}
        if len(previous_monomers["topologies"]) > 0:
            topology_update["_delete"] = list(previous_monomers["topologies"].keys())
        if len(previous_monomers["particles"]) > 0:
            particle_update["_delete"] = list(previous_monomers["particles"].keys())

        add_topologies = []
        add_particles = []
        for fiber_id in fiber_monomers["topologies"]:
            particle_ids = fiber_monomers["topologies"][fiber_id]["particle_ids"]
            add_topologies.append({
                "key": fiber_id,
                "state": {
                    "type_name": fiber_monomers["topologies"][fiber_id]["type_name"],
                    "particle_ids": particle_ids,
                },
            })
            for monomer_id in particle_ids:
                monomer = fiber_monomers["particles"][monomer_id]
                add_particles.append({
                    "key": monomer_id,
                    "state": {
                        "type_name": monomer["type_name"],
                        "position": monomer["position"],
                        "neighbor_ids": monomer["neighbor_ids"],
                    },
                })
        topology_update["_add"] = add_topologies
        particle_update["_add"] = add_particles

        return {
            "monomers": {
                "topologies": topology_update,
                "particles": particle_update}}


def get_initial_filament_data():
    fibers = ActinTestData.linear_actin_fiber()
    filaments_dict = {}
    for fiber in fibers:
        filaments_dict[fiber.fiber_id] = dict(fiber)
    return filaments_dict
        

def test_filament_to_monomer():
    filament_data = get_initial_filament_data()
    filament_to_monomer = FilamentToMonomer()

    engine = Engine({
        "processes": {
            "filament_to_monomer": filament_to_monomer},
        "topology": {
            "filament_to_monomer": {
                "filaments": ("filaments",),
                "monomers": ("monomers",)}},
        "initial_state": {
            "filaments": filament_data,
            "monomers": {}}})

    engine.update(1.0)

    output = engine.emitter.get_data()
    print(pf(output))


if __name__ == "__main__":
    test_filament_to_monomer()
