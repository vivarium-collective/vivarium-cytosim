import numpy as np
import uuid

from vivarium.core.process import Deriver, Process
from vivarium.core.engine import Engine
from simularium_models_util.actin import ActinGenerator, ActinTestData


class FiberToMonomer(Deriver):
    defaults = {}

    def __init__(self, parameters=None):
        super().__init__(parameters)

    def ports_schema(self):
        # add a level here for topologies to match ReaddyActinProcess
        return {
            'fibers': {
                '_default': []},
            'monomers': {
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
            }
        }

    def next_update(self, timestep, states):
        fibers = states['fibers']
        previous_monomers = states['monomers']

        fiber_monomers = ActinGenerator.get_monomers(fibers)

        monomer_update = {}
        if len(previous_monomers) > 0:
            monomer_update['_delete'] = list(previous_monomers.keys())

        add_monomers = []
        for fiber_id in fiber_monomers["topologies"]:
            monomer_ids = fiber_monomers["topologies"][fiber_id]["particle_ids"]
            for monomer_id in monomer_ids:
                monomer = fiber_monomers["particles"][monomer_id]
                add_monomers.append({
                    'key': monomer_id,
                    'state': {
                        "type": monomer["type_name"],
                        "position": monomer["position"],
                        "neighbors": monomer["neighbor_ids"],
                    }
                })

        monomer_update['_add'] = add_monomers

        return {
            'monomers': monomer_update}


def get_initial_fiber_data():
    return ActinTestData.linear_actin_fiber()

def test_fiber_to_monomer():
    fiber_data = get_initial_fiber_data()
    fiber_to_monomer = FiberToMonomer()

    engine = Engine({
        'processes': {
            'fiber_to_monomer': fiber_to_monomer},
        'topology': {
            'fiber_to_monomer': {
                'fibers': ('fibers',),
                'monomers': ('monomers',)}},
        'initial_state': {
            'fibers': fiber_data,
            'monomers': {}}})

    engine.update(1.0)


    # update = fiber_to_monomer.next_update(0, )



if __name__ == '__main__':
    test_fiber_to_monomer()
