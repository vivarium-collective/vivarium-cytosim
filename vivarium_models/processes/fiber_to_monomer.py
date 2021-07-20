import numpy as np
import uuid

from vivarium.core.process import Deriver, Process
from vivarium.core.engine import Engine
from simularium_models_util.actin import FiberData, CurvePointData, ActinAnalyzer, ActinGenerator, ActinUtil


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

        # ActinUtil.add_fibers_from_data(simulation, fibers_data)
        fiber_monomers = ActinGenerator.get_monomers(fibers)
        # TODO: handle the ends of the fiber, don't make neighbors off the end of the fiber

        monomer_update = {}
        if len(previous_monomers) > 0:
            monomer_update['_delete'] = list(previous_monomers.keys())

        add_monomers = []
        for fiber in fiber_monomers:
            monomers = list(zip(*fiber))
            monomer_ids = [str(uuid.uuid4()) for _ in range(len(monomers) + 1)]

            for index, monomer in enumerate(monomers):
                monomer_type, position, neighbors = monomer
                id = monomer_ids[index]
                add_monomers.append({
                    'key': id,
                    'state': {
                        'type': monomer_type,
                        'position': position,
                        'neighbors': [monomer_ids[neighbor] for neighbor in neighbors]}})

        monomer_update['_add'] = add_monomers

        return {
            'monomers': monomer_update}


def get_initial_fiber_data():
    return [
        FiberData(
            0,
            [
                CurvePointData(
                    np.array([-75, 0, 0]),
                    np.array([1, 0, 0]),
                    0,
                    ),
                    CurvePointData(
                        np.array([10, 0, 0]),
                        np.array([1, 0, 0]),
                        50,
                        ),
            ],
            )
    ]

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
