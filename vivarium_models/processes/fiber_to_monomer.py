import numpy as np

from vivarium.core.process import Deriver
from simularium_models_util.actin import FiberData, CurvePointData, ActinAnalyzer, ActinGenerator, ActinUtil


class FiberToMonomer(Deriver):
    defaults = {}

    def __init__(self, parameters=None):
        super().__init__(parameters)

    def ports_schema(self):
        return {
            'fibers': {
                },
            'monomers': {
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
            }
        }

    def next_update(self, timestep, states):
        fibers = states['fibers']
        monomers = states['monomers']

        # ActinUtil.add_fibers_from_data(simulation, fibers_data)
        fiber_monomers = ActinGenerator.get_monomers(fibers)

        monomers = {}
        for fiber in fiber_monomers:
            for monomer_type, position, neighbors in enumerate(zip(**fiber)):
                monomers.append({
                    'type': monomer_type,
                    'position': position,
                    'neighbors': neighbors})
            


        import ipdb; ipdb.set_trace()

        return {
            'monomers': fiber_monomers}


def get_initial_fiber_data():
    return [
        FiberData(
            0,
            [
                CurvePointData(
                    np.array([-90, 0, 0]),
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
    update = fiber_to_monomer.next_update(0, {'fibers': fiber_data})

    import ipdb; ipdb.set_trace()


if __name__ == '__main__':
    test_fiber_to_monomer()
