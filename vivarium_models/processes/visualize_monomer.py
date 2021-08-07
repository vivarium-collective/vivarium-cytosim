import numpy as np

from vivarium.core.process import Deriver
from vivarium.core.engine import Engine, pf

from simularium_models_util.actin import ActinTestData
from simulariumio import (
    TrajectoryConverter,
    TrajectoryData,
    AgentData,
    MetaData,
    UnitData,
)


class VisualizeMonomer(Deriver):
    defaults = {}

    def __init__(self, parameters=None):
        super().__init__(parameters)

    def ports_schema(self):
        return {
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
            "simularium_json": {
                "_default": "",
                "_updater": "set",
                "_emit": True,
            },
        }

    def next_update(self, timestep, states):
        monomers = states["monomers"]

        n_agents = 0
        unique_ids = []
        type_names = []
        positions = []
        for particle_id in monomers["particles"]:
            particle = monomers["particles"][particle_id]
            unique_ids.append(particle_id)
            type_names.append(particle["type_name"])
            positions.append(particle["position"])
            n_agents += 1

        box_size = 150.0
        actin_radius = 3.0
        simularium_json = TrajectoryConverter(
            TrajectoryData(
                meta_data=MetaData(
                    box_size=np.array([box_size, box_size, box_size]),
                ),
                agent_data=AgentData(
                    times=np.array([0]),
                    n_agents=np.array([n_agents]),
                    viz_types=1000.0 * np.ones((1, n_agents)),
                    unique_ids=np.array([unique_ids]),
                    types=[type_names],
                    positions=np.array([positions]),
                    radii=actin_radius * np.ones((1, n_agents)),
                ),
                time_units=UnitData("ns"),  # nanoseconds
                spatial_units=UnitData("nm"),  # nanometers
            )
        ).to_JSON()

        return {"simularium_json": simularium_json}


def test_visualize_monomer():
    monomer_data = ActinTestData.linear_actin_monomers()
    visualize_monomer = VisualizeMonomer()

    engine = Engine(
        {
            "processes": {"visualize_monomer": visualize_monomer},
            "topology": {
                "visualize_monomer": {
                    "monomers": ("monomers",),
                    "simularium_json": ("simularium_json",),
                }
            },
            "initial_state": {"monomers": monomer_data, "simularium_json": ""},
        }
    )

    engine.update(1.0)

    output = engine.emitter.get_data()
    print(pf(output))
    simularium_json = output[0.0]["simularium_json"]
    with open("monomers.simularium", "w+") as outfile:
        outfile.write(simularium_json)
    print("Saved Simularium output to monomers.simularium")


if __name__ == "__main__":
    test_visualize_monomer()
