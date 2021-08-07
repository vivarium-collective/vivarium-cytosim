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


class VisualizeFilament(Deriver):
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
                    },
                }
            },
            "simularium_json": {
                "_default": "",
                "_updater": "set",
                "_emit": True,
            },
        }

    def next_update(self, timestep, states):
        filaments = states["filaments"]

        n_agents = 0
        unique_ids = []
        type_names = []
        n_subpoints = []
        subpoints = []
        for filament_id in filaments:
            filament = filaments[filament_id]
            unique_ids.append(filament_id)
            type_names.append(filament["type_name"])
            n_subpoints.append(len(filament["points"]))
            subpoints.append(filament["points"])
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
                    viz_types=1001.0 * np.ones((1, n_agents)),
                    unique_ids=np.array([unique_ids]),
                    types=[type_names],
                    positions=np.zeros((1, n_agents, 3)),
                    radii=actin_radius * np.ones((1, n_agents)),
                    n_subpoints=np.array([n_subpoints]),
                    subpoints=np.array([subpoints]),
                ),
                time_units=UnitData("ns"),  # nanoseconds
                spatial_units=UnitData("nm"),  # nanometers
            )
        ).to_JSON()

        return {"simularium_json": simularium_json}


def get_initial_filament_data():
    fibers = ActinTestData.linear_actin_fiber()
    filaments_dict = {}
    for fiber in fibers:
        filaments_dict[fiber.fiber_id] = dict(fiber)
    return filaments_dict


def test_visualize_filament():
    filament_data = get_initial_filament_data()
    visualize_filament = VisualizeFilament()

    engine = Engine(
        {
            "processes": {"visualize_filament": visualize_filament},
            "topology": {
                "visualize_filament": {
                    "filaments": ("filaments",),
                    "simularium_json": ("simularium_json",),
                }
            },
            "initial_state": {"filaments": filament_data, "simularium_json": ""},
        }
    )

    engine.update(1.0)

    output = engine.emitter.get_data()
    print(pf(output))
    simularium_json = output[0.0]["simularium_json"]
    with open("filaments.simularium", "w+") as outfile:
        outfile.write(simularium_json)
    print("Saved Simularium output to filaments.simularium")


if __name__ == "__main__":
    test_visualize_filament()
