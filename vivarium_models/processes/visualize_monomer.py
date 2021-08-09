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

        box_size = 150.0
        actin_radius = 3.0
        unique_ids = []
        type_names = []
        positions = []
        edge_ids = []
        edge_positions = []
        for particle_id in monomers["particles"]:
            particle = monomers["particles"][particle_id]
            unique_ids.append(particle_id)
            type_names.append(particle["type_name"])
            positions.append(particle["position"])
            for neighbor_id in particle["neighbor_ids"]:
                edge = (particle_id, neighbor_id)
                reverse_edge = (
                    neighbor_id,
                    particle_id,
                )
                if edge not in edge_ids and reverse_edge not in edge_ids:
                    edge_ids.append(edge)
                    edge_positions.append(
                        [
                            particle["position"],
                            monomers["particles"][neighbor_id]["position"],
                        ]
                    )
        # add fiber agents for edges
        n_agents = len(unique_ids)
        n_edges = len(edge_ids)
        viz_types = 1000.0 * np.ones((1, n_agents + n_edges))
        viz_types[:, n_agents:] += 1
        unique_ids = np.array([unique_ids + [1000 + i for i, e in enumerate(edge_ids)]])
        type_names += ["edge" for edge in edge_ids]
        positions += n_edges * [np.zeros(3)]
        radii = actin_radius * np.ones((1, n_agents + n_edges))
        radii[:, n_agents:] = 1.0
        n_subpoints = np.zeros(n_agents + n_edges)
        n_subpoints[n_agents:] += 2
        subpoints = np.zeros((n_agents + n_edges, 2, 3))
        subpoints[n_agents:] = np.array(edge_positions)

        simularium_json = TrajectoryConverter(
            TrajectoryData(
                meta_data=MetaData(
                    box_size=np.array([box_size, box_size, box_size]),
                ),
                agent_data=AgentData(
                    times=np.array([0]),
                    n_agents=np.array([n_agents + n_edges]),
                    viz_types=viz_types,
                    unique_ids=unique_ids,
                    types=[type_names],
                    positions=np.array([positions]),
                    radii=radii,
                    n_subpoints=np.array([n_subpoints]),
                    subpoints=np.array([subpoints]),
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
