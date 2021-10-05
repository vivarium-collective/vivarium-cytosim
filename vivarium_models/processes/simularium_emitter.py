import json
from typing import Any, Dict

from vivarium.core.emitter import Emitter
        
import numpy as np
from simulariumio import (
    TrajectoryConverter,
    TrajectoryData,
    AgentData,
    MetaData,
    UnitData,
)

class SimulariumEmitter(Emitter):

    def __init__(self, config: Dict[str, str]) -> None:
        super().__init__(config)
        self.configuration_data = None
        self.saved_data: Dict[float, Dict[str, Any]] = {}

    def emit(self, data: Dict[str, Any]) -> None:
        """
        Emit the timeseries history portion of ``data``, which is
        ``data['data'] if data['table'] == 'history'`` and put it at
        ``data['data']['time']`` in the history.
        """
        if data['table'] == 'configuration':
            self.configuration_data = data['data']
            assert 'processes' in self.configuration_data, 'please emit processes'
        if data['table'] == 'history':
            emit_data = data['data']
            time = emit_data['time']
            self.saved_data[time] = {
                key: value for key, value in emit_data.items()
                if key not in ['time']}

    def get_simularium_fibers(self, filaments):
        box_size = 150.0 # get from configuration data
        actin_radius = 3.0 # get from configuration data
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
        subpoints = np.array([subpoints], dtype=float)

        return TrajectoryConverter(
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
                    subpoints=subpoints,
                ),
                time_units=UnitData("ns"),  # nanoseconds
                spatial_units=UnitData("nm"),  # nanometers
            )
        ).to_JSON()

    def get_simularium_monomers(self, monomers):
        box_size = 150.0 # get from configuration data
        actin_radius = 3.0 # get from configuration data
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

        return TrajectoryData(
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

    def get_data(self) -> dict:
        """ Return the accumulated timeseries history of "emitted" data. """
        import ipdb; ipdb.set_trace()

        results = []
        for time, state in self.saved_data.items():
            results = self.get_simularium_fibers(state['fibers'])

        # write to file
