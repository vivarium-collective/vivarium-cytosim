import os
import numpy as np

from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process,
    PROCESS_OUT_DIR,
)
from vivarium.plots.simulation_output import plot_simulation_output

from tqdm import tqdm
from simularium_models_util.actin import ActinSimulation


NAME = "ReaDDy_actin"


class ReaddyActinProcess(Process):
    """
    This process uses ReaDDy to model the dynamics
    of branched actin networks made of spatially explicit monomers
    """

    name = NAME

    defaults = {
        "name": "actin",
        "total_steps": 1e3,
        "timestep": 0.1,
        "box_size": 150.0,  # nm
        "temperature_C": 22.0,  # from Pollard experiments
        "viscosity": 8.1,  # cP, viscosity in cytoplasm
        "force_constant": 250.0,
        "reaction_distance": 1.0,  # nm
        "n_cpu": 4,
        "actin_concentration": 200.0,  # uM
        "arp23_concentration": 10.0,  # uM
        "cap_concentration": 0.0,  # uM
        "n_fibers": 0,
        "fiber_length": 0.0,
        "actin_radius": 2.0,  # nm
        "arp23_radius": 2.0,  # nm
        "cap_radius": 3.0,  # nm
        "dimerize_rate": 2.1e-2,  # 1/ns
        "dimerize_reverse_rate": 1.4e-1,  # 1/ns
        "trimerize_rate": 2.1e-2,  # 1/ns
        "trimerize_reverse_rate": 1.4e-1,  # 1/ns
        "pointed_growth_ATP_rate": 2.4e5,  # 1/ns
        "pointed_growth_ADP_rate": 3.0e4,  # 1/ns
        "pointed_shrink_ATP_rate": 8.0e-15,  # 1/ns
        "pointed_shrink_ADP_rate": 3.0e-15,  # 1/ns
        "barbed_growth_ATP_rate": 2.1e6,  # 1/ns
        "barbed_growth_ADP_rate": 7.0e5,  # 1/ns
        "nucleate_ATP_rate": 2.1e6,  # 1/ns
        "nucleate_ADP_rate": 7.0e5,  # 1/ns
        "barbed_shrink_ATP_rate": 1.4e-14,  # 1/ns
        "barbed_shrink_ADP_rate": 8.0e-14,  # 1/ns
        "arp_bind_ATP_rate": 2.1e6,  # 1/ns
        "arp_bind_ADP_rate": 7.0e5,  # 1/ns
        "arp_unbind_ATP_rate": 1.4e-14,  # 1/ns
        "arp_unbind_ADP_rate": 8.0e-14,  # 1/ns
        "barbed_growth_branch_ATP_rate": 2.1e6,  # 1/ns
        "barbed_growth_branch_ADP_rate": 7.0e5,  # 1/ns
        "debranching_ATP_rate": 1.4e-14,  # 1/ns
        "debranching_ADP_rate": 8.0e-14,  # 1/ns
        "cap_bind_rate": 2.1e6,  # 1/ns
        "cap_unbind_rate": 1.4e-14,  # 1/ns
        "hydrolysis_actin_rate": 3.5e-15,  # 1/ns
        "hydrolysis_arp_rate": 3.5e-15,  # 1/ns
        "nucleotide_exchange_actin_rate": 1e-10,  # 1/ns
        "nucleotide_exchange_arp_rate": 1e-10,  # 1/ns
        "verbose": False,
    }

    def __init__(self, parameters=None):
        super(ReaddyActinProcess, self).__init__(parameters)
        actin_simulation = ActinSimulation(self.parameters)
        self.readdy_system = actin_simulation.system
        self.readdy_simulation = actin_simulation.simulation

    def ports_schema(self):
        """
        ports_schema returns a dictionary that declares how each state will behave.
        Each key can be assigned settings for the schema_keys declared in Store:

        * `_default`
        * `_updater`
        * `_divider`
        * `_value`
        * `_properties`
        * `_emit`
        * `_serializer`
        """
        return {
            "topologies": {
                "*": {
                    "type": {
                        "_default": "",
                        "_updater": "set",
                        "_emit": True,
                    },
                    "particles": {
                        "*": {
                            "type": {
                                "_default": "",
                                "_updater": "set",
                                "_emit": True,
                            },
                            "position": {
                                "_default": np.zeros(3),
                                "_updater": "set",
                                "_emit": True,
                            },
                            "neighbors": {
                                "_default": [],
                                "_updater": "set",
                                "_emit": True,
                            },
                        },
                    },
                }
            }
        }

    @staticmethod
    def add_topologies_to_readdy(readdy_simulation, topologies):
        """
        Add the given topologies of particles to the ReaDDy simulation
        """
        for t in topologies:
            particles = topologies[t]["particles"]
            types = []
            positions = []
            for p in particles:
                types.append(particles[p]["type"])
                positions.append(particles[p]["position"])
            top = readdy_simulation.add_topology(
                topologies[t]["type"], types, np.array(positions)
            )
            pid0 = -1
            for pid, p in particles.items():
                if pid0 < 0:
                    pid0 = int(pid)
                for n in range(len(particles[pid]["neighbors"])):
                    top.get_graph().add_edge(
                        int(pid) - pid0, int(p["neighbors"][n]) - pid0
                    )

    @staticmethod
    def get_readdy_particle_edges(readdy_topologies):
        """
        get all the edges in the ReaDDy topologies
        as (particle1 id, particle2 id)
        """
        result = []
        for top in readdy_topologies:
            for v1, v2 in top.graph.edges:
                p1_id = top.particle_id_of_vertex(v1)
                p2_id = top.particle_id_of_vertex(v2)
                if p1_id <= p2_id:
                    result.append((p1_id, p2_id))
        return result

    @staticmethod
    def get_topologies_from_readdy(readdy_topologies):
        """
        get data for topologies of particles from ReaDDy topologies
        """
        edges = ReaddyActinProcess.get_readdy_particle_edges(readdy_topologies)
        result = {}
        for t in range(len(readdy_topologies)):
            particles = {}
            for p in readdy_topologies[t].particles:
                neighbor_ids = []
                for edge in edges:
                    if p.id == edge[0]:
                        neighbor_ids.append(edge[1])
                    elif p.id == edge[1]:
                        neighbor_ids.append(edge[0])
                particles[p.id] = {
                    "type": p.type,
                    "position": p.pos,
                    "neighbors": neighbor_ids,
                }
            result[t] = {"type": readdy_topologies[t].type, "particles": particles}
        return result

    def next_update(self, timestep, states):

        # set state of ReaDDy particles
        ReaddyActinProcess.add_topologies_to_readdy(
            self.readdy_simulation, states["topologies"]
        )

        # simulate in ReaDDy for the given timestep
        def loop():

            readdy_actions = self.readdy_simulation._actions
            init = readdy_actions.initialize_kernel()
            diffuse = readdy_actions.integrator_euler_brownian_dynamics(
                self.parameters["timestep"]
            )
            calculate_forces = readdy_actions.calculate_forces()
            create_nl = readdy_actions.create_neighbor_list(
                self.readdy_system.calculate_max_cutoff().magnitude
            )
            update_nl = readdy_actions.update_neighbor_list()
            react = readdy_actions.reaction_handler_uncontrolled_approximation(
                self.parameters["timestep"]
            )
            observe = readdy_actions.evaluate_observables()

            init()
            create_nl()
            calculate_forces()
            update_nl()
            observe(0)

            n_steps = int(timestep * 1e9 / self.parameters["timestep"])
            for t in tqdm(range(1, n_steps + 1)):
                diffuse()
                update_nl()
                react()
                update_nl()
                calculate_forces()
                observe(t)

        self.readdy_simulation._run_custom_loop(loop)

        # return an update that mirrors the ports structure
        topologies = ReaddyActinProcess.get_topologies_from_readdy(
            self.readdy_simulation.current_topologies
        )
        return {
            "topologies": topologies,
        }

    # functions to configure and run the process
    def run_readdy_actin_process():
        """
        Run a simulation of the process.

        Returns:
            The simulation output.
        """
        # initialize the process
        readdy_actin_process = ReaddyActinProcess({})

        # declare the initial state, mirroring the ports structure
        initial_state = {
            "topologies": {
                0: {
                    "type": "Actin-Monomer",
                    "particles": {
                        0: {
                            "type": "actin#free_ATP",
                            "position": np.array([2, 0, 0]),
                            "neighbors": [],
                        }
                    },
                },
                1: {
                    "type": "Arp23-Dimer",
                    "particles": {
                        1: {
                            "type": "arp2",
                            "position": np.array([0, 0, 0]),
                            "neighbors": [2],
                        },
                        2: {
                            "type": "arp3#ATP",
                            "position": np.array([0, 0, 4]),
                            "neighbors": [1],
                        },
                    },
                },
            }
        }

        # run the simulation
        sim_settings = {
            "total_time": 0.000000005,  # 1e3 steps
            "initial_state": initial_state,
        }
        output = simulate_process(readdy_actin_process, sim_settings)
        return output


def main():
    """Simulate the process and plot results."""
    # make an output directory to save plots
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    output = ReaddyActinProcess.run_readdy_actin_process()

    # plot the simulation output
    plot_settings = {}
    plot_simulation_output(output, plot_settings, out_dir)


if __name__ == "__main__":
    main()
