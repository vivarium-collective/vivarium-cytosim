import numpy as np

from vivarium.core.process import Process
from vivarium.core.engine import Engine, pf
from vivarium.core.composition import simulate_process

from tqdm import tqdm
from simularium_models_util.actin import ActinSimulation, ActinUtil, ActinTestData
from simularium_models_util import ReaddyUtil
from vivarium_models.util import create_monomer_update


NAME = "ReaDDy_actin"

test_monomer_data = {
    "monomers": {
        "topologies": {
            1: {
                "type_name": "Arp23-Dimer",
                "particle_ids": [1, 2],
            },
            0: {
                "type_name": "Actin-Monomer",
                "particle_ids": [0],
            },
        },
        "particles": {
            0: {
                "type_name": "actin#free_ATP",
                "position": np.array([2, 0, 0]),
                "neighbor_ids": [],
            },
            1: {
                "type_name": "arp2",
                "position": np.array([0, 0, 0]),
                "neighbor_ids": [2],
            },
            2: {
                "type_name": "arp3#ATP",
                "position": np.array([0, 0, 4]),
                "neighbor_ids": [1],
            },
        },
    }
}


class ReaddyActinProcess(Process):
    """
    This process uses ReaDDy to model the dynamics
    of branched actin networks made of spatially explicit monomers
    """

    name = NAME

    defaults = {
        "name": "actin",
        "total_steps": 1e3,
        "time_step": 0.0000001,
        "internal_timestep": 0.1,
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
        "use_box_actin": False,
        "use_box_arp": False,
        "use_box_cap": False,
        "implicit_actin_concentration": 0,
        "nonspatial_polymerization": False,
        "verbose": False,
    }

    def __init__(self, parameters=None):
        super(ReaddyActinProcess, self).__init__(parameters)
        actin_simulation = ActinSimulation(self.parameters)
        self.readdy_system = actin_simulation.system
        self.readdy_simulation = actin_simulation.simulation

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
            }
        }

    def initial_state(self, config=None):
        # TODO: make this more general
        return test_monomer_data

    def simulate_readdy(self, timestep):
        """
        Simulate in ReaDDy for the given timestep
        """

        def loop():
            readdy_actions = self.readdy_simulation._actions
            init = readdy_actions.initialize_kernel()
            diffuse = readdy_actions.integrator_euler_brownian_dynamics(
                self.parameters["internal_timestep"]
            )
            calculate_forces = readdy_actions.calculate_forces()
            create_nl = readdy_actions.create_neighbor_list(
                self.readdy_system.calculate_max_cutoff().magnitude
            )
            update_nl = readdy_actions.update_neighbor_list()
            react = readdy_actions.reaction_handler_uncontrolled_approximation(
                self.parameters["internal_timestep"]
            )
            observe = readdy_actions.evaluate_observables()
            init()
            create_nl()
            calculate_forces()
            update_nl()
            observe(0)
            n_steps = int(timestep * 1e9 / self.parameters["internal_timestep"])
            for t in tqdm(range(1, n_steps + 1)):
                diffuse()
                update_nl()
                react()
                update_nl()
                calculate_forces()
                observe(t)

        self.readdy_simulation._run_custom_loop(loop)

    def next_update(self, timestep, states):
        print("in readdy actin process next update")

        ActinUtil.add_monomers_from_data(self.readdy_simulation, states["monomers"])
        self.simulate_readdy(timestep)
        readdy_monomers = ReaddyUtil.get_current_monomers(
            self.readdy_simulation.current_topologies
        )

        return create_monomer_update(states["monomers"], readdy_monomers)

    # functions to configure and run the process
    def run_readdy_actin_process():
        """
        Run a simulation of the process.

        Returns:
            The simulation output.
        """
        # initialize the process
        readdy_actin_process = ReaddyActinProcess({})

        # run the simulation
        sim_settings = {
            "total_time": 0.000000005,  # 50 steps
            "initial_state": test_monomer_data,
        }
        output = simulate_process(readdy_actin_process, sim_settings)
        return output


def test_readdy_actin_process():
    monomer_data = ActinTestData.linear_actin_monomers()
    readdy_actin_process = ReaddyActinProcess()

    engine = Engine(
        {
            "processes": {"readdy_actin_process": readdy_actin_process},
            "topology": {"readdy_actin_process": {"monomers": ("monomers",)}},
            "initial_state": {"monomers": monomer_data},
        }
    )

    engine.update(0.0000001)  # 1e3 steps

    output = engine.emitter.get_data()
    print(pf(output))


if __name__ == "__main__":
    test_readdy_actin_process()
