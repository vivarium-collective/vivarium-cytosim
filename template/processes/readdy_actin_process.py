import os

from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process_in_experiment,
    PROCESS_OUT_DIR,
)
from vivarium.plots.simulation_output import plot_simulation_output

import readdy
from readdy_util import ReaddyUtil
from actin_util import *


NAME = 'ReaDDy_actin'

class ReaddyActinProcess(Process):
    '''
    This process uses ReaDDy to model the dynamics 
    of branched actin networks made of spatially explicit monomers 
    '''
    name = READDY_ACTIN

    defaults = {
        'box_size' : 150., #nm
        'actin_concentration' : 250., #uM
        'arp23_concentration' : 10., #uM
        'cap_concentration' : 0., #uM
        'dimerize_rate' : 2.1e-2, #1/ns
        'dimerize_reverse_rate' : 1.4e-1, #1/ns
        'trimerize_rate' : 2.1e-2, #1/ns
        'trimerize_reverse_rate' : 1.4e-1, #1/ns
        'pointed_growth_ATP_rate' : 2.4e5, #1/ns
        'pointed_growth_ADP_rate' : (0.16 / 1.3) * 2.4e5, #1/ns
        'pointed_shrink_ATP_rate' : 0.8e-14, #1/ns
        'pointed_shrink_ADP_rate' : (0.3 / 0.8) * 0.8e-14, #1/ns
        'barbed_growth_ATP_rate' : 2.1e6, #1/ns
        'barbed_growth_ADP_rate' : (4. / 12.) * 2.1e6, #1/ns
        'nucleate_ATP_rate' : 2.1e6, #1/ns
        'nucleate_ADP_rate' : (4. / 12.) * 2.1e6, #1/ns
        'barbed_shrink_ATP_rate' : 1.4e-14, #1/ns
        'barbed_shrink_ADP_rate' : (8. / 1.4) * 1.4e-14, #1/ns
        'arp_bind_ATP_rate' : 2.1e6, #1/ns 
        'arp_bind_ADP_rate' : (4. / 12.) * 2.1e6, #1/ns
        'arp_unbind_ATP_rate' : 1.4e-14, #1/ns
        'arp_unbind_ADP_rate' : (8. / 1.4) * 1.4e-14, #1/ns
        'barbed_growth_branch_ATP_rate' : 2.1e6, #1/ns
        'barbed_growth_branch_ADP_rate' : (4. / 12.) * 2.1e6, #1/ns
        'debranching_ATP_rate' : 1.4e-14, #1/ns
        'debranching_ADP_rate' : (8. / 1.4) * 1.4e-14, #1/ns
        'cap_bind_rate' : 2.1e6, #1/ns
        'cap_unbind_rate' : 1.4e-14, #1/ns
        'hydrolysis_actin_rate' : 3.5e-15, #1/ns
        'hydrolysis_arp_rate' : 3.5e-15, #1/ns
        'nucleotide_exchange_actin_rate' : 1e-10, #1/ns
        'nucleotide_exchange_arp_rate' : 1e-10, #1/ns
        'temperature_C' : 22., # from Pollard experiments
        'eta' : 8.1, #cP, viscosity in cytoplasm
        'force_constant' : 9000.,
        'reaction_distance' : 1., #nm
        'actin_radius' : 2., #nm
        'arp23_radius' : 2., #nm
        'cap_radius' : 3., #nm
        'timestep' : 0.005, #ns
        'n_cpu' : 4,
        'init_particles' : None  # TODO should this go here?
    }

    def __init__(self, parameters=None):
        super(ReaddyActinProcess, self).__init__(parameters)
        actin_util = ActinUtil()
        self.readdy_simulation, self.readdy_system = actin_util.create_actin_simulation(
            self.parameters)
        if self.parameters['init_particles'] is None:
            # add random monomers
            actin_util.add_free_monomers(self.parameters, self.readdy_simulation)
        else:
            # TODO add specific initial particles

    def ports_schema(self):
        '''
        ports_schema returns a dictionary that declares how each state will behave.
        Each key can be assigned settings for the schema_keys declared in Store:

        * `_default`
        * `_updater`
        * `_divider`
        * `_value`
        * `_properties`
        * `_emit`
        * `_serializer`
        '''
        return {
            'internal': {
                # TODO how to represent data shaped like this? or shape it differently?
                # {p_id : (p_type, neighbor_ids, np.array([p_pos[0], p_pos[1], p_pos[2]]))}
                'particles': {
                    '_default': {},
                    '_updater': 'set',
                },
            },
        }

    def derivers(self):
        '''
        declare which derivers are needed for this process
        '''
        return {}

    def next_update(self, timestep, states):

        # get the states
        particles = states['internal']['particles']
        # TODO update particles in ReaDDy as needed to match states
        
        # calculate timestep-dependent updates
        def loop():

            init = self.readdy_simulation._actions.initialize_kernel()
            diffuse = self.readdy_simulation._actions.integrator_euler_brownian_dynamics(
                self.parameters['timestep'])
            calculate_forces = self.readdy_simulation._actions.calculate_forces()
            create_nl = self.readdy_simulation._actions.create_neighbor_list(
                self.readdy_system.calculate_max_cutoff())
            update_nl = self.readdy_simulation._actions.update_neighbor_list()
            react = self.readdy_simulation._actions.reaction_handler_uncontrolled_approximation(
                self.parameters['timestep'])
            observe = self.readdy_simulation._actions.evaluate_observables()

            init()
            create_nl()
            calculate_forces()
            update_nl()
            observe(0)

            n_steps = timestep * 1e9 / self.parameters['timestep']
            for t in tqdm(range(1, n_steps + 1)):
                diffuse()
                update_nl()
                react()
                update_nl()
                calculate_forces()
                observe(t)

        self.readdy_simulation._run_custom_loop(loop)

        # return an update that mirrors the ports structure
        observables = self.readdy_simulation.observe()
        # TODO update the particles port using observables from ReaDDy
        return {
            'internal': {
                'particles': particles,
            }
        }


# functions to configure and run the process
def run_readdy_actin_process():
    '''
    Run a simulation of the process.

    Returns:
        The simulation output.
    '''
    # initialize the process by passing in parameters
    parameters = {}
    readdy_actin_process = ReaddyActinProcess(parameters)

    # declare the initial state, mirroring the ports structure
    initial_state = {
        'internal': {
            'particles': {}
        },
    }

    # run the simulation
    sim_settings = {
        'total_time': 0.0000005, # 1e5 steps
        'initial_state': initial_state
    }
    output = simulate_process_in_experiment(readdy_actin_process, sim_settings)

    return output


def test_readdy_actin_process():
    '''Test that the process runs correctly.

    This will be executed by pytest.
    '''
    output = run_readdy_actin_process()
    # TODO: Add assert statements to ensure correct performance.


def main():
    '''Simulate the process and plot results.'''
    # make an output directory to save plots
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    output = run_readdy_actin_process()

    # plot the simulation output
    plot_settings = {}
    plot_simulation_output(output, plot_settings, out_dir)


if __name__ == '__main__':
    main()
