import os
import numpy as np

from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process,
    PROCESS_OUT_DIR,
)
from vivarium.plots.simulation_output import plot_simulation_output

from tqdm import tqdm

from jinja2 import Environment, PackageLoader, select_autoescape
env = Environment(
    loader=PackageLoader("vivarium_models"),
    autoescape=select_autoescape()
)

from pathlib import Path
import subprocess

NAME = "MEDYAN"


def filament_to_string(type, points):
    begin = ' '.join([
        str(element) for element in points[0]])
    end = ' '.join([
        str(element) for element in points[-1]])

    line = ' '.join(['FILAMENT', str(type), begin, end])

    return line


class MedyanProcess(Process):
    """
    MEDYAN
    """

    name = NAME

    defaults = {
        'input_directory': 'in/filaments',
        'output_directory': 'out/filaments',
        'system_template': 'filament-system.txt',
        'system_file': 'filament-system.txt',
        'filament_file': 'filaments.txt',
    }

    def __init__(self, parameters=None):
        super(MedyanProcess, self).__init__(parameters)

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
            'filaments': {
                '*': {
                    'type': {
                        '_default': 0},
                    'points': {
                        '_default': [] # list of shape (3) numpy arrays
                        }}}}

    def initial_state(self, config):
        # TODO: make this more general

        initial_state = {
            id: {
                'type': 0,
                'points': points}
            for id, points in config.get('points', {}).items()
        }

        return initial_state

    def next_update(self, timestep, state):
        initial_filaments = state['filaments']

        filament_lines = [
            filament_to_string(filament['type'], filament['points'])
            for filament in initial_filaments.values()]

        input_directory = Path(self.parameters['input_directory'])

        filament_text = '\n'.join(filament_lines)

        filament_path = input_directory / self.parameters[
            'filament_file']
        with open(filament_path, 'w') as filament_file:
            filament_file.write(filament_text)
        
        system_template = self.parameters['system_template']
        template = env.get_template(system_template)
        system_text = template.render(filament_file=filament_path)

        system_path = input_directory / self.parameters[
            'system_file']
        with open(system_path, 'w') as system_file:
            system_file.write(system_text)

        medyan_command = [
            self.parameters['medyan_executable'],
            '-s',
            system_file,
            '-i',
            input_directory,
            '-o',
            self.parameters['output_directory']]

        medyan_process = subprocess.Popen(
            medyan_command, stdout=subprocess.PIPE)
        output, error = medyan_process.communicate()

        


def main():
    # """Simulate the process and plot results."""
    # # make an output directory to save plots
    # out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    # if not os.path.exists(out_dir):
    #     os.makedirs(out_dir)

    # output = MedyanProcess.next_update()

    # # plot the simulation output
    # plot_settings = {}
    # plot_simulation_output(output, plot_settings, out_dir)


if __name__ == "__main__":
    main()
