import os
import numpy as np
import argparse

from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process,
    PROCESS_OUT_DIR,
)
from vivarium.plots.simulation_output import plot_simulation_output

from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
    loader=PackageLoader("vivarium_models"), autoescape=select_autoescape()
)

from pathlib import Path
import subprocess
from vivarium_models.library.schema import fibers_schema
from vivarium_models.data.fibers import initial_fibers

NAME = "CYTOSIM"

def fiber_section(id, fiber):
    points = fiber['points']
    relative_micron = 0.001
    # convert to microns
    point_strs = [" ".join([str(element * relative_micron) for element in point]) for point in points]
    point_line = ', '.join(point_strs)

    lines = [
        f'new actin',
        '{',
        f'    mark = {id}',
        f'    shape = {point_line}',
        '}',
    ]

    return '\n'.join(lines)


class CytosimProcess(Process):
    defaults = {
        'internal_timestep': 0.001,
        'cytosim_template': "cytosim-buckling.cym",
        'cell_radius': 5,
        "input_directory": "in/",
        "output_directory": "out/",
        'cytosim_executable': "../cytosim/bin/sim"}

    def __init__(self, parameters=None):
        super().__init__(parameters)

    def ports_schema(self):
        return fibers_schema()

    def initial_state(self, config):
        return {}

    def next_update(self, timestep, state):
        initial_fibers = state["fibers"]
        fiber_ids = list(initial_fibers.keys())

        fiber_sections = [fiber_section(id, fiber) for id, fiber in initial_fibers.items()]

        template = env.get_template(self.parameters['cytosim_template'])
        cytosim_config = template.render(
            internal_timestep=self.parameters['internal_timestep'],
            radius=self.parameters['cell_radius'],
            filament_section='\n\n\n'.join(fiber_sections),
            simulation_time=int(timestep/self.parameters['internal_timestep']),
        )

        import ipdb; ipdb.set_trace()

        return {}

def main():
    cytosim = CytosimProcess({})
    output = simulate_process(
        cytosim, {
            'initial_state': initial_fibers,
            'total_time': 1,
            'return_raw_data': True
        }
    )

if __name__ == '__main__':
    main()
