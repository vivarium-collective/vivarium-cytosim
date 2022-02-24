import os
import shutil
import numpy as np
import argparse

from vivarium.core.process import Process
from vivarium.core.engine import Engine
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
import vivarium_models.library.go_sim as go_sim

from simulariumio.cytosim import CytosimConverter, CytosimData, CytosimObjectInfo
from simulariumio import InputFileData, DisplayData

NAME = "CYTOSIM"

RELATIVE_MICRON = 0.001

def fiber_section(id, fiber):
    points = fiber['points']
    # convert to microns
    point_strs = [" ".join([str(element * RELATIVE_MICRON) for element in point]) for point in points]
    return {
        'id': id,
        'points': point_strs}


    # point_line = ', '.join(point_strs)

    # lines = [
    #     f'new actin',
    #     '{',
    #     f'    mark = {id}',
    #     f'    shape = {point_line}',
    #     '}',
    # ]

    # return '\n'.join(lines)


def load_report(output):
    data = CytosimData(
        object_info={
            'fibers': CytosimObjectInfo(
                cytosim_file=InputFileData(
                    file_contents=output),
                display_data={
                    1: DisplayData(name='Actin-Polymer')})})
    converter = CytosimConverter(data)
    all_ids = converter._data.agent_data.unique_ids[-1]
    all_fibers = converter._data.agent_data.subpoints[-1]
    all_n_points = converter._data.agent_data.n_subpoints[-1]
    all_types = converter._data.agent_data.types[-1]

    return {
        str(int(id)): {
            'type_name': fiber_type,
            'points': [points / RELATIVE_MICRON for points in fiber_points[:int(n_points)]]}
        for id, fiber_points, n_points, fiber_type in zip(all_ids, all_fibers, all_n_points, all_types)}


class CytosimProcess(Process):
    defaults = {
        'model_name': 'cytosim_model',
        'internal_timestep': 0.001,
        'cytosim_template': "cytosim-buckling.cym",
        'cell_radius': 5,
        'confine': None,
        "input_directory": "in/",
        "output_directory": "out/",
        'cytosim_sim': "../cytosim/bin/sim",
        'cytosim_report': "../cytosim/bin/report"}

    def __init__(self, parameters=None):
        super().__init__(parameters)

        self.input_path = Path(self.parameters['input_directory'])
        self.output_path = Path(self.parameters['output_directory'])

    def ports_schema(self):
        ports = fibers_schema()
        ports['choices'] = {
            'medyan_active': {
                '_default': True,
                '_emit': True},
            'readdy_active': {
                '_default': False,
                '_emit': True}}

        return ports

    def calculate_timestep(self, state):
        return 0.1

    def initial_state(self, config):
        return {}

    def next_update(self, timestep, state):
        initial_fibers = state["fibers"]
        fiber_ids = list(initial_fibers.keys())

        fiber_sections = [fiber_section(id, fiber) for id, fiber in initial_fibers.items()]

        box_extent = state['fibers_box_extent'] * RELATIVE_MICRON

        template = env.get_template(self.parameters['cytosim_template'])
        cytosim_config = template.render(
            internal_timestep=self.parameters['internal_timestep'],
            # radius=self.parameters['cell_radius'],
            confine=self.parameters['confine'],
            bounds_x=box_extent[0],
            bounds_y=box_extent[1],
            bounds_z=box_extent[2],
            filaments=fiber_sections,
            # filament_section='\n\n\n'.join(fiber_sections),
            simulation_time=int(timestep/self.parameters['internal_timestep']),
        )

        config_path = self.input_path / self.parameters['cytosim_template']
        with open(config_path, "w") as cytosim_file:
            cytosim_file.write(cytosim_config)

        go_sim.local_run(
            os.path.abspath(self.parameters['cytosim_sim']),
            config_path,
            self.parameters['model_name'],
            self.output_path,
        )

        previous_dir = os.getcwd()
        report_base = self.output_path / self.parameters['model_name']
        report_exec = os.path.abspath(self.parameters['cytosim_report'])
        os.chdir(report_base)

        report_command = [
            report_exec,
            "fiber:points",
        ]

        cytosim_process = subprocess.Popen(report_command, stdout=subprocess.PIPE)
        output, error = cytosim_process.communicate()

        os.chdir(previous_dir)
        shutil.rmtree(report_base)

        print(output.decode("utf-8"))
        report = load_report(output.decode("utf-8"))

        import ipdb; ipdb.set_trace()

        return {'fibers': report}


def main():
    cytosim = CytosimProcess({
        'confine': {
            'side': 'inside',
            'force': 100,
            'space': 'cell'}})

    engine = Engine(
        processes={'cytosim': cytosim},
        topology={
            'cytosim': {
                'fibers': ('fibers',),
                'fibers_box_extent': ('fibers_box_extent',),
                'choices': ('choices',)}},
        initial_state=initial_fibers,
        emitter='simularium')

    engine.update(10.0)
    output = engine.emitter.get_data()

    import ipdb; ipdb.set_trace()

if __name__ == '__main__':
    main()
