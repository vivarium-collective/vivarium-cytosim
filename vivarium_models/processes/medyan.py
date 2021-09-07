import os
import numpy as np

from vivarium.core.engine import pf, pp
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

NAME = "MEDYAN"


def filament_to_string(type, points):
    point_strs = [
        " ".join([str(element) for element in point])
        for point in points]
    line = " ".join(["FILAMENT", str(type)] + point_strs)
    return line

def read_coordinates(coordinates_line):
    coordinate_strs = coordinates_line.split(' ')
    coordinates = []
    for n in range(0, len(coordinate_strs), 3):
        point = coordinate_strs[n:n+3]
        coordinates.append(np.array([
            float(p)
            for p in point]))
    return coordinates

def read_filament(filament_line, coordinates_line):
    _, id, type, length, delta_l, delta_r = filament_line.split(' ')
    coordinates = read_coordinates(coordinates_line)
    return {
        id: {
            'type_name': type,
            'points': coordinates}}


class MedyanProcess(Process):
    """
    MEDYAN
    """

    name = NAME

    defaults = {
        "input_directory": "in/filaments",
        "output_directory": "out/filaments",
        "system_template": "filament-system.txt",
        "system_file": "filament-system.txt",
        "filament_file": "filaments.txt",
        "medyan_executable": "medyan",
        "snapshot": 1.0,
        "tranform_bounds": np.array([0, 0, 0]),
        # TODO: provide a way to parameterize type name,
        #    translating between simulation type names and MEDYAN type indexes
    }

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
                    "points": {  # list of shape (3) numpy arrays
                        "_default": [],
                        "_updater": "set",
                        "_emit": True,
                    },
                }
            }
        }

    def initial_state(self, config):
        initial_filaments = {
            unique_id: {"type_name": "A", "points": points}
            for unique_id, points in config.get("filaments", {}).items()
        }
        return {
            'filaments': initial_filaments}

    def transform_points(self, points, inverse=False):
        transform = np.array(self.parameters['transform_points'])
        if inverse:
            transform = -transform
        return [
            transform + point
            for point in points]

    def transform_filament(self, filament, inverse=False):
        filament['points'] = self.transform_points(filament['points'], inverse)
        return filament

    def read_snapshot(self, snapshot_path):
        # TODO: read only the last timepoint for each filament

        with open(snapshot_path, "r") as snapshot:
            snapshot_lines = snapshot.read().split("\n")

        index = 0
        filaments = {}
        while index < len(snapshot_lines):
            line = snapshot_lines[index]
            if "FILAMENT" in line:
                # TODO: parse line and pull out filament id and type
                coordinates_line = snapshot_lines[index + 1]
                filament = read_filament(line, coordinates_line)
                filaments.update(filament)
                index += 2
            else:
                index += 1

        return filaments

    def next_update(self, timestep, state):
        initial_filaments = state["filaments"]
        filament_ids = list(initial_filaments.keys())

        filament_types = set()
        for filament in initial_filaments.values():
            filament_types.add(filament['type_name'])
        num_filament_types = len(filament_types)

        filament_lines = [
            filament_to_string(
                filament["type_name"],
                self.transform_points(filament["points"]))
            for filament in initial_filaments.values()
        ]

        input_directory = Path(self.parameters["input_directory"])

        filament_text = "\n".join(filament_lines)

        filament_path = input_directory / self.parameters["filament_file"]
        with open(filament_path, "w") as filament_file:
            filament_file.write(filament_text)

        system_template = self.parameters["system_template"]
        template = env.get_template(system_template)
        system_text = template.render(
            filament_file=self.parameters['filament_file'],
            num_filament_types=num_filament_types,
            timestep=timestep,
            snapshot_time=self.parameters['snapshot'])

        system_path = input_directory / self.parameters["system_file"]
        with open(system_path, "w") as system_file:
            system_file.write(system_text)

        medyan_command = [
            self.parameters["medyan_executable"],
            "-s",
            system_file.name,
            "-i",
            str(input_directory),
            "-o",
            self.parameters["output_directory"],
        ]

        medyan_process = subprocess.Popen(medyan_command, stdout=subprocess.PIPE)
        output, error = medyan_process.communicate()

        print(output.decode('utf-8'))

        # TODO: perform the reverse transform for output points

        output_directory = Path(self.parameters['output_directory'])
        filaments = self.read_snapshot(
            output_directory / "snapshot.traj")

        filaments = {
            filament_ids[int(id)]: self.transform_filament(filament, inverse=True)
            for id, filament in filaments.items()}

        return {
            'filaments': filaments}
        



def main():
    """Simulate the process and plot results."""
    # make an output directory to save plots
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    medyan = MedyanProcess({
        'medyan_executable': '/home/youdonotexist/Downloads/medyan-4.2.0/build/medyan',
        'transform_points': [500, 500, 500],
        'time_step': 10.0})
    initial_state = {
        'filaments': {
            '1': {
                'type_name': '0', # 'A',
                'points': [
                    np.array([-70.0, 0.0, 100.0]),
                    np.array([10.0, 100.0, 0.0])]},
            '2': {
                'type_name': '0', # 'B',
                'points': [
                    np.array([-70.0, 100.0, 0.0]),
                    np.array([10.0, 0.0, 100.0])]}}}

    output = simulate_process(medyan, {
        'initial_state': initial_state,
        'total_time': 100,
        'return_raw_data': True})

    import ipdb; ipdb.set_trace()

    # plot the simulation output
    plot_settings = {}
    plot_simulation_output(output, plot_settings, out_dir)


if __name__ == "__main__":
    main()
