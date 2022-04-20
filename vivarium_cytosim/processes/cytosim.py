import os
import shutil

from vivarium.core.process import Process
from vivarium.core.engine import Engine

from jinja2 import Environment, FileSystemLoader, select_autoescape

from pathlib import Path
import subprocess
from vivarium_cytosim.library.schema import fibers_schema
from vivarium_cytosim.data.fibers import initial_fibers
import vivarium_cytosim.library.go_sim as go_sim

from simulariumio.cytosim import CytosimConverter, CytosimData, CytosimObjectInfo
from simulariumio import InputFileData, DisplayData

NAME = "CYTOSIM"

RELATIVE_MICRON = 0.001
BOUNDARY_BUFFER = 0.95


def fiber_section(id, fiber):
    points = fiber["points"]
    # convert to microns
    point_strs = [
        " ".join([str(element * RELATIVE_MICRON) for element in point])
        for point in points
    ]
    return {"id": id, "points": point_strs}

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
            "fibers": CytosimObjectInfo(
                cytosim_file=InputFileData(file_contents=output),
                display_data={1: DisplayData(name="Actin-Polymer")},
            )
        }
    )
    converter = CytosimConverter(data)
    all_ids = converter._data.agent_data.unique_ids[-1]
    all_fibers = converter._data.agent_data.subpoints[-1]
    all_n_points = converter._data.agent_data.n_subpoints[-1]
    all_types = converter._data.agent_data.types[-1]

    return {
        str(int(id)): {
            "type_name": fiber_type,
            "points": [
                points / RELATIVE_MICRON for points in fiber_points[: int(n_points)]
            ],
        }
        for id, fiber_points, n_points, fiber_type in zip(
            all_ids, all_fibers, all_n_points, all_types
        )
    }


def temperature_to_kT(temp):
    return 1.3806e-5 * (temp + 273.15)


# def random_actin_fiber(bounds):


class CytosimProcess(Process):
    defaults = {
        "model_name": "cytosim-buckling",
        "internal_timestep": 0.001,
        "viscosity": 1,
        "actin_segmentation": 0.1,
        "temperature": 37,  # in Celsius
        "cytosim_template": "cytosim-buckling.cym",
        "cell_radius": 5,
        "confine": None,
        "input_directory": "in/",
        "output_directory": "out/",
        "template_directory": "vivarium_cytosim/templates/",
        "cytosim_sim": "../cytosim/bin/sim",
        "cytosim_report": "../cytosim/bin/report",
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)

        self._jinja_environment = None
        self.input_path = Path(self.parameters["input_directory"])
        if not os.path.exists(self.input_path):
            os.makedirs(self.input_path)
        self.output_path = Path(self.parameters["output_directory"])
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)

    def jinja_environment(self):
        if self._jinja_environment is None:
            self._jinja_environment = Environment(
                loader=FileSystemLoader(self.parameters["template_directory"]),
                autoescape=select_autoescape(),
            )
        return self._jinja_environment

    def ports_schema(self):
        ports = fibers_schema()
        return ports

    def initial_state(self, config):
        return {}

    def next_update(self, timestep, state):
        initial_fibers = state["fibers"]

        fiber_sections = [
            fiber_section(id, fiber) for id, fiber in initial_fibers.items()
        ]

        box_extent = state["fibers_box_extent"] * RELATIVE_MICRON * BOUNDARY_BUFFER

        system_template = self.parameters["model_name"] + ".cym"
        template = self.jinja_environment().get_template(system_template)
        cytosim_config = template.render(
            internal_timestep=self.parameters["internal_timestep"],
            # radius=self.parameters['cell_radius'],
            kT=temperature_to_kT(self.parameters["temperature"]),
            viscosity=self.parameters["viscosity"],
            actin_segmentation=self.parameters["actin_segmentation"],
            confine=self.parameters["confine"],
            bounds_x=box_extent[0],
            bounds_y=box_extent[1],
            bounds_z=box_extent[2],
            filaments=fiber_sections,
            simulation_time=int(timestep / self.parameters["internal_timestep"]),
        )

        config_path = self.input_path / self.parameters["cytosim_template"]
        with open(config_path, "w") as cytosim_file:
            cytosim_file.write(cytosim_config)

        go_sim.local_run(
            os.path.abspath(self.parameters["cytosim_sim"]),
            config_path,
            self.parameters["model_name"],
            self.output_path,
        )

        previous_dir = os.getcwd()
        report_base = self.output_path / self.parameters["model_name"]
        report_exec = os.path.abspath(self.parameters["cytosim_report"])
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

        # import ipdb; ipdb.set_trace()

        return {"fibers": report}


def main():
    cytosim = CytosimProcess(
        {"confine": {"side": "inside", "force": 100, "space": "cell"}}
    )

    engine = Engine(
        processes={"cytosim": cytosim},
        topology={
            "cytosim": {
                "fibers": ("fibers",),
                "fibers_box_extent": ("fibers_box_extent",),
            }
        },
        initial_state=initial_fibers,
    )

    engine.update(10.0)
    engine.emitter.get_data()

    # import ipdb; ipdb.set_trace()


if __name__ == "__main__":
    main()
