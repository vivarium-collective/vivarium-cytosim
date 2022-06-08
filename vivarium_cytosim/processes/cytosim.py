import os

from vivarium.core.process import Process
from vivarium.core.engine import Engine
import docker

from jinja2 import Environment, FileSystemLoader, select_autoescape

from pathlib import Path
from vivarium_cytosim.library.schema import fibers_schema
from vivarium_cytosim.data.fibers import initial_fibers

from simulariumio.cytosim import CytosimConverter, CytosimData, CytosimObjectInfo
from simulariumio import InputFileData, DisplayData, DISPLAY_TYPE

NAME = "CYTOSIM"

RELATIVE_MICRON = 0.001
BOUNDARY_BUFFER = 0.95


class CytosimProcess(Process):
    defaults = {
        "model_name": "cytosim-buckling",
        "internal_timestep": 0.001,
        "viscosity": 1,
        "actin_segmentation": 0.1,
        "temperature": 37,  # in Celsius
        "confine": None,
        "working_directory": "working/",
        "template_directory": "vivarium_cytosim/templates/",
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)
        self._jinja_environment = None
        self.working_path = Path(self.parameters["working_directory"]) / Path(
            self.parameters["model_name"]
        )
        if not os.path.exists(self.working_path):
            os.makedirs(self.working_path)
        self.check_pull_docker_image()

    def jinja_environment(self):
        if self._jinja_environment is None:
            self._jinja_environment = Environment(
                loader=FileSystemLoader(self.parameters["template_directory"]),
                autoescape=select_autoescape(),
            )
        return self._jinja_environment

    def check_pull_docker_image(self):
        client = docker.from_env()
        images = client.images.list("simularium/cytosim")
        if len(images) < 1:
            print("Downloading simularium/cytosim:latest from Docker Hub...")
            client.images.pull("simularium/cytosim")
        else:
            print("cytosim docker image already exists, skipping download")

    def ports_schema(self):
        ports = fibers_schema()
        return ports

    def initial_state(self, config):
        return {}

    def next_update(self, timestep, state):
        print("in cytosim process next update")

        self._render_template(state["fibers"], state["fibers_box_extent"], timestep)
        CytosimProcess._run_cytosim(self.working_path)
        report = CytosimProcess._load_report(self.working_path)

        return {"fibers": report}

    @staticmethod
    def _fiber_section(id, fiber):
        points = fiber["points"]
        # convert to microns
        point_strs = [
            " ".join([str(element * RELATIVE_MICRON) for element in point])
            for point in points
        ]
        return {"id": id, "points": point_strs}

    def _render_template(self, init_fibers, fibers_box_extent, timestep):
        fiber_sections = [
            CytosimProcess._fiber_section(id, fiber)
            for id, fiber in init_fibers.items()
        ]
        box_extent = fibers_box_extent * RELATIVE_MICRON * BOUNDARY_BUFFER
        system_template = self.parameters["model_name"] + ".cym"
        template = self.jinja_environment().get_template(system_template)
        cytosim_config = template.render(
            internal_timestep=self.parameters["internal_timestep"],
            kT=CytosimProcess._temperature_to_kT(self.parameters["temperature"]),
            viscosity=self.parameters["viscosity"],
            actin_segmentation=self.parameters["actin_segmentation"],
            confine=self.parameters["confine"],
            bounds_x=box_extent[0],
            bounds_y=box_extent[1],
            bounds_z=box_extent[2],
            filaments=fiber_sections,
            simulation_time=int(timestep / self.parameters["internal_timestep"]),
        )
        config_path = self.working_path / "config.cym"
        with open(config_path, "w") as cytosim_file:
            cytosim_file.write(cytosim_config)

    @staticmethod
    def _run_cytosim(working_path):
        abs_working_path = os.path.abspath(working_path)
        client = docker.from_env()
        container = client.containers.run(
            image="simularium/cytosim:latest",
            name="cytosim-container",
            volumes=[
                f"{abs_working_path}:/home/",
            ],
            detach=True,
        )
        container.wait()  # block until container run is complete
        logs = container.logs().decode("utf-8")  # get container logs
        container.remove()  # remove the container
        print(logs)

    @staticmethod
    def _load_report(working_path):
        report_path = working_path / "fiber_points.txt"
        data = CytosimData(
            object_info={
                "fibers": CytosimObjectInfo(
                    cytosim_file=InputFileData(file_path=report_path),
                    display_data={
                        1: DisplayData(
                            name="Actin-Polymer", display_type=DISPLAY_TYPE.FIBER
                        )
                    },
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

    @staticmethod
    def _temperature_to_kT(temp):
        return 1.3806e-5 * (temp + 273.15)


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
    engine.update(3.0)
    engine.emitter.get_data()


if __name__ == "__main__":
    main()
