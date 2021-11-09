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

NAME = "MEDYAN"


class MedyanProcess(Process):
    """
    MEDYAN
    """

    name = NAME

    defaults = {
        "model_name": "medyan_Chandrasekaran_2019_no_tread_2mUNI_alphaA_0.1_MA_0.675",
        "input_directory": "in/",
        "output_directory": "out/",
        "medyan_executable": "medyan",
        "snapshot": 1.0,
        "tranform_bounds": np.array([0, 0, 0]),
        # TODO: provide a way to parameterize type name,
        #    translating between simulation type names and MEDYAN type indexes
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)

        assert self.parameters["time_step"] > self.parameters["snapshot"]

    def ports_schema(self):
        return {
            "fibers_box_extent": {
                "_default": np.array([4000.0, 2000.0, 2000.0]),
                "_updater": "set",
                "_emit": True,
            },
            "fibers": {
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
            },
        }

    def initial_state(self, config):
        initial_fibers = {
            unique_id: {"type_name": "A", "points": points}
            for unique_id, points in config.get("fibers", {}).items()
        }
        return {
            "fibers_box_extent": np.array([4000.0, 2000.0, 2000.0]),
            "fibers": initial_fibers,
        }

    def transform_points(self, points, inverse=False):
        transform = np.array(self.parameters["transform_points"])
        if inverse:
            transform = -transform
        return [transform + point for point in points]

    def transform_fiber(self, fiber, inverse=False):
        fiber["points"] = self.transform_points(fiber["points"], inverse)
        return fiber

    def read_snapshot(self, snapshot_path):
        # TODO: read only the last timepoint for each fiber

        with open(snapshot_path, "r") as snapshot:
            snapshot_lines = snapshot.read().split("\n")

        index = 0
        fibers = {}
        while index < len(snapshot_lines):
            line = snapshot_lines[index]
            if "FILAMENT" in line:
                # TODO: parse line and pull out fiber id and type
                coordinates_line = snapshot_lines[index + 1]
                fiber = MedyanProcess.read_fiber(line, coordinates_line)
                fibers.update(fiber)
                index += 2
            else:
                index += 1

        return fibers

    def next_update(self, timestep, state):
        print("in medyan process next update")

        initial_fibers = state["fibers"]
        fiber_ids = list(initial_fibers.keys())

        fiber_lines = [
            MedyanProcess.fiber_to_string(
                fiber["type_name"], self.transform_points(fiber["points"])
            )
            for fiber in initial_fibers.values()
        ]

        input_directory = Path(self.parameters["input_directory"]) / Path(
            self.parameters["model_name"]
        )

        output_directory = Path(self.parameters["output_directory"]) / Path(
            self.parameters["model_name"]
        )
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        fiber_text = "\n".join(fiber_lines)

        fiber_path = input_directory / "filaments.txt"
        with open(fiber_path, "w") as fiber_file:
            fiber_file.write(fiber_text)

        system_template = self.parameters["model_name"] + ".txt"
        template = env.get_template(system_template)
        system_text = template.render(
            timestep=timestep,
            snapshot_time=self.parameters["snapshot"],
        )
        box_extent = MedyanProcess.read_box_extent(system_text)

        system_path = input_directory / (self.parameters["model_name"] + ".txt")
        with open(system_path, "w") as system_file:
            system_file.write(system_text)

        medyan_command = [
            self.parameters["medyan_executable"],
            "-s",
            system_file.name,
            "-i",
            str(input_directory),
            "-o",
            Path(self.parameters["output_directory"])
            / Path(self.parameters["model_name"]),
        ]

        medyan_process = subprocess.Popen(medyan_command, stdout=subprocess.PIPE)
        output, error = medyan_process.communicate()

        print(output.decode("utf-8"))

        # TODO: perform the reverse transform for output points

        fibers = self.read_snapshot(output_directory / "snapshot.traj")

        fibers = {
            fiber_ids[int(id)]: self.transform_fiber(fiber, inverse=True)
            for id, fiber in fibers.items()
        }

        import ipdb; ipdb.set_trace()
        
        return {"fibers_box_extent": box_extent, "fibers": fibers}

    @staticmethod
    def fiber_to_string(type_name, points):
        point_strs = [" ".join([str(element) for element in point]) for point in points]
        line = " ".join(["FILAMENT", str(type_name)] + point_strs)
        return line

    @staticmethod
    def read_coordinates(coordinates_line):
        coordinate_strs = coordinates_line.strip().split(" ")
        coordinates = []
        for n in range(0, len(coordinate_strs), 3):
            point = coordinate_strs[n : n + 3]
            coordinates.append(np.array([float(p) for p in point]))
        return coordinates

    @staticmethod
    def read_fiber(fiber_line, coordinates_line):
        _, id, type, length, delta_l, delta_r = fiber_line.split(" ")
        coordinates = MedyanProcess.read_coordinates(coordinates_line)
        return {id: {"type_name": type, "points": coordinates}}

    @staticmethod
    def read_box_extent(system_text):
        lines = system_text.split("\n")
        n_compartments = np.zeros(3)
        compartment_size = np.zeros(3)
        coords = ["X", "Y", "Z"]
        for line in lines:
            for dim in range(len(coords)):
                coord = coords[dim]
                if f"N{coord}:" in line:
                    n_compartments[dim] = int(line.split()[1])
                if f"COMPARTMENTSIZE{coord}:" in line:
                    compartment_size[dim] = float(line.split()[1])
        return np.multiply(n_compartments, compartment_size)


def main():
    """Simulate the process and plot results."""
    parser = argparse.ArgumentParser(description="Run a MEDYAN simulation")
    parser.add_argument(
        "medyan_executable_path",
        help="the file path to the MEDYAN executable",
    )
    args = parser.parse_args()

    # make an output directory to save plots
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    medyan = MedyanProcess(
        {
            "medyan_executable": args.medyan_executable_path,
            "transform_points": [500, 500, 500],
            "time_step": 10.0,
        }
    )
    initial_state = {
        "fibers_box_extent": np.array([4000.0, 2000.0, 2000.0]),
        "fibers": {
            "1": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 912.50000000, 1000.00000000]),
                    np.array([3160.00000000, 912.50000000, 1000.00000000]),
                ],
            },
            "2": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 947.50000000, 939.37822174]),
                    np.array([3160.00000000, 947.50000000, 939.37822174]),
                ],
            },
            "3": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 930.00000000, 969.68911087]),
                    np.array([3160.00000000, 930.00000000, 969.68911087]),
                ],
            },
            "4": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 947.50000000, 1000.00000000]),
                    np.array([3160.00000000, 947.50000000, 1000.00000000]),
                ],
            },
            "5": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 930.00000000, 1030.31088913]),
                    np.array([3160.00000000, 930.00000000, 1030.31088913]),
                ],
            },
            "6": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 947.50000000, 1060.62177826]),
                    np.array([3160.00000000, 947.50000000, 1060.62177826]),
                ],
            },
            "7": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 965.00000000, 909.06733260]),
                    np.array([3160.00000000, 965.00000000, 909.06733260]),
                ],
            },
            "8": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 982.50000000, 939.37822174]),
                    np.array([3160.00000000, 982.50000000, 939.37822174]),
                ],
            },
            "9": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 965.00000000, 969.68911087]),
                    np.array([3160.00000000, 965.00000000, 969.68911087]),
                ],
            },
            "10": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 982.50000000, 1000.00000000]),
                    np.array([3160.00000000, 982.50000000, 1000.00000000]),
                ],
            },
            "11": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 965.00000000, 1030.31088913]),
                    np.array([3160.00000000, 965.00000000, 1030.31088913]),
                ],
            },
            "12": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 982.50000000, 1060.62177826]),
                    np.array([3160.00000000, 982.50000000, 1060.62177826]),
                ],
            },
            "13": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 965.00000000, 1090.93266740]),
                    np.array([3160.00000000, 965.00000000, 1090.93266740]),
                ],
            },
            "14": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1000.00000000, 909.06733260]),
                    np.array([3160.00000000, 1000.00000000, 909.06733260]),
                ],
            },
            "15": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1017.50000000, 939.37822174]),
                    np.array([3160.00000000, 1017.50000000, 939.37822174]),
                ],
            },
            "16": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1000.00000000, 969.68911087]),
                    np.array([3160.00000000, 1000.00000000, 969.68911087]),
                ],
            },
            "17": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1017.50000000, 1000.00000000]),
                    np.array([3160.00000000, 1017.50000000, 1000.00000000]),
                ],
            },
            "18": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1000.00000000, 1030.31088913]),
                    np.array([3160.00000000, 1000.00000000, 1030.31088913]),
                ],
            },
            "19": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1017.50000000, 1060.62177826]),
                    np.array([3160.00000000, 1017.50000000, 1060.62177826]),
                ],
            },
            "20": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1000.00000000, 1090.93266740]),
                    np.array([3160.00000000, 1000.00000000, 1090.93266740]),
                ],
            },
            "21": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1035.00000000, 909.06733260]),
                    np.array([3160.00000000, 1035.00000000, 909.06733260]),
                ],
            },
            "22": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1052.50000000, 939.37822174]),
                    np.array([3160.00000000, 1052.50000000, 939.37822174]),
                ],
            },
            "23": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1035.00000000, 969.68911087]),
                    np.array([3160.00000000, 1035.00000000, 969.68911087]),
                ],
            },
            "24": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1052.50000000, 1000.00000000]),
                    np.array([3160.00000000, 1052.50000000, 1000.00000000]),
                ],
            },
            "25": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1035.00000000, 1030.31088913]),
                    np.array([3160.00000000, 1035.00000000, 1030.31088913]),
                ],
            },
            "26": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1052.50000000, 1060.62177826]),
                    np.array([3160.00000000, 1052.50000000, 1060.62177826]),
                ],
            },
            "27": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1035.00000000, 1090.93266740]),
                    np.array([3160.00000000, 1035.00000000, 1090.93266740]),
                ],
            },
            "28": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1070.00000000, 969.68911087]),
                    np.array([3160.00000000, 1070.00000000, 969.68911087]),
                ],
            },
            "29": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1087.50000000, 1000.00000000]),
                    np.array([3160.00000000, 1087.50000000, 1000.00000000]),
                ],
            },
            "30": {
                "type_name": "Actin-Polymer",
                "points": [
                    np.array([1000.00000000, 1070.00000000, 1030.31088913]),
                    np.array([3160.00000000, 1070.00000000, 1030.31088913]),
                ],
            },
        },
    }

    output = simulate_process(
        medyan,
        {"initial_state": initial_state, "total_time": 100, "return_raw_data": True},
    )

    # plot the simulation output
    plot_settings = {}
    plot_simulation_output(output, plot_settings, out_dir)


if __name__ == "__main__":
    main()
