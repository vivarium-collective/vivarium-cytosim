import numpy as np
import argparse

from vivarium.core.composer import Composer
from vivarium.core.engine import Engine
from vivarium.processes.alternator import Alternator, PeriodicEvent

from vivarium_models.processes.readdy_actin_process import ReaddyActinProcess
from vivarium_models.processes.medyan import MedyanProcess
from vivarium_models.processes.monomer_to_fiber import MonomerToFiber
from vivarium_models.processes.fiber_to_monomer import FiberToMonomer

READDY_TIMESTEP = 0.0000001
ALTERNATOR_PERIODS = [10.0, READDY_TIMESTEP]


class ActinFiber(Composer):
    defaults = {
        "periodic_event": {"periods": ALTERNATOR_PERIODS},
        "readdy_actin": {
            "time_step": READDY_TIMESTEP,
            "_condition": ("choices", "readdy_active"),
        },
        "medyan": {"time_step": 5.0, "_condition": ("choices", "medyan_active")},
        "fiber_to_monomer": {"_condition": ("choices", "medyan_active")},
        "monomer_to_fiber": {"_condition": ("choices", "readdy_active")},
        "alternator": {"choices": ["medyan_active", "readdy_active"]},
    }

    def __init__(self, config=None):
        super().__init__(config)

    def generate_processes(self, config):
        periodic_event = PeriodicEvent(config["periodic_event"])
        readdy = ReaddyActinProcess(config["readdy_actin"])
        medyan = MedyanProcess(config["medyan"])
        fiber_to_monomer = FiberToMonomer(config["fiber_to_monomer"])
        monomer_to_fiber = MonomerToFiber(config["monomer_to_fiber"])
        alternator = Alternator(config["alternator"])

        return {
            "periodic_event": periodic_event,
            "readdy_actin": readdy,
            "medyan": medyan,
            "fiber_to_monomer": fiber_to_monomer,
            "monomer_to_fiber": monomer_to_fiber,
            "alternator": alternator,
        }

    def generate_topology(self, config):
        return {
            "periodic_event": {
                "event_trigger": ("alternate_trigger",),
                "period_index": ("period_index",),
            },
            "readdy_actin": {"monomers": ("monomers",)},
            "medyan": {"fibers": ("fibers",)},
            "fiber_to_monomer": {"fibers": ("fibers",), "monomers": ("monomers",)},
            "monomer_to_fiber": {"fibers": ("fibers",), "monomers": ("monomers",)},
            "alternator": {
                "alternate_trigger": ("alternate_trigger",),
                "choices": {
                    "medyan_active": (
                        "choices",
                        "medyan_active",
                    ),
                    "readdy_active": (
                        "choices",
                        "readdy_active",
                    ),
                },
            },
        }


def test_actin_fiber():
    parser = argparse.ArgumentParser(description="Run a MEDYAN simulation")
    parser.add_argument(
        "medyan_executable_path",
        help="the file path to the MEDYAN executable",
    )
    args = parser.parse_args()
    initial_state = {
        "choices": {"medyan_active": True, "readdy_active": False},
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

    medyan_config = {
        "medyan_executable": args.medyan_executable_path,  # "...../medyan/build/medyan"
        "transform_points": [0, 0, 0],
    }

    actin_fiber_config = {"medyan": medyan_config}
    actin_fiber = ActinFiber(actin_fiber_config)

    composite = actin_fiber.generate()
    composite["initial_state"] = initial_state

    engine = Engine(
        processes=composite["processes"],
        topology=composite["topology"],
        initial_state=composite["initial_state"],
        emitter="simularium",
        emit_processes=True,
    )

    engine.update(15)
    engine.emitter.get_data()


if __name__ == "__main__":
    test_actin_fiber()
