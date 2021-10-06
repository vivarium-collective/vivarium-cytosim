import numpy as np

from vivarium.core.composer import Composer
from vivarium.core.engine import Engine
from vivarium.processes.alternator import Alternator, PeriodicEvent

from vivarium_models.processes.readdy_actin_process import ReaddyActinProcess
from vivarium_models.processes.medyan import MedyanProcess
from vivarium_models.processes.monomer_to_fiber import MonomerToFiber
from vivarium_models.processes.fiber_to_monomer import FiberToMonomer
from vivarium_models.processes.simularium_emitter import SimulariumEmitter

READDY_TIMESTEP = 0.0000001
ALTERNATOR_PERIODS = [10.0, READDY_TIMESTEP]

class ActinFiber(Composer):
    defaults = {
        'periodic_event': {
            'periods': ALTERNATOR_PERIODS},
        'readdy_actin': {
            'time_step': READDY_TIMESTEP,
            '_condition': ('choices', 'readdy_active')},
        'medyan': {
            'time_step': 5.0,
            '_condition': ('choices', 'medyan_active')},
        'fiber_to_monomer': {
            '_condition': ('choices', 'medyan_active')},
        'monomer_to_fiber': {
            '_condition': ('choices', 'readdy_active')},
        'alternator': {
            'choices': ['medyan_active', 'readdy_active']},
    }

    def __init__(self, config=None):
        super().__init__(config)

    def generate_processes(self, config):
        periodic_event = PeriodicEvent(config['periodic_event'])
        readdy = ReaddyActinProcess(config['readdy_actin'])
        medyan = MedyanProcess(config['medyan'])
        fiber_to_monomer = FiberToMonomer(config['fiber_to_monomer'])
        monomer_to_fiber = MonomerToFiber(config['monomer_to_fiber'])
        alternator = Alternator(config['alternator'])

        return {
            'periodic_event': periodic_event,
            'readdy_actin': readdy,
            'medyan': medyan,
            'fiber_to_monomer': fiber_to_monomer,
            'monomer_to_fiber': monomer_to_fiber,
            'alternator': alternator,
        }

    def generate_topology(self, config):
        return {
            'periodic_event': {
                'event_trigger': ('alternate_trigger',),
                'period_index': ('period_index',)},
            'readdy_actin': {
                'monomers': ('monomers',)},
            'medyan': {
                'fibers': ('fibers',)},
            'fiber_to_monomer': {
                'fibers': ('fibers',),
                'monomers': ('monomers',)},
            'monomer_to_fiber': {
                'fibers': ('fibers',),
                'monomers': ('monomers',)},
            'alternator': {
                'alternate_trigger': ('alternate_trigger',),
                'choices': {
                    'medyan_active': ('choices', 'medyan_active',),
                    'readdy_active': ('choices', 'readdy_active',)}}}


def test_actin_fiber():
    initial_state = {
        'choices': {
            'medyan_active': True,
            'readdy_active': False},
        'fibers': {
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

    medyan_config = {
        'medyan_executable': '/home/youdonotexist/Downloads/medyan-4.2.0/build/medyan',
        'transform_points': [500, 500, 500]}

    actin_fiber_config = {
        'medyan': medyan_config}
    actin_fiber = ActinFiber(actin_fiber_config)

    composite = actin_fiber.generate()
    composite['initial_state'] = initial_state

    engine = Engine(
        processes=composite['processes'],
        topology=composite['topology'],
        initial_state=composite['initial_state'],
        emitter='simularium',
        emit_processes=True)

    engine.update(15)
    engine.emitter.get_data()


if __name__ == '__main__':
    test_actin_fiber()
