from vivarium.core.composer import Composer
from vivarium.core.engine import Engine

from vivarium_models.processes.readdy_actin_process import ReaddyActinProcess
from vivarium_models.processes.fiber_controls import FiberControls

class ActinFiber(Composer):
    defaults = {}

    def __init__(self, config=None):
        super().__init__(config)

    def generate_processes(self, config):
        readdy = ReaddyActinProcess(config.get('readdy_actin', {}))
        fiber_controls = FiberControls(config.get('fiber_controls', {}))

        return {
            'readdy_actin': readdy,
            'fiber_controls': fiber_controls}

    def generate_topology(self, config):
        return {
            'readdy_actin': {
                'topologies': ('topologies',)},
            'fiber_controls': {
                'topologies': ('topologies',),
                'fibers': ('fibers',)}}

def test_actin_fiber():
    actin_fiber_config = {}
    actin_fiber = ActinFiber(actin_fiber_config)

    composite = actin_fiber.generate()
    composite['initial_state'] = actin_fiber.initial_state()
    engine = Engine(composite)

    engine.update(10)

    import ipdb; ipdb.set_trace()

if __name__ == '__main__':
    test_actin_fiber()
