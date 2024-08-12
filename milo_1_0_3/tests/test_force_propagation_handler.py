import unittest
from milo_1_0_3 import enumerations as enums
from milo_1_0_3.force_propagation_handler import (
    get_propagation_handler,
    Verlet,
    VelocityVerlet,
)


class TestForcePropagationHandler(unittest.TestCase):
    def test_get_propagation_handler_verlet(self):
        program_state = MockProgramState(
            propagation_algorithm=enums.PropagationAlgorithm.VERLET
        )
        handler = get_propagation_handler(program_state)
        self.assertEqual(handler, Verlet)

    def test_get_propagation_handler_velocity_verlet(self):
        program_state = MockProgramState(
            propagation_algorithm=enums.PropagationAlgorithm.VELOCITY_VERLET
        )
        handler = get_propagation_handler(program_state)
        self.assertEqual(handler, VelocityVerlet)

    def test_get_propagation_handler_unknown_algorithm(self):
        program_state = MockProgramState(propagation_algorithm="unknown")
        with self.assertRaises(ValueError):
            get_propagation_handler(program_state)


class MockProgramState:
    def __init__(self, propagation_algorithm):
        self.propagation_algorithm = propagation_algorithm


if __name__ == "__main__":
    unittest.main()
