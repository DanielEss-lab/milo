import unittest
from milo_1_0_3.containers import Positions, enums


class TestPositions(unittest.TestCase):
    def setUp(self):
        self.positions = Positions()
        self.positions.append(1, 2, 3, enums.DistanceUnits.ANGSTROM)
        self.positions.append(4, 5, 6, enums.DistanceUnits.ANGSTROM)

    def test_length(self):
        self.assertEqual(len(self.positions), 2)

    def test_alter_position(self):
        self.positions.alter_position(0, 10, 20, 30, enums.DistanceUnits.ANGSTROM)
        self.assertEqual(self.positions.as_angstrom(0), (10, 20, 30))

    def test_append(self):
        self.positions.append(7, 8, 9, enums.DistanceUnits.ANGSTROM)
        self.assertEqual(len(self.positions), 3)
        self.assertEqual(self.positions.as_angstrom(2), (7, 8, 9))

    def test_as_angstrom(self):
        self.assertEqual(self.positions.as_angstrom(), [(1, 2, 3), (4, 5, 6)])
        self.assertEqual(self.positions.as_angstrom(0), (1, 2, 3))

    def test_as_bohr(self):
        self.assertEqual(
            self.positions.as_bohr(),
            [
                (1.8897261246229133, 3.7794522492458267, 5.66917837386874),
                (7.558904498491653, 9.448630623114568, 11.33835674773748),
            ],
        )
        self.assertEqual(
            self.positions.as_bohr(0),
            (1.8897261246229133, 3.7794522492458267, 5.66917837386874),
        )
        self.assertEqual(
            self.positions.as_bohr(1),
            (7.558904498491653, 9.448630623114568, 11.33835674773748),
        )

    def test_as_meter(self):
        self.assertEqual(
            self.positions.as_meter(), [(1e-10, 2e-10, 3e-10), (4e-10, 5e-10, 6e-10)]
        )
        self.assertEqual(self.positions.as_meter(0), (1e-10, 2e-10, 3e-10))

    def test_add(self):
        other = Positions()
        other.append(2, 3, 4, enums.DistanceUnits.ANGSTROM)
        other.append(1, 1, 1, enums.DistanceUnits.ANGSTROM)
        result = self.positions + other
        self.assertEqual(len(result), 2)
        self.assertEqual(result.as_angstrom(0), (3, 5, 7))
        self.assertEqual(result.as_angstrom(1), (5, 6, 7))

    def test_sub(self):
        other = Positions()
        other.append(1, 1, 1, enums.DistanceUnits.ANGSTROM)
        other.append(2, 3, 4, enums.DistanceUnits.ANGSTROM)
        result = self.positions - other
        self.assertEqual(len(result), 2)
        self.assertEqual(result.as_angstrom(0), (0, 1, 2))
        self.assertEqual(result.as_angstrom(1), (2, 2, 2))

    def test_mul(self):
        result = self.positions * 2
        self.assertEqual(len(result), 2)
        self.assertEqual(result.as_angstrom(0), (2, 4, 6))
        self.assertEqual(result.as_angstrom(1), (8, 10, 12))

    def test_str(self):
        expected = "  1.000000   2.000000   3.000000\n  4.000000   5.000000   6.000000"
        self.assertEqual(str(self.positions), expected)

    def test_repr(self):
        expected = "<Positions object with 2 atoms.>"
        self.assertEqual(repr(self.positions), expected)


if __name__ == "__main__":
    unittest.main()
