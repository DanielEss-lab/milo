import unittest
from milo_1_0_3.atom import Atom


class AtomTests(unittest.TestCase):
    def test_atom_creation(self):
        atom = Atom("H", 1, 1, 1.00782503223)
        self.assertEqual(atom.symbol, "H")
        self.assertEqual(atom.atomic_number, 1)
        self.assertEqual(atom.mass_number, 1)
        self.assertEqual(atom.mass, 1.00782503223)

    def test_atom_from_symbol(self):
        atom = Atom.from_symbol("H")
        self.assertEqual(atom.symbol, "H")
        self.assertEqual(atom.atomic_number, 1)
        self.assertEqual(atom.mass_number, 1)
        self.assertEqual(atom.mass, 1.00782503223)

    def test_atom_from_atomic_number(self):
        atom = Atom.from_atomic_number(1)
        self.assertEqual(atom.symbol, "H")
        self.assertEqual(atom.atomic_number, 1)
        self.assertEqual(atom.mass_number, 1)
        self.assertEqual(atom.mass, 1.00782503223)

    def test_str_representation(self):
        atom = Atom("H", 1, 1, 1.00782503223)
        self.assertEqual(str(atom), "H    1.0078250 amu")

    def test_change_mass_exact(self):
        atom = Atom("H", 1, 1, 1.00782503223)
        atom.change_mass("1.008")
        self.assertEqual(atom.mass, 1.008)
        self.assertEqual(atom.mass_number, 1)

    def test_change_mass_integer(self):
        atom = Atom("H", 1, 1, 1.00782503223)
        atom.change_mass("2")
        self.assertEqual(atom.mass, 2.01410177812)
        self.assertEqual(atom.mass_number, 2)

    def test_change_mass_invalid(self):
        atom = Atom("H", 1, 1, 1.00782503223)
        with self.assertRaises(ValueError):
            atom.change_mass("abc")
        self.assertEqual(atom.mass, 1.00782503223)
        self.assertEqual(atom.mass_number, 1)


if __name__ == "__main__":
    unittest.main()
