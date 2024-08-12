import unittest
from milo_1_0_3 import enumerations as en


class TestAccelerationUnits(unittest.TestCase):
    def test_meter_per_sec_squared(self):
        self.assertEqual(en.AccelerationUnits.METER_PER_SEC_SQRD.value, 1)


class TestAngleUnits(unittest.TestCase):
    def test_radians(self):
        self.assertEqual(en.AngleUnits.RADIANS.value, 1)

    def test_degrees(self):
        self.assertEqual(en.AngleUnits.DEGREES.value, 2)


class TestDistanceUnits(unittest.TestCase):
    def test_angstrom(self):
        self.assertEqual(en.DistanceUnits.ANGSTROM.value, 1)

    def test_bohr(self):
        self.assertEqual(en.DistanceUnits.BOHR.value, 2)

    def test_meter(self):
        self.assertEqual(en.DistanceUnits.METER.value, 3)


class TestEnergyUnits(unittest.TestCase):
    def test_joule(self):
        self.assertEqual(en.EnergyUnits.JOULE.value, 1)

    def test_kcal_per_mole(self):
        self.assertEqual(en.EnergyUnits.KCAL_PER_MOLE.value, 2)

    def test_millidyne_angstrom(self):
        self.assertEqual(en.EnergyUnits.MILLIDYNE_ANGSTROM.value, 3)

    def test_hartree(self):
        self.assertEqual(en.EnergyUnits.HARTREE.value, 4)


class TestGeometryDisplacement(unittest.TestCase):
    def test_none(self):
        self.assertEqual(en.GeometryDisplacement.NONE.value, 1)

    def test_edge_weighted(self):
        self.assertEqual(en.GeometryDisplacement.EDGE_WEIGHTED.value, 2)

    def test_gaussian_distribution(self):
        self.assertEqual(en.GeometryDisplacement.GAUSSIAN_DISTRIBUTION.value, 3)

    def test_uniform(self):
        self.assertEqual(en.GeometryDisplacement.UNIFORM.value, 4)


class TestMassUnits(unittest.TestCase):
    def test_amu(self):
        self.assertEqual(en.MassUnits.AMU.value, 1)

    def test_kilogram(self):
        self.assertEqual(en.MassUnits.KILOGRAM.value, 2)

    def test_gram(self):
        self.assertEqual(en.MassUnits.GRAM.value, 3)


class TestOscillatorType(unittest.TestCase):
    def test_quasiclassical(self):
        self.assertEqual(en.OscillatorType.QUASICLASSICAL.value, 1)

    def test_classical(self):
        self.assertEqual(en.OscillatorType.CLASSICAL.value, 2)


class TestPhaseDirection(unittest.TestCase):
    def test_bring_together(self):
        self.assertEqual(en.PhaseDirection.BRING_TOGETHER.value, 1)

    def test_push_apart(self):
        self.assertEqual(en.PhaseDirection.PUSH_APART.value, 2)

    def test_random(self):
        self.assertEqual(en.PhaseDirection.RANDOM.value, 3)


class TestProgramID(unittest.TestCase):
    def test_gaussian_16(self):
        self.assertEqual(en.ProgramID.GAUSSIAN_16.value, 1)

    def test_g16(self):
        self.assertEqual(en.ProgramID.G16.value, 1)

    def test_gaussian_09(self):
        self.assertEqual(en.ProgramID.GAUSSIAN_09.value, 2)

    def test_g09(self):
        self.assertEqual(en.ProgramID.G09.value, 2)


class TestPropagationAlgorithm(unittest.TestCase):
    def test_verlet(self):
        self.assertEqual(en.PropagationAlgorithm.VERLET.value, 1)

    def test_velocity_verlet(self):
        self.assertEqual(en.PropagationAlgorithm.VELOCITY_VERLET.value, 2)


class TestRotationalEnergy(unittest.TestCase):
    def test_yes(self):
        self.assertEqual(en.RotationalEnergy.YES.value, 1)

    def test_no(self):
        self.assertEqual(en.RotationalEnergy.NO.value, 2)


class TestTimeUnits(unittest.TestCase):
    def test_second(self):
        self.assertEqual(en.TimeUnits.SECOND.value, 1)

    def test_femtosecond(self):
        self.assertEqual(en.TimeUnits.FEMTOSECOND.value, 2)


class TestVelocityUnits(unittest.TestCase):
    def test_meter_per_sec(self):
        self.assertEqual(en.VelocityUnits.METER_PER_SEC.value, 1)

    def test_angstrom_per_fs(self):
        self.assertEqual(en.VelocityUnits.ANGSTROM_PER_FS.value, 2)

    def test_angstrom_per_sec(self):
        self.assertEqual(en.VelocityUnits.ANGSTROM_PER_SEC.value, 3)
