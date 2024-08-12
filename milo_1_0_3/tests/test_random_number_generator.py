import unittest
from milo_1_0_3.random_number_generator import RandomNumberGenerator


class RandomNumberGeneratorTests(unittest.TestCase):
    def setUp(self):
        self.rng = RandomNumberGenerator()

    def test_uniform(self):
        result = self.rng.uniform()
        self.assertGreaterEqual(result, 0)
        self.assertLessEqual(result, 1)

    def test_edge_weighted(self):
        result = self.rng.edge_weighted()
        self.assertGreaterEqual(result, -1)
        self.assertLessEqual(result, 1)

    def test_gaussian(self):
        result = self.rng.gaussian()
        self.assertGreaterEqual(result, -1)
        self.assertLessEqual(result, 1)

    def test_one_or_neg_one(self):
        result = self.rng.one_or_neg_one()
        self.assertIn(result, [-1, 1])


if __name__ == "__main__":
    unittest.main()
