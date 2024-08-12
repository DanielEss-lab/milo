import os
import unittest
from milo_1_0_3.tools import parse_xyz as px


class TestParseXYZ(unittest.TestCase):
    def test_main(self):
        # Create temporary test files
        with open("test1.out", mode="w") as f:
            f.write("  Coordinates:\n")
            f.write("  1.0  2.0  3.0\n")
            f.write("  4.0  5.0  6.0\n")
            f.write("  SCF Energy:\n")
            f.write("  Normal termination.\n")

        with open("test2.out", mode="w") as f:
            f.write("  Coordinates:\n")
            f.write("  7.0  8.0  9.0\n")
            f.write("  10.0  11.0  12.0\n")
            f.write("  SCF Energy:\n")
            f.write("  Normal termination.\n")

        # Run the main function
        px.main()

        # Check if .xyz files are created
        self.assertTrue(os.path.isfile("test1.xyz"))
        self.assertTrue(os.path.isfile("test2.xyz"))

        # Clean up temporary test files
        os.remove("test1.out")
        os.remove("test2.out")
        os.remove("test1.xyz")
        os.remove("test2.xyz")


if __name__ == "__main__":
    unittest.main()
