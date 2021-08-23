#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Creates a .xyz file for each .out file in the current directory."""

import os


def main():
    """Serve as main."""
    out_files = [f for f in os.listdir('.') if os.path.isfile(f) and
                 f.endswith('.out')]

    for out_file in out_files:
        with open(out_file, mode="r") as out_reader:
            num_atoms = None
            final_xyz_lines = list()
            current_xyz = list()

            in_coordinates_section = False
            for line in out_reader:
                if "  Coordinates:" in line:
                    in_coordinates_section = True
                elif "  SCF Energy:" in line or "Normal termination." in line:
                    if num_atoms is None:
                        num_atoms = str(len(current_xyz))
                    final_xyz_lines.append(num_atoms)
                    final_xyz_lines.append("")
                    final_xyz_lines.extend(current_xyz)
                    current_xyz = list()
                    in_coordinates_section = False
                elif in_coordinates_section and line.strip() != "":
                    current_xyz.append(line.strip())

        with open(out_file[:-4] + ".xyz", mode="w") as xyz_writer:
            for line in final_xyz_lines:
                xyz_writer.write(line + "\n")


if __name__ == "__main__":
    main()
