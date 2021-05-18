#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Serves as the main."""

import sys
import traceback

from milo_1_0_1 import electronic_structure_program_handler as esph
from milo_1_0_1 import force_propagation_handler as fph
from milo_1_0_1 import initial_energy_sampler
from milo_1_0_1 import input_parser
from milo_1_0_1 import program_state as ps


def main():
    """Run Milo."""
    try:
        print_header()

        program_state = ps.ProgramState()
        program_handler = esph.get_program_handler(program_state)
        propagation_handler = fph.get_propagation_handler(program_state)

        input_parser.parse_input(sys.stdin, program_state)

        if len(program_state.velocities) == 0:
            initial_energy_sampler.generate(program_state)

        print_trajectory_units_header()

        print("### Step 0: 0.0 fs ".ljust(66, '-'))
        print_structure(program_state)

        while not end_conditions_met(program_state):
            program_handler.generate_forces(program_state)
            propagation_handler.run_next_step(program_state)

            print_energy(program_state)
            print_forces(program_state)
            print_accelerations(program_state)
            print_velocities(program_state)

            # Change to next time step
            print()
            program_state.current_step += 1
            current_time = round(program_state.current_step *
                            program_state.step_size.as_femtosecond(), 10)
            print(f"### Step {program_state.current_step}: "
                  f"{current_time} fs ".ljust(66, '-'))

            print_structure(program_state)

        print_footer()

        if program_state.output_xyz_file:
            output_xyz_file(program_state)

    except Exception as e:
        print("\n\n")
        print("Oh no! It looks like there was an error! Error message:")
        print(e)
        print("\n\nPython error details:")
        print(traceback.format_exc())
        raise


def end_conditions_met(program_state):
    """Check conditions to see if the program should abort."""
    # Step count conditions
    if (program_state.max_steps is not None
            and program_state.current_step >= program_state.max_steps):
        return True
    return False


def output_xyz_file(program_state):
    """Write .xyz file."""
    with open(f"{program_state.job_name}.xyz", 'w') as file:
        for current_step, structure in enumerate(program_state.structures):
            file.write(f"{program_state.number_atoms}\n")
            current_time = round(current_step *
                            program_state.step_size.as_femtosecond(), 10)
            file.write(f"Step {current_step}: {current_time} fs\n")
            for atom, (x, y, z) in zip(program_state.atoms,
                                       structure.as_angstrom()):
                file.write(f"{atom.symbol} {x:15.6f} {y:15.6f} {z:15.6f}\n")


def print_structure(program_state):
    """Print the last structure list in program_state."""
    print("  Coordinates:")
    for atom, position in zip(program_state.atoms,
                              program_state.structures[-1].as_angstrom()):
        print(f"    {atom.symbol.ljust(2)} {position[0]:15.6f} "
              f"{position[1]:15.6f} {position[2]:15.6f}")


def print_velocities(program_state):
    """Print the last set of velocities in program_state."""
    print("  Velocities:")
    for atom, velocity in zip(program_state.atoms,
                              program_state.velocities[-1].as_meter_per_sec()):
        print(f"    {atom.symbol.ljust(2)} {velocity[0]:15.6e} "
              f"{velocity[1]:15.6e} {velocity[2]:15.6e}")


def print_accelerations(program_state):
    """Print the last set of accelerations in program_state."""
    print("  Accelerations:")
    for atom, acceleration in zip(program_state.atoms,
                                  program_state.accelerations[-1].
                                  as_meter_per_sec_sqrd()):
        print(f"    {atom.symbol.ljust(2)} {acceleration[0]:15.6e} "
              f"{acceleration[1]:15.6e} {acceleration[2]:15.6e}")


def print_forces(program_state):
    """Print the last set of forces in program_state."""
    print("  Forces:")
    for atom, force in zip(program_state.atoms,
                           program_state.forces[-1].as_newton()):
        print(f"    {atom.symbol.ljust(2)} {force[0]:15.6e} "
              f"{force[1]:15.6e} {force[2]:15.6e}")


def print_energy(program_state):
    """Print the last energy in program_state in Hartrees."""
    print("  SCF Energy:")
    print(f"    {program_state.energies[-1].as_hartree(0):.8f}")


def print_header():
    """Print the output file header."""
    print("Thank you for using")
    print()
    print("               ___   ___   ___   ___       _______ ")
    print("              |   |_|   | |   | |   |     |       |")
    print("              |         | |   | |   |     |   _   |")
    print("              |         | |   | |   |     |  | |  |")
    print("              |  ||_||  | |   | |   |___  |  |_|  |")
    print("              |  |   |  | |   | |       | |       |")
    print("              |__|   |__| |___| |_______| |_______|")
    print()
    print("Milo is free software; you can redistribute it and/or modify it")
    print("under the terms of the GNU General Public License as published by")
    print("the Free Software Foundation, either version 3 of the License, or")
    print("(at your option) any later version.")
    print()
    print("Milo is distributed in the hope that it will be useful,")
    print("but WITHOUT ANY WARRANTY; without even the implied warranty of")
    print("MERCHNATABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the")
    print("GNU General Public License for more details.")
    print()
    print("You should have received a copy of the GNU General Public License")
    print("along with Milo. If not, see <http://www.gnu.org/licenses/>.")
    print()
    print("Copyright 2020 Brigham Young University")
    print()
    print("You are using Milo 1.0.1 (18 May 2021 Update)")
    print()
    print("Authors:")
    print("Matthew S. Teynor, Nathan Wohlgemuth, Lily Carlson, Johnny Huang,")
    print("Samuel L. Pugh, Benjamin O. Grant, R. Spencer Hamilton, Ryan")
    print("Carlsen, Daniel H. Ess")
    print()
    print("Please cite Milo as:")
    print("Milo, Revision 1.0.1, M. S. Teynor, N. Wohlgemuth, L. Carlson,")
    print("J. Huang, S. L. Pugh, B. O. Grant, R. S. Hamilton, R. Carlsen,")
    print("D. H. Ess, Brigham Young University, Provo UT, 2021.")
    print()


def print_trajectory_units_header():
    """Print the units that the data from each trajectory step will use."""
    print("### Starting Trajectory ".ljust(66, '-'))
    print("  Units for trajectory output:")
    print("    Coordinates    angstrom")
    print("    SCF Energy     hartree")
    print("    Forces         newton")
    print("    Accelerations  meter/second^2")
    print("    Velocities     meter/second")
    print()


def print_footer():
    """Print the closing message to the output file."""
    print()
    print()
    print("Normal termination.")
    print("Thank you for using Milo!")


if __name__ == "__main__":
    main()
