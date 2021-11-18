#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a Milo input file from a frequency calculation.

It must be a Gaussian 09 or 16 high-precision frequency calculation. You
request this with '# freq=(hpmodes) ... '.
"""

import argparse
import sys

from milo_1_0_3 import atom
from milo_1_0_3 import containers
from milo_1_0_3 import enumerations as enums
from milo_1_0_3 import exceptions
from milo_1_0_3 import program_state as ps


def main():
    """Parse frequency file and print to new Milo input."""
    parser = argparse.ArgumentParser(description="Make a Milo input file "
                                     "from a high-precision Gaussian frequency"
                                     " calculation.\n")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin, help="Frequency calculation file. "
                        "<stdin> by default.")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout, help="New Milo input file. "
                        "<stdout> by default.")
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="Print other parameters in $job section. "
                        "-v for common parameters, -vv for all parameters")
    args = parser.parse_args()
    program_state = ps.ProgramState()
    try:
        parse_gaussian_header(args.infile, program_state)
        parse_gaussian_charge_spin(args.infile, program_state)
        parse_gaussian_molecule_data(args.infile, program_state)
        parse_gaussian_frequency_data(args.infile, program_state)
        parse_gaussian_isotope_data(args.infile, program_state)
        print_job_section(args.outfile, program_state, args.verbose)
        print_output_comment(args.infile, args.outfile)
        print_molecule_section(args.outfile, program_state)
        print_frequency_data_section(args.outfile, program_state)
    except Exception as e:
        print("Oh no! It looks like there was an error!")
        print("Error message:", e)
        print("\nPython error details:")
        raise


def parse_gaussian_header(input_iterable, program_state):
    """
    Parse gaussian_header from frequency file.

    Looking for:
        ******************************************
        ------------------------------
        # opt freq=hpmodes m062x/3-21g
        ------------------------------
    Result:
        gaussian_header = 'm062x/3-21g'
    """
    past_warning = False
    lines = list()
    for line in input_iterable:
        if "*****" in line:
            past_warning = True
        if past_warning and "-----" in line:
            for next_line in input_iterable:
                if "-----" in next_line:
                    break
                lines.append(next_line[1:].strip("\n"))
            clean_line = "".join(lines).strip()
            if "hpmodes" not in clean_line.casefold():
                raise exceptions.InputError("Must be high-precision frequency "
                                            "calculation. Use 'freq=hpmodes'.")
            tokens = clean_line.split()
            tokens = [x for x in tokens if "#" not in x
                      and "opt" not in x.casefold()
                      and "freq" not in x.casefold()]
            program_state.gaussian_header = " ".join(tokens)
            return
    raise exceptions.InputError("Error parsing gaussian_header.")


def parse_gaussian_charge_spin(input_iterable, program_state):
    """
    Parse charge and spin multiplicity from frequency file.

    Looking for:
         ---------------------------------------------
         Symbolic Z-matrix:
         Charge =  0 Multiplicity = 1
         O                    -0.19334  -0.19871   0.
    """
    for line in input_iterable:
        if "Charge =" in line:
            program_state.charge = int(line.split()[2])
            program_state.spin = int(line.split()[5])
            return
    raise exceptions.InputError("Error parsing charge and spin multiplicity.")


def parse_gaussian_molecule_data(input_iterable, program_state):
    """
    Parse molecule data from frequency file.

    Will pull the last "Standard orientation:" in the log file, or the last
    "Input orientation:" if there is no "Standard orientation:" (for example,
    if the nosymm keyword is used).

    Looking for:
                                 Standard orientation:
         ---------------------------------------------------------------------
         Center     Atomic      Atomic             Coordinates (Angstroms)
         Number     Number       Type             X           Y           Z
         ---------------------------------------------------------------------
    """
    for line in input_iterable:
        if "Harmonic frequencies (cm**-1)" in line:
            return
        if "Input orientation:" in line or "Standard orientation:" in line:
            positions = containers.Positions()
            for coordinate_line in input_iterable:
                if ("Rotational constants" in coordinate_line or
                        "Distance matrix" in coordinate_line):
                    break
                coordinates = coordinate_line.split()
                if coordinates[0].isnumeric():
                    x = float(coordinates[3])
                    y = float(coordinates[4])
                    z = float(coordinates[5])
                    positions.append(x, y, z, enums.DistanceUnits.ANGSTROM)
            program_state.input_structure = positions
    raise exceptions.InputError("Error parsing molecule data.")


def parse_gaussian_frequency_data(input_iterable, program_state):
    """
    Parse frequency data from frequency file.

    Will pull the first time they are listed (with high-precision).

    Looking for:
               Frequencies ---  1682.1354 3524.4296 3668.7401
            Reduced masses ---     1.0895    1.0389    1.0827
           Force constants ---     1.8163    7.6032    8.5864
            IR Intensities ---    52.8486    4.2243    0.3831
         Coord Atom Element:
           1     1     8         -0.00000   0.00000  -0.00000
           2     1     8          0.00000  -0.00000  -0.07070
           3     1     8         -0.07382   0.04553  -0.00000
           1     2     1          0.00000   0.00000   0.00000
           2     2     1          0.39258   0.60700   0.56106
           3     2     1          0.58580  -0.36126  -0.42745
           1     3     1          0.00000  -0.00000   0.00000
           2     3     1         -0.39258  -0.60700   0.56106
           3     3     1          0.58580  -0.36126   0.42745
           Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scatt
           activities (A**4/AMU), depolarization ratios for plane and unpolariz
    """
    has_started = False
    for line in input_iterable:
        if "Frequencies ---" in line:
            has_started = True
            for frequency in line.split()[2:]:
                program_state.frequencies.append(float(frequency),
                                                 enums.FrequencyUnits
                                                 .RECIP_CM)
        elif "Reduced masses ---" in line:
            for reduced_mass in line.split()[3:]:
                program_state.reduced_masses\
                    .append(float(reduced_mass), enums.MassUnits.AMU)
        elif "Force constants ---" in line:
            for force_constant in line.split()[3:]:
                program_state.force_constants\
                    .append(float(force_constant), enums.ForceConstantUnits
                            .MILLIDYNE_PER_ANGSTROM)
        elif "Coord Atom Element:" in line:
            data_in_columns = list()
            for coordinate_line in input_iterable:
                if ("Harmonic frequencies (cm**-1)" in coordinate_line
                        or "                    " in coordinate_line):
                    break
                data_in_columns.append(coordinate_line.split()[3:])
            data_in_rows = list(zip(*data_in_columns))
            for frequency in data_in_rows:
                program_state.mode_displacements.append(containers.Positions())
                for x, y, z in zip(*[iter(frequency)] * 3):
                    program_state.mode_displacements[-1].append(float(x),
                        float(y), float(z), enums.DistanceUnits.ANGSTROM)
        elif has_started and "activities (A**4/AMU)" in line:
            return
    raise exceptions.InputError("Error parsing frequency data.")


def parse_gaussian_isotope_data(input_iterable, program_state):
    """
    Parse isotope and atomic number data from frequency file.

    Looking for:
         -------------------
         - Thermochemistry -
         -------------------
         Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
         Atom     1 has atomic number  8 and mass  15.99491
         Atom     2 has atomic number  1 and mass   1.00783
         Atom     3 has atomic number  1 and mass   1.00783
         Molecular mass:    18.01056 amu.
    """
    for line in input_iterable:
        if "Thermochemistry" in line:
            atoms = list()
            for mass_line in input_iterable:
                if "Molecular mass" in mass_line:
                    break
                split_line = mass_line.split()
                if split_line[0] == "Atom":
                    atomic_number = int(split_line[5])
                    atoms.append(atom.Atom.from_atomic_number(atomic_number))
                    atoms[-1].change_mass(split_line[8])
            program_state.atoms = atoms
            return
    raise exceptions.InputError("Error parsing isotope data.")


def print_section(output_iterable, section_name, inside):
    """Print a section to output_iterable."""
    stdout = sys.stdout
    sys.stdout = output_iterable

    print(f"${section_name}")
    print(inside)
    print("$end")
    print()

    sys.stdout = stdout


def print_job_section(output_iterable, program_state, verbose):
    """
    Print the $job section with gaussian_header from program_state.

    verbose controls how other job parameters are printed.
    """
    section = list()
    section.append("    gaussian_header         "
                   f"{program_state.gaussian_header}")
    if verbose >= 1:
        section.append("    # step_size               1.00  # in femtoseconds")
        section.append("    # max_steps               100  # or no_limit")
        section.append("    # temperature             298.15  # in kelvin")
        section.append("    # phase                   bring_together n m"
                       "  #  or  push_apart n m")
        section.append("    # memory                  24  # in GB")
        section.append("    # processors              24")
        section.append("    # random_seed             generate  # or an "
                       "integer")
    if verbose >= 2:
        section.append("    # oscillator_type         quasiclassical")
        section.append("    # geometry_displacement   off")
        section.append("    # rotational_energy       off")
        section.append("    # energy_boost            off")
        section.append("    # integration_algorithm   verlet")
        section.append("    # program                 gaussian16")
        section.append("    # fixed_mode_direction    n 1  # or n -1")
    print_section(output_iterable, "job", "\n".join(section))


def print_molecule_section(output_iterable, program_state):
    """Print $molecule section with data from program_state."""
    section = list()
    section.append(f"    {program_state.charge} {program_state.spin}")
    for _atom, (x, y, z) in zip(program_state.atoms,
                               program_state.input_structure.as_angstrom()):
        section.append(f"    {_atom.symbol} {x:12.6f} {y:12.6f} {z:12.6f}")
    print_section(output_iterable, "molecule", "\n".join(section))

    section = list()
    for i, _atom in enumerate(program_state.atoms, 1):
        section.append(f"    {i:< 3d} {_atom.mass:10.5f}")
    print_section(output_iterable, "isotope", "\n".join(section))


def print_frequency_data_section(output_iterable, program_state):
    """Print $frequencies section with data from program_state."""
    section = list()
    for frequency, reduced_mass, force_constant, mode_displacement in zip(
            program_state.frequencies.as_recip_cm(),
            program_state.reduced_masses.as_amu(),
            program_state.force_constants.as_millidyne_per_angstrom(),
            program_state.mode_displacements):
        section.append(f"   {frequency:10.4f} {reduced_mass:7.4f} "
                       f"{force_constant:7.4f}")
        for x, y, z in mode_displacement.as_angstrom():
            section.append(f"  {x:8.5f} {y:8.5f} {z:8.5f}")
        section.append("\n")
    section.pop()
    print_section(output_iterable, "frequency_data", "".join(section))


def print_output_comment(input_iterable, output_iterable):
    """Print comment with frequency file name and date of parsing."""
    from datetime import datetime
    import os

    comment = list()
    comment.append("    Frequency and molecule data parsed ")
    if input_iterable != sys.stdin:
        comment.append("from ")
        comment.append(os.path.basename(input_iterable.name))
        comment.append(" ")
    else:
        try:
            name = os.readlink('/proc/self/fd/0').split('/')[-1].split('.')[0]
            comment.append("from ")
            comment.append(name)
            comment.append(" ")
        except FileNotFoundError:
            comment.append("from <stdin> ")
    comment.append(datetime.now().strftime("on %d-%b-%Y at %X"))
    print_section(output_iterable, "comment", "".join(comment))


if __name__ == "__main__":
    main()
