#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create an input files that will continue a job from the last completed step.

The script will also create a SLURM submission script, unless turned off with
the --no_script flag.
"""

import argparse
import os
import sys

from milo_1_0_3 import input_parser
from milo_1_0_3 import program_state as ps


def main():
    """Parse frequency file and print to new Milo input."""
    parser = argparse.ArgumentParser(description="Make a Milo input file "
                                     "to restart a job from the last completed"
                                     " step.")
    parser.add_argument('job_to_restart', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin,
                        help="Milo output file. <stdin> by default.")
    parser.add_argument('new_input_filename', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="New Milo input file. <stdout> by default.")
    parser.add_argument("--no_script", action="store_true",
                        help="Stops SLURM submission scripts from being made.")
    args = parser.parse_args()

    # Get ### Input File section from output file
    in_input_section = False
    old_input = list()
    for line in args.job_to_restart:
        if "### Input File" in line:
            in_input_section = True
        elif in_input_section:
            if "### Default Parameters Being Used" in line:
                break
            old_input.append(line)

    # Parse the old input into program_state, w/o printing input_parser output
    program_state = ps.ProgramState()
    stdout = sys.stdout
    null_output = open(os.devnull, 'w')
    sys.stdout = null_output
    input_parser.parse_input(old_input, program_state)
    sys.stdout = stdout
    null_output.close()

    # Get random seed
    for line in args.job_to_restart:
        if "### Random Seed" in line:
            break
    random_seed = next(args.job_to_restart, "").strip()

    # Jump to first completed geometry
    for line in args.job_to_restart:
        if "###" in line and "Step" in line:
            break

    # Continue through input file to get most recent step
    current_coordinates = list()
    completed_coordinates = list()
    current_velocities = list()
    completed_velocities = list()
    completed_step_number = 0

    for line in args.job_to_restart:
        if "###" in line and "Step" in line:
            completed_step_number = int(line.split(":")[0].split()[-1]) - 1
            completed_coordinates = current_coordinates
            current_coordinates = list()
            completed_velocities = current_velocities
            current_velocities = list()
        elif "Coordinates:" in line:
            for _ in range(program_state.number_atoms):
                current_coordinates.append(next(args.job_to_restart, None))
        elif "Velocities:" in line:
            for _ in range(program_state.number_atoms):
                current_velocities.append(next(args.job_to_restart, None))

    # Remove atom labels from velocities section
    completed_velocities = [line.split(maxsplit=1)[1].rjust(50)
                            for line in completed_velocities]

    # Output
    print_section(args.new_input_filename, "comment",
        get_output_comment(args.job_to_restart, completed_step_number))
    print_section(args.new_input_filename, "job", get_job_section(old_input,
        completed_step_number, random_seed))
    print_section(args.new_input_filename, "molecule",
        "".join([f"    {program_state.charge} {program_state.spin}\n"]
        + completed_coordinates))
    print_section(args.new_input_filename, "isotope",
        get_isotope_section(program_state))
    print_section(args.new_input_filename, "velocities",
        "".join(completed_velocities))
    if program_state.gaussian_footer is not None:
        print_section(args.new_input_filename, "gaussian_footer",
            program_state.gaussian_footer)


def get_job_section(old_input, current_step, random_seed):
    """Create job section."""
    section = list()
    section.append(f"    current_step            {current_step}\n")
    in_job_section = False
    for line in old_input:
        if "$job" in line:
            in_job_section = True
        elif in_job_section and "$end" in line:
            in_job_section = False
        elif in_job_section:
            if "random_seed" in line:
                continue
            if "current_step" in line:
                continue
            section.append(line)
    section.append(f"    random_seed             {random_seed}\n")
    return "".join(section)


def get_isotope_section(program_state):
    """Create isotope section."""
    section = list()
    for i in range(program_state.number_atoms):
        section.append(f"    {(i + 1):< 3d} "
                       f"{program_state.atoms[i].mass:10.5f}\n")
    return "".join(section)


def get_output_comment(input_iterable, current_step):
    """Return comment with parsed output file name and date of parsing."""
    from datetime import datetime
    import os

    line = [f"    Input file to restart Milo job from step {current_step} of "]
    if input_iterable != sys.stdin:
        line.append(os.path.basename(input_iterable.name))
    else:
        try:
            name = os.readlink('/proc/self/fd/0').split('/')[-1].split('.')[0]
            line.append(name)
        except FileNotFoundError:
            line.append("unknown_job")
    line.append(datetime.now().strftime(". Generated on %d-%b-%Y at %X."))
    return "".join(line)


def print_section(output_iterable, section_name, inside):
    """Print a section to output_iterable."""
    stdout = sys.stdout
    sys.stdout = output_iterable

    print(f"${section_name}")
    print(inside.rstrip())
    print("$end")
    print()

    sys.stdout = stdout


if __name__ == "__main__":
    main()
