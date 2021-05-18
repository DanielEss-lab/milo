#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parses input file and populates a ProgramState object."""

from copy import deepcopy
import os

from milo_1_0_1 import atom
from milo_1_0_1 import containers
from milo_1_0_1 import enumerations as enums
from milo_1_0_1 import exceptions

required_sections = ("$job", "$molecule", )

no_duplicate_sections = ("$molecule", "$isotope", "$velocities",
                         "$frequency_data", "$gaussian_footer")

mutually_exclusive_sections = (
    ("$velocities", "$frequency_data"),
)

required_job_parameters = ("gaussian_header", )

allowed_duplicate_parameters = ("fixed_mode_direction",
                                "fixed_vibrational_quanta")

mutually_exclusive_job_parameters = (
    ("gaussian_header", "qchem_options"),  # as an example
)

# This is only used to print what defaults are used, not to acutally set the
# values. The key needs to be the name of parameter used in the input file. For
# example, use 'program' not 'program_id'. The value is just a text descriptor
# of the default value.
parameters_with_defaults = {
    "max_steps": "no_limit",
    "phase": "random",
    "program": "gaussian16",
    "integration_algorithm": "verlet",
    "step_size": "1.00 fs",
    "temperature": "298.15 K",
    "energy_boost": "off",
    "oscillator_type": "quasiclassical",
    "geometry_displacement": "off",
    "rotational_energy": "off",
}


def parse_input(input_iterable, program_state):
    """Populate a ProgramState object from an input file."""
    # Break input into lines
    input = input_iterable.readlines()

    # Print entire input file
    print("### Input File ---------------------------------------------------")
    print("".join(input))
    print()

    # TOKENIZE
    # Example - input: ["$job", "  parameter a b c # comment", "$end"]
    #         - tokenized_lines: [["$job"], ["parameter", "a b c"], ["$end"]]
    # Remove in-line comments
    tokenized_lines = [line.split('#', maxsplit=1)[0] for line in input]
    # Strip whitespace
    tokenized_lines = [line.strip() for line in tokenized_lines]
    # Removes blank lines
    tokenized_lines = [line for line in tokenized_lines if line]
    # Splits line into first word and rest
    tokenized_lines = [line.split(maxsplit=1) for line in tokenized_lines]

    # Check for required sections
    token_keys = [tokens[0].casefold() for tokens in tokenized_lines]
    for section in required_sections:
        if section not in token_keys:
            raise exceptions.InputError(f"Could not find {section} section.")

    # Check against duplicate sections
    for section in no_duplicate_sections:
        if token_keys.count(section) > 1:
            raise exceptions.InputError(f"Multiple {section} sections.")

    # Check against mutually exclusive sections
    for m_e_tuple in mutually_exclusive_sections:
        has_m_e_section = [(section in token_keys) for section in m_e_tuple]
        if has_m_e_section.count(True) > 1:
            raise exceptions.InputError(f"{str(m_e_tuple)[1:-1]} are mutually "
                                        "exclusive.")

    # Break into different sections and remove block comments
    tokenized_lines_it = iter(tokenized_lines)

    job_tokens = list()
    molecule_tokens = list()
    isotope_tokens = list()
    velocities_tokens = list()
    frequency_data_tokens = list()

    for tokens in tokenized_lines_it:
        if tokens[0].casefold() == '$comment':
            while tokens[0].casefold() != '$end' and tokens is not None:
                tokens = next(tokenized_lines_it, None)

        elif tokens[0].casefold() == '$job':
            tokens = next(tokenized_lines_it, None)  # Advance past '$job'
            while tokens[0].casefold() != '$end' and tokens is not None:
                job_tokens.append(tokens)
                tokens = next(tokenized_lines_it, None)

        elif tokens[0].casefold() == '$molecule':
            tokens = next(tokenized_lines_it, None)  # Advance past '$molecule'
            while tokens[0].casefold() != '$end' and tokens is not None:
                molecule_tokens.append(tokens)
                tokens = next(tokenized_lines_it, None)

        elif tokens[0].casefold() == '$isotope':
            tokens = next(tokenized_lines_it, None)  # Advance past '$isotope'
            while tokens[0].casefold() != '$end' and tokens is not None:
                isotope_tokens.append(tokens)
                tokens = next(tokenized_lines_it, None)

        elif tokens[0].casefold() == '$velocities':
            tokens = next(tokenized_lines_it, None)
            while tokens[0].casefold() != '$end' and tokens is not None:
                velocities_tokens.append(tokens)
                tokens = next(tokenized_lines_it, None)

        elif tokens[0].casefold() == '$frequency_data':
            tokens = next(tokenized_lines_it, None)
            while tokens[0].casefold() != '$end' and tokens is not None:
                frequency_data_tokens.append(tokens)
                tokens = next(tokenized_lines_it, None)

        elif "$" in tokens[0].casefold() and ("$gaussian_footer" !=
                tokens[0].casefold() and "$end" != tokens[0].casefold()):
            raise exceptions.InputError(f"Could not interpret {tokens[0]} "
                                        "section.")

    # CHECK LEGALITY OF INPUT
    # Check for required job parameters
    job_parameters = [tokens[0].casefold() for tokens in job_tokens]
    for parameter in required_job_parameters:
        if parameter not in job_parameters:
            raise exceptions.InputError(f"Could not find the required "
                                        f"{parameter} parameter in the $job "
                                        "section.")

    # Check against mutually exclusive job parameters
    for m_e_tuple in mutually_exclusive_job_parameters:
        has_m_e_parameter = [parameter in job_parameters
                             for parameter in m_e_tuple]
        if has_m_e_parameter.count(True) > 1:
            raise exceptions.InputError(f"{str(m_e_tuple)[1:-1]} are mutually "
                                        "exclusive.")

    # Check against duplicate parameters
    for parameter in job_parameters:
        if (parameter not in allowed_duplicate_parameters and
                job_parameters.count(parameter) > 1):
            raise exceptions.InputError(f"The '{parameter}' parameter can "
                                        "only be listed once.")

    # POPULATE DATA
    # Populate program_state with data from $molecule and $isotope
    try:
        charge_and_spin = molecule_tokens.pop(0)
        program_state.charge = int(charge_and_spin[0])
        program_state.spin = int(charge_and_spin[1])
    except (IndexError, ValueError):
        raise exceptions.InputError("Could not find charge and/or spin "
                                    "multiplicity in the $molecule section.")
    program_state.number_atoms = len(molecule_tokens)
    for atom_token in molecule_tokens:
        try:
            program_state.atoms.append(atom.Atom.from_symbol(atom_token[0]))
            x, y, z = atom_token[1].split()
            program_state.input_structure.append(float(x), float(y), float(z),
                                                 enums.DistanceUnits.ANGSTROM)
        except (IndexError, KeyError, ValueError):
            raise exceptions.InputError("Could not interpret "
                                        f"'{'  '.join(atom_token)}' in the "
                                        "$molecule section.")
    for mass_token in isotope_tokens:
        try:
            index = int(mass_token[0]) - 1  # '- 1' to bring to 0-based index
            program_state.atoms[index].change_mass(mass_token[1])
        except (IndexError, KeyError):
            raise exceptions.InputError("Could not interpret "
                                        f"'{'  '.join(mass_token)}' in the "
                                        "$isotope section.")
    program_state.structures.append(deepcopy(program_state.input_structure))

    # Populate program_state with job parameters
    for tokens in job_tokens:
        parameter = tokens[0].casefold()
        if parameter in parameters_with_defaults:
            del parameters_with_defaults[parameter]
        try:
            job_function = getattr(JobSection, parameter)
        except AttributeError:
            raise exceptions.InputError(f"Invalid parameter '{tokens[0]}' in "
                                        "$job section.")
        options = tokens[1] if len(tokens) > 1 else ""
        job_function(options, program_state)

    # Populate gaussian_footer from $gaussian_footer section
    # This is done from the raw input, so comments won't be taken out
    if "$gaussian_footer" in token_keys:
        in_section = False
        for i, line in enumerate(input):
            if "$gaussian_footer" in line:
                in_section = True
                start_index = i + 1
            elif in_section and "$end" in line:
                end_index = i
                break
        program_state.gaussian_footer = "".join(input[start_index:end_index])

    # Populate program_state with frequency data
    try:
        for frequency_token in frequency_data_tokens:
            program_state.frequencies.append(float(frequency_token[0]),
                                             enums.FrequencyUnits.RECIP_CM)
            data = frequency_token[1].split()
            reduced_mass = float(data.pop(0))
            program_state.reduced_masses.append(reduced_mass,
                                                enums.MassUnits.AMU)
            force_constant = float(data.pop(0))
            program_state.force_constants.append(force_constant,
                            enums.ForceConstantUnits.MILLIDYNE_PER_ANGSTROM)
            current_displacements = containers.Positions()
            for _ in range(program_state.number_atoms):
                x = float(data.pop(0))
                y = float(data.pop(0))
                z = float(data.pop(0))
                current_displacements.append(x, y, z,
                                             enums.DistanceUnits.ANGSTROM)
            program_state.mode_displacements.append(current_displacements)
    except (IndexError, ValueError):
        raise exceptions.InputError("Could not interpret $frequency_data"
                                    "section.")

    # Populate program_state with velocity data
    if "$velocities" in token_keys:
        program_state.velocities.append(containers.Velocities())
        try:
            for velocities_token in velocities_tokens:
                x = float(velocities_token[0])
                y, z = [float(component) for component in
                    velocities_token[1].split()]
                program_state.velocities[0].append(x, y, z,
                                             enums.VelocityUnits.METER_PER_SEC)
        except (IndexError, ValueError):
            raise exceptions.InputError("Could not interpret $velocities "
                                        "section.")
        if len(velocities_tokens) != program_state.number_atoms:
            raise exceptions.InputError("Number of atoms in $velocities and "
                                        "$molecule sections does not match.")

    # Pull job name from output file name
    try:
        name = os.readlink('/proc/self/fd/1').split('/')[-1].split('.')[0]
    except FileNotFoundError:
        name = "MiloJob"
    program_state.job_name = name

    # Print results from parsing
    print("### Default Parameters Being Used --------------------------------")
    for parameter in parameters_with_defaults:
        print("  ", parameter, ": ", parameters_with_defaults[parameter],
              sep="")
    if len(parameters_with_defaults) == 0:
        print("  (No defaults used.)")
    print()
    print("### Random Seed --------------------------------------------------")
    print(f"  {program_state.random.seed}")
    print()
    print("### Atomic Mass Data ---------------------------------------------")
    for i, atom_ in enumerate(program_state.atoms, 1):
        print(f"  {i:< 3d}  {atom_}")
    print()


class JobSection():
    """
    Namespace for all the parameter functions in $job.

    The name of each function must match the name of the input parameter. These
    functions are automatically called by:
            job_function = getattr(JobSection, tokens[0].casefold())
            options = tokens[1] if len(tokens) > 1 else ""
            job_function(options, program_state)

    Function parameters:
        options - a list containing the input split into the first word and
    """

    @staticmethod
    def energy_boost(options, program_state):
        """Populate program_state.energy_boost from options."""
        err_msg = (f"Could not interpret parameter 'energy_boost {options}'. "
                   "Expected 'energy_boost on min max' or 'energy_boost off'.")
        if options.casefold() == "off":
            program_state.energy_boost = enums.EnergyBoost.OFF
        elif options.split()[0].casefold() == "on":
            try:
                program_state.energy_boost = enums.EnergyBoost.ON
                min, max = float(options.split()[1]), float(options.split()[2])
                if min > max:
                    min, max = max, min
                program_state.energy_boost_min = min
                program_state.energy_boost_max = max
            except (IndexError, ValueError):
                raise exceptions.InputError(err_msg)
        else:
            raise exceptions.InputError(err_msg)

    @staticmethod
    def fixed_mode_direction(options, program_state):
        """Add mode direction (1 or -1) to fixed_mode_directions."""
        err_msg = ("Could not interpret parameter 'fixed_mode_direction "
                   f"{options}'. Expected 'fixed_mode_direction n 1', or "
                   "'fixed_mode_direction n -1', where n is the mode index.")
        try:
            mode, direction = [int(x) for x in options.split()]
        except ValueError:
            raise exceptions.InputError(err_msg)
        if mode < 1 or (direction != -1 and direction != 1):
            raise exceptions.InputError(err_msg)
        program_state.fixed_mode_directions[mode] = direction

    @staticmethod
    def fixed_vibrational_quanta(options, program_state):
        """Add vibrational quantum to fixed_vibrational_quanta."""
        err_msg = ("Could not interpret parameter 'fixed_vibrational_quanta "
                   f"{options}'. Expected 'fixed_vibrational_quanta n m', "
                   "where n is the mode index and m is the vibrational "
                   "quantum number (integer >= 0).")
        try:
            mode, quantum_number = [int(x) for x in options.split()]
        except ValueError:
            raise exceptions.InputError(err_msg)
        if mode < 1 or quantum_number < 0:
            raise exceptions.InputError(err_msg)
        program_state.fixed_vibrational_quanta[mode] = quantum_number

    @staticmethod
    def gaussian_header(options, program_state):
        """Populate program_state.gaussian_header from options."""
        program_state.gaussian_header = options

    @staticmethod
    def gaussian_footer(options, program_state):
        """Populate program_state.gaussian_header from options."""
        program_state.gaussian_footer = options.replace("\\n", "\n")

    @staticmethod
    def geometry_displacement(options, program_state):
        """Populate program_state.geometry_displacement_type from options."""
        err_msg = ("Could not interpret parameter 'geometry_displacement "
                   f"{options}'. Expected 'geometry_displacement edge_weighted"
                   "', 'geometry_displacement gaussian', 'geometry_displacemen"
                   "t uniform' or 'geometry_displacement off'.")
        if options.casefold() == "edge_weighted":
            program_state.geometry_displacement_type = \
                enums.GeometryDisplacement.EDGE_WEIGHTED
        elif options.casefold() == "gaussian":
            program_state.geometry_displacement_type = \
                enums.GeometryDisplacement.GAUSSIAN_DISTRIBUTION
        elif options.casefold() == "uniform":
            program_state.geometry_displacement_type = \
                enums.GeometryDisplacement.UNIFORM
        elif options.casefold() == "off":
            program_state.geometry_displacement_type = \
                (enums.GeometryDisplacement.NONE)
        else:
            raise exceptions.InputError(err_msg)

    @staticmethod
    def integration_algorithm(options, program_state):
        """Populate program_state.propagation_algorithm from options."""
        err_msg = (f"Could not interpret parameter 'integration_algorithm "
                   f"{options}'. Expected 'verlet' or 'velocity_verlet'.")
        if options.casefold() == "verlet":
            program_state.propagation_algorithm = \
                enums.PropagationAlgorithm.VERLET
        elif options.casefold() == "velocity_verlet":
            program_state.propagation_algorithm = \
                enums.PropagationAlgorithm.VELOCITY_VERLET
        else:
            raise exceptions.InputError(err_msg)

    @staticmethod
    def max_steps(options, program_state):
        """Populate program_state.max_steps from options."""
        err_msg = (f"Could not interpret parameter 'max_steps {options}'. "
                   "Expected 'max_steps integer' or 'no_limit'.")
        if options.casefold() == "no_limit":
            program_state.max_steps = None
        else:
            try:
                program_state.max_steps = int(options)
            except ValueError:
                raise exceptions.InputError(err_msg)

    @staticmethod
    def memory(options, program_state):
        """Populate program_state.memory_amount from options."""
        err_msg = (f"Could not interpret parameter 'memory {options}'. "
                   "Expected 'memory integer'.")
        try:
            program_state.memory_amount = int(options)
        except ValueError:
            raise exceptions.InputError(err_msg)

    @staticmethod
    def oscillator_type(options, program_state):
        """Populate program_state.oscillator_type from options."""
        err_msg = (f"Could not interpret parameter 'oscillator_type {options}'"
                   ". Expected 'oscillator_type classical' or 'oscillator_type"
                   " quasiclassical'.")
        if options.casefold() == "classical":
            program_state.oscillator_type = enums.OscillatorType.CLASSICAL
        elif options.casefold() == "quasiclassical":
            program_state.oscillator_type = enums.OscillatorType.QUASICLASSICAL
        else:
            raise exceptions.InputError(err_msg)

    @staticmethod
    def phase(options, program_state):
        """Populate program_state.phase_direction from options."""
        err_msg = (f"Could not interpret parameter 'phase {options}'. "
                   "Expected 'phase bring_together index1 index2', 'phase "
                   "push_apart index1 index2' or 'phase random'.")
        if options.casefold() == "random":
            program_state.phase_direction = enums.PhaseDirection.RANDOM
        elif options.split()[0].casefold() == "bring_together":
            program_state.phase_direction = enums.PhaseDirection.BRING_TOGETHER
            try:
                program_state.phase = (int(options.split()[1]),
                                       int(options.split()[2]))
            except (IndexError, ValueError):
                raise exceptions.InputError(err_msg)
        elif options.split()[0].casefold() == "push_apart":
            program_state.phase_direction = enums.PhaseDirection.PUSH_APART
            try:
                program_state.phase = (int(options.split()[1]),
                                       int(options.split()[2]))
            except (IndexError, ValueError):
                raise exceptions.InputError(err_msg)
        else:
            raise exceptions.InputError(err_msg)

    @staticmethod
    def processors(options, program_state):
        """Populate program_state.processor_count from options."""
        err_msg = (f"Could not interpret parameter 'processors {options}'. "
                   "Expected 'processors integer'.")
        try:
            program_state.processor_count = int(options)
        except ValueError:
            raise exceptions.InputError(err_msg)

    @staticmethod
    def program(options, program_state):
        """Populate program_state.program_id from options."""
        err_msg = (f"Could not interpret parameter 'program {options}'. "
                   "Expected 'program gaussian16' or 'program gaussian09'.")
        if options.casefold() == "gaussian16":
            program_state.program_id = enums.ProgramID.GAUSSIAN_16
        elif options.casefold() == "gaussian09":
            program_state.program_id = enums.ProgramID.GAUSSIAN_09
        else:
            raise exceptions.InputError(err_msg)

    @staticmethod
    def random_seed(options, program_state):
        """Reset program_state.random seed from options."""
        err_msg = (f"Could not interpret parameter 'random_seed {options}'. "
                   "Expected 'random_seed integer' or 'random_seed generate'.")
        if options.casefold() == "generate":
            program_state.random.reset_seed()
        else:
            try:
                program_state.random.reset_seed(int(options))
            except ValueError:
                raise exceptions.InputError(err_msg)

    @staticmethod
    def rotational_energy(options, program_state):
        """Populate program_state.add_rotational_energy from options."""
        err_msg = ("Could not interpret parameter 'rotational_energy "
                   f" {options}. Expected 'rotational_energy on' or "
                   "'rotational_energy off'.")
        if options.casefold() == "on":
            program_state.add_rotational_energy = enums.RotationalEnergy.YES
        elif options.casefold() == "off":
            program_state.add_rotational_energy = enums.RotationalEnergy.NO
        else:
            raise exceptions.InputError(err_msg)

    @staticmethod
    def step_size(options, program_state):
        """Populate program_state.step_size from options."""
        err_msg = (f"Could not interpret parameter 'step_size {options}'. "
                   "Expected 'step_size floating-point'.")
        try:
            program_state.step_size = \
                containers.Time(float(options), enums.TimeUnits.FEMTOSECOND)
        except ValueError:
            raise exceptions.InputError(err_msg)

    @staticmethod
    def temperature(options, program_state):
        """Populate program_state.temperature from options."""
        err_msg = (f"Could not interpret parameter 'temperature {options}'. "
                   "Expected 'temperature floating-point'.")
        try:
            program_state.temperature = float(options)
        except ValueError:
            raise exceptions.InputError(err_msg)


def main():
    """Parse input from stdin to check input file validity."""
    import os
    import sys
    from milo_1_0_1 import program_state as ps

    stdout = sys.stdout
    null_output = open(os.devnull, 'w')
    sys.stdout = null_output

    program_state = ps.ProgramState()

    try:
        parse_input(sys.stdin, program_state)
    except Exception as e:
        sys.stdout = stdout
        print("Input file is NOT valid.")
        print(e)
    else:
        sys.stdout = stdout
        print("Input file is valid.")

    null_output.close()

    return program_state


if __name__ == "__main__":
    main()
