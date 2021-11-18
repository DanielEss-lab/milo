#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate initial velocities based on frequency data."""

import math

from milo_1_0_3 import containers
from milo_1_0_3 import enumerations as enums
from milo_1_0_3 import exceptions
from milo_1_0_3 import scientific_constants as sc


def _calculate_zero_point_energies(program_state):
    """
    Calculate the zero point energy for each mode.

    When oscillator_type is set to CLASSICAL, it isn't technically ZPE, but it
    is treated the same, so for simplicity the variable is still named ZPE

    Reference C++ - lines 59.
    """
    zero_point_energies = containers.Energies()
    zpe_sum = 0

    for frequency in program_state.frequencies.as_recip_cm():
        if frequency < 0:
            frequency = 2
        if program_state.oscillator_type == enums.OscillatorType.CLASSICAL:
            energy = 0.5 * sc.h * sc.c * sc.CLASSICAL_SPACING
        else:
            energy = 0.5 * sc.h * sc.c * frequency
        zero_point_energies.append(energy, enums.EnergyUnits.JOULE)
        zpe_sum += energy

    total_zpe = containers.Energies()
    total_zpe.append(zpe_sum, enums.EnergyUnits.JOULE)
    return zero_point_energies, total_zpe


def _sample(zero_point_energies, program_state):
    """
    Determine the vibrational excitation quanta numbers.

    Reference: C++ - line 59.
    """
    vibrational_quantum_numbers = list()

    for f in range(len(program_state.frequencies)):
        if program_state.temperature == 0:
            vibrational_quantum_numbers.append(0)
            continue
        else:
            # zpe_ratio = e^(-2 * zpe / (RT))
            zpe_ratio = (math.exp((-2 *
                zero_point_energies.as_kcal_per_mole(f)) /
                (sc.GAS_CONSTANT_KCAL * program_state.temperature)))
            if zpe_ratio >= 1:
                zpe_ratio = 0.99999999999

            # c_e_fraction is 'cumulative excitation fraction'
            c_e_fraction = 1 - zpe_ratio
            random_number = program_state.random.uniform()

            i = 1
            while i <= 4000 * zpe_ratio + 2:
                if random_number > c_e_fraction:
                    c_e_fraction += (math.pow(zpe_ratio, i) * (1 - zpe_ratio))
                    i += 1
                else:
                    break
            vibrational_quantum_numbers.append(i - 1)

    for mode in program_state.fixed_vibrational_quanta:
        # This needs to rewrite the random calls from above so that the random
        # number generator will produce the same results for everything else.
        vibrational_quantum_numbers[mode - 1] = \
            program_state.fixed_vibrational_quanta[mode]

    return vibrational_quantum_numbers


def _calculate_displacement(zero_point_energies, vibrational_quantum_numbers,
                            program_state):
    """
    Calculate the displacements (shifts) and vibrational energy per mode.

    Reference: C++ - line 95.
    """
    mode_energies = containers.Energies()
    mode_energies_sum = 0

    shifts = list()

    for f in range(len(program_state.frequencies)):
        if (program_state.oscillator_type is
                enums.OscillatorType.QUASICLASSICAL
                and program_state.frequencies.as_recip_cm(f) > 10):
            energy = (zero_point_energies.as_joules(f) *
                      (2 * vibrational_quantum_numbers[f] + 1))
        else:
            energy = (zero_point_energies.as_joules(f) *
                      (2 * vibrational_quantum_numbers[f]))
        mode_energies.append(energy, enums.EnergyUnits.JOULE)
        mode_energies_sum += energy
        max_shift = (math.sqrt(2 * mode_energies.as_millidyne_angstrom(f)
            / program_state.force_constants.as_millidyne_per_angstrom(f)))

        random_weight = 0  # For zero displacement or imaginary/low frequencies
        if program_state.frequencies.as_recip_cm(f) > 10:
            if (program_state.geometry_displacement_type is
                    enums.GeometryDisplacement.EDGE_WEIGHTED):
                random_weight = program_state.random.edge_weighted()
            elif (program_state.geometry_displacement_type is
                    enums.GeometryDisplacement.GAUSSIAN_DISTRIBUTION):
                random_weight = program_state.random.gaussian()
            elif (program_state.geometry_displacement_type is
                    enums.GeometryDisplacement.UNIFORM):
                random_weight = 2 * (program_state.random.uniform() - 0.5)
        shifts.append(max_shift * random_weight)

    total_mode_energy = containers.Energies()
    total_mode_energy.append(mode_energies_sum,
                             enums.EnergyUnits.JOULE)
    return total_mode_energy, shifts, mode_energies


def _energy_boost(total_mode_energy, program_state):
    """
    Boost/drop the energy of the system by changing temperature.

    Returns true if the strucutre needs to be resampled (because the
    temperature changed). Returns false if sampling does not need to be done
    again.

    Reference: C++ - line 138.
    """
    if total_mode_energy.as_kcal_per_mole(0) <= program_state.energy_boost_min:
        program_state.temperature += 5.0
        return True
    elif total_mode_energy.as_kcal_per_mole(0) >= \
            program_state.energy_boost_max:
        program_state.temperature -= 2.0
        return True
    return False


def _geometry_displacement(shifts, program_state):
    """
    Apply the displacements (shifts) to the coordinates in program_state.

    Referece: C++ - line 167.
    """
    for f in range(len(program_state.frequencies)):
        for j in range(program_state.number_atoms):
            current_x, current_y, current_z = \
                program_state.structures[0].as_angstrom(j)
            mode_x, mode_y, mode_z = \
                program_state.mode_displacements[f].as_angstrom(j)
            new_x = current_x + (mode_x * shifts[f])
            new_y = current_y + (mode_y * shifts[f])
            new_z = current_z + (mode_z * shifts[f])
            program_state.structures[0].alter_position(j, new_x, new_y, new_z,
                enums.DistanceUnits.ANGSTROM)
    return


def _check_if_mode_pushes_apart(program_state):
    """
    Check if the first mode increases the distance between the phase atoms.

    Reference: C++ line 185.
    """
    # Go from 1-based index to 0-based index
    atom1_id = program_state.phase[0] - 1
    atom2_id = program_state.phase[1] - 1

    atom1_modes = program_state.mode_displacements[0].as_angstrom(atom1_id)
    atom2_modes = program_state.mode_displacements[0].as_angstrom(atom2_id)
    atom1_position = program_state.structures[0].as_angstrom(atom1_id)
    atom2_position = program_state.structures[0].as_angstrom(atom2_id)

    before_distance = (pow(atom1_position[0] - atom2_position[0], 2)
                       + pow(atom1_position[1] - atom2_position[1], 2)
                       + pow(atom1_position[2] - atom2_position[2], 2))

    atom1_new_position = (atom1_position[0] + atom1_modes[0],
                          atom1_position[1] + atom1_modes[1],
                          atom1_position[2] + atom1_modes[2])
    atom2_new_position = (atom2_position[0] + atom2_modes[0],
                          atom2_position[1] + atom2_modes[1],
                          atom2_position[2] + atom2_modes[2])

    after_distance = (pow(atom1_new_position[0] - atom2_new_position[0], 2)
                      + pow(atom1_new_position[1] - atom2_new_position[1], 2)
                      + pow(atom1_new_position[2] - atom2_new_position[2], 2))

    return after_distance - before_distance > 0


def _calculate_mode_velocities(mode_energy, shift, program_state):
    """
    Calculate the modal velocity, in a random direction, for each frequency.

    Reference: C++ line 181.
    """
    mode_velocities = list()
    mode_directions = list()

    #          10^-3      *    10^-2      *    10^10             = 10^5
    units = sc.FROM_MILLI * sc.FROM_CENTI * sc.METER_TO_ANGSTROM

    for f in range(len(program_state.frequencies)):
        # kinetic energy units: gram angstrom**2 / s**2
        kinetic_energy = (units * (mode_energy.as_millidyne_angstrom(f) -
            (0.5 * program_state.force_constants.as_millidyne_per_angstrom(f)
            * math.pow(shift[f], 2))))
        if f == 0 and program_state.frequencies.as_recip_cm(f) < 0:
            if program_state.phase_direction is enums.PhaseDirection.RANDOM:
                direction = program_state.random.one_or_neg_one()
            else:
                if _check_if_mode_pushes_apart(program_state):
                    direction = 1
                else:
                    direction = -1
        else:
            direction = program_state.random.one_or_neg_one()
        if (program_state.phase_direction is
                enums.PhaseDirection.BRING_TOGETHER):
            direction *= -1
        if f + 1 in program_state.fixed_mode_directions:
            # This needs to rewrite the random call from above so the random
            # number generator will give the same results for everything else.
            direction = program_state.fixed_mode_directions[f + 1]
        mode_directions.append(direction)
        mode_velocities.append(direction * math.sqrt(2 * kinetic_energy
            / (program_state.reduced_masses.as_amu(f) / sc.AVOGADROS_NUMBER)))
    return mode_velocities, mode_directions


def _calculate_atomic_velocities(mode_velocities, program_state):
    """
    Calculate atomic velocities from modal velocities.

    Referece: C++ line 230.
    """
    atomic_velocities = [[0, 0, 0] for j in range(program_state.number_atoms)]
    for f in range(len(program_state.frequencies)):
        for j in range(program_state.number_atoms):
            for k in range(3):
                atomic_velocities[j][k] += \
                    (program_state.mode_displacements[f].as_angstrom(j)[k]
                     * mode_velocities[f])
    return atomic_velocities


def _calculate_kinetic_energy(atomic_velocities, program_state):
    """
    Calculate the total kinetic energy of the system.

    Referece: C++ line - 308
    """
    kinetic_energy_sum = 0
    units = (sc.AMU_TO_KG * math.pow(sc.ANGSTROM_TO_METER, 2)
            * sc.JOULE_TO_KCAL_PER_MOLE)

    for j in range(program_state.number_atoms):
        kinetic_energy_sum += (0.5 * program_state.atoms[j].mass
                               * (math.pow(atomic_velocities[j][0], 2)
                                  + math.pow(atomic_velocities[j][1], 2)
                                  + math.pow(atomic_velocities[j][2], 2))
                               * units)

    kinetic_energy = containers.Energies()
    kinetic_energy.append(kinetic_energy_sum, enums.EnergyUnits.KCAL_PER_MOLE)
    return kinetic_energy


def _add_rotational_energy(atomic_velocities, program_state):
    """
    Calculate and add the rotational energy.

    Reference: C++ line - 251
    """
    rotateX = []
    rotateY = []
    rotateZ = []
    for x, y, z in range(program_state.structures[0].as_angstrom()):
        rotateX.append([0, (-1 * z), y])
        rotateY.append([z, 0, (-1 * x)])
        rotateZ.append([(-1 * y), x, 0])

    eRotX, eRotY, eRotZ = 0.0, 0.0, 0.0
    step_size = program_state.step_size.as_second()
    units = (sc.AMU_TO_KG * pow(sc.ANGSTROM_TO_METER, 2) *
        sc.JOULE_TO_KCAL_PER_MOLE)

    for j in range(program_state.number_atoms):
        for k in range(3):
            eRotX += (0.5 * program_state.atoms[j].mass * pow(
                rotateX[j][k], 2) / pow(step_size, 2) * units)
            eRotY += (0.5 * program_state.atoms[j].mass * pow(
                rotateY[j][k], 2) / pow(step_size, 2) * units)
            eRotZ += (0.5 * program_state.atoms[j].mass * pow(
                rotateZ[j][k], 2) / pow(step_size, 2) * units)

    kinetic_rotational_X, kinetic_rotational_Y, kinetic_rotational_Z = 0, 0, 0
    if eRotX >= 1:
        kinetic_rotational_X = (math.log(1 - program_state.random.uniform()) *
            -0.5 * sc.GAS_CONSTANT_KCAL * program_state.temperature)
    if eRotY >= 1:
        kinetic_rotational_Y = (math.log(1 - program_state.random.uniform()) *
            -0.5 * sc.GAS_CONSTANT_KCAL * program_state.temperature)
    if eRotZ >= 1:
        kinetic_rotational_Z = (math.log(1 - program_state.random.uniform()) *
            -0.5 * sc.GAS_CONSTANT_KCAL * program_state.temperature)
    rotational_kinetic_energy = containers.Energies()
    rotational_kinetic_energy.append(kinetic_rotational_X +
                                     kinetic_rotational_Y +
                                     kinetic_rotational_Z,
                                     enums.EnergyUnits.KCAL_PER_MOLE)

    signX = program_state.random.one_or_neg_one()
    signY = program_state.random.one_or_neg_one()
    signZ = program_state.random.one_or_neg_one()

    scaleX = math.sqrt(kinetic_rotational_X / eRotX)
    scaleY = math.sqrt(kinetic_rotational_Y / eRotY)
    scaleZ = math.sqrt(kinetic_rotational_Z / eRotZ)

    for j in range(program_state.number_atoms):
        for k in range(3):
            rotateX[j][k] *= scaleX * signX / step_size
            rotateY[j][k] *= scaleY * signY / step_size
            rotateZ[j][k] *= scaleZ * signZ / step_size
    for j in range(program_state.number_atoms):
        for k in range(3):
            atomic_velocities[j][k] += (rotateX[j][k] + rotateY[j][k] +
                                        rotateZ[j][k])
    return rotational_kinetic_energy


def _add_velocities_to_program_state(atomic_velocities, program_state):
    """
    Append initial velocities to program_state.

    Referece: C++ line - 318.
    """
    velocities = containers.Velocities()
    for i in range(len(atomic_velocities)):
        x = atomic_velocities[i][0]
        y = atomic_velocities[i][1]
        z = atomic_velocities[i][2]
        velocities.append(x, y, z, enums.VelocityUnits.ANGSTROM_PER_SEC)
    program_state.velocities.append(velocities)


def generate(program_state):
    """Sample initial energy (kinetic and potential)."""
    zero_point_energies, total_zpe \
        = _calculate_zero_point_energies(program_state)

    vibrational_quantum_numbers = _sample(zero_point_energies, program_state)

    total_mode_energy, shifts, mode_energies = _calculate_displacement(
        zero_point_energies, vibrational_quantum_numbers, program_state)

    # Check energy boost and if needed, re-sample and redo displacements
    print("### Energy Boost -------------------------------------------------")
    if program_state.energy_boost is enums.EnergyBoost.ON:
        print("  Energy boost on")
        print("  Changing temperature and resampling until the vibrational ")
        print(f"  energy is between {program_state.energy_boost_min} and "
              f"{program_state.energy_boost_max} kcal/mol.")
        print()
        if program_state.energy_boost_max < total_zpe.as_kcal_per_mole(0):
            raise exceptions.InputError("Energy Boost max energy is less than "
                                        "ZPE.")
        print("  Attempt   Vibrational Energy (kcal/mol)   Temperature (K)")
        print("  ---------------------------------------------------------")
        i = 1
        print(f"  {i:>7}   {total_mode_energy.as_kcal_per_mole(0):18.6f}"
              f"              {program_state.temperature:11.2f}")
        while _energy_boost(total_mode_energy, program_state):
            i += 1
            vibrational_quantum_numbers = _sample(zero_point_energies,
                                                  program_state)
            total_mode_energy, shifts, mode_energies = \
                _calculate_displacement(zero_point_energies,
                    vibrational_quantum_numbers, program_state)
            print(f"  {i:>7}   {total_mode_energy.as_kcal_per_mole(0):18.6f}"
                  f"              {program_state.temperature:11.2f}")
        print("  Energy boost criteria met")
    else:
        print("  Energy boost off")
    print()

    # Stretch coordinates
    print("### Initial Geometry Displacement --------------------------------")
    if (program_state.geometry_displacement_type is not
            enums.GeometryDisplacement.NONE):
        _geometry_displacement(shifts, program_state)
        print("  Modified initial structure")
        for atom, position in zip(program_state.atoms, program_state
                                  .structures[0].as_angstrom()):
            print(f"    {atom.symbol.ljust(2)} {position[0]:10.6f} "
                  f"{position[1]:10.6f} {position[2]:10.6f}")
    else:
        print("  Geometry displacement turned off. Using input structure for")
        print("  starting geometry.")
    print()

    mode_velocities, mode_directions = _calculate_mode_velocities(
        mode_energies, shifts, program_state)

    print("### Vibrational Quantum Numbers ----------------------------------")
    print("  Mode  Wavenumber  Quantum No.  Energy (kcal/mol)  Mode Direction")
    print("  ----------------------------------------------------------------")
    for i, (mode_energy, quantum_n, frequency, direction) in enumerate(zip(
            mode_energies.as_kcal_per_mole(), vibrational_quantum_numbers,
            program_state.frequencies.as_recip_cm(), mode_directions),
            1):
        print(f"  {i:>4}  {frequency:10.3f}  {quantum_n:>11}  "
              f"{mode_energy:17.6f}  {direction:>14}")
    print()

    print("### Mode Velocities (meters/second) ------------------------------")
    for mode_velocity in mode_velocities:
        mode_velocity *= sc.ANGSTROM_TO_METER
        print(f"  {mode_velocity:15.6e}")
    print()

    atomic_velocities = _calculate_atomic_velocities(mode_velocities,
                                                     program_state)

    vibrational_kinetic_energy = _calculate_kinetic_energy(
        atomic_velocities, program_state)

    print("### Rotational Energy --------------------------------------------")
    if program_state.add_rotational_energy is enums.RotationalEnergy.YES:
        rotational_kinetic_energy = _add_rotational_energy(atomic_velocities,
            program_state)
        print(f"  {rotational_kinetic_energy.as_kcal_per_mole(0):.6f} kcal/"
              "mol rotational energy added.")
        total_kinetic_energy = _calculate_kinetic_energy(atomic_velocities,
            program_state)
    else:
        print("  Rotational energy turned off.")
        rotational_kinetic_energy = containers.Energies()
        rotational_kinetic_energy.append(0, enums.EnergyUnits.KCAL_PER_MOLE)
        total_kinetic_energy = vibrational_kinetic_energy
    print()

    _add_velocities_to_program_state(atomic_velocities, program_state)

    print("### Initial Velocities (meters/second) ---------------------------")
    for atom, velocity in zip(program_state.atoms,
                              program_state.velocities[-1].as_meter_per_sec()):
        print(f"  {atom.symbol.ljust(2)} {velocity[0]:15.6e} "
              f"{velocity[1]:15.6e} {velocity[2]:15.6e}")
    print()

    excitation_energy = containers.Energies()
    excitation_energy.append(total_mode_energy.as_kcal_per_mole(0)
        - total_zpe.as_kcal_per_mole(0), enums.EnergyUnits.KCAL_PER_MOLE)

    print("### Initial Energy Sampling Summary (kcal/mol) -------------------")
    print("  Zero point energy:")
    print(f"  {total_zpe.as_kcal_per_mole(0):11.6f}")
    print("  Excitation energy:")
    print(f"  {excitation_energy.as_kcal_per_mole(0):11.6f}")
    print("  Quantum vibrational energy (zpe + excitation):")
    print(f"  {total_mode_energy.as_kcal_per_mole(0):11.6f}")
    print("  Vibrational component of kinetic energy:")
    print(f"  {vibrational_kinetic_energy.as_kcal_per_mole(0):11.6f}")
    print("  Rotation component of kinetic energy:")
    print(f"  {rotational_kinetic_energy.as_kcal_per_mole(0):11.6f}")
    print("  Total kinetic energy:")
    print(f"  {total_kinetic_energy.as_kcal_per_mole(0):11.6f}")
    print()
