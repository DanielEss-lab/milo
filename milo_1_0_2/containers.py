#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Define classes used to extend list functionality to make units easy."""

from milo_1_0_2 import enumerations as enums
from milo_1_0_2 import scientific_constants as sc
from milo_1_0_2 import atom


class Positions:
    """
    Extend a list to make x, y, z positional data easy to use.

    Always stored as angstrom
    """

    def __init__(self):
        """Create the internal list that holds positions."""
        self._positions = list()

    def __len__(self):
        """Return the length of _positions."""
        return len(self._positions)

    def alter_position(self, index, x, y, z, units):
        """Set the xyz position at index to new tuple."""
        if index >= self.__len__():
            raise IndexError("Position to be altered is outside range.")
        if units is enums.DistanceUnits.ANGSTROM:
            self._positions[index] = (x, y, z)
        elif units is enums.DistanceUnits.BOHR:
            x *= sc.BOHR_TO_ANGSTROM
            y *= sc.BOHR_TO_ANGSTROM
            z *= sc.BOHR_TO_ANGSTROM
            self._positions[index] = (x, y, z)
        elif units is enums.DistanceUnits.METER:
            x *= sc.METER_TO_ANGSTROM
            y *= sc.METER_TO_ANGSTROM
            z *= sc.METER_TO_ANGSTROM
            self._positions[index] = (x, y, z)
        else:
            raise ValueError(f"Unknown positions unit: {units}")

    def append(self, x, y, z, units):
        """Append x, y, z as a tuple to the end of the list."""
        if units is enums.DistanceUnits.ANGSTROM:
            self._positions.append((x, y, z))
        elif units is enums.DistanceUnits.BOHR:
            x *= sc.BOHR_TO_ANGSTROM
            y *= sc.BOHR_TO_ANGSTROM
            z *= sc.BOHR_TO_ANGSTROM
            self._positions.append((x, y, z))
        elif units is enums.DistanceUnits.METER:
            x *= sc.METER_TO_ANGSTROM
            y *= sc.METER_TO_ANGSTROM
            z *= sc.METER_TO_ANGSTROM
            self._positions.append((x, y, z))
        else:
            raise ValueError(f"Unknown positions unit: {units}")

    def as_angstrom(self, index=None):
        """Return the entire list or specific index in angstroms."""
        if index is None:
            return self._positions
        else:
            return self._positions[index]

    def as_bohr(self, index=None):
        """Return the entire list or specific index in bohr radii."""
        if index is None:
            conversion = list()
            for x, y, z in self._positions:
                conversion.append((x * sc.ANGSTROM_TO_BOHR,
                                   y * sc.ANGSTROM_TO_BOHR,
                                   z * sc.ANGSTROM_TO_BOHR))
            return conversion
        else:
            return self._positions[index] * sc.ANGSTROM_TO_BOHR

    def as_meter(self, index=None):
        """Return the entire list or specific index in meters."""
        if index is None:
            conversion = list()
            for x, y, z in self._positions:
                conversion.append((x * sc.ANGSTROM_TO_METER,
                                   y * sc.ANGSTROM_TO_METER,
                                   z * sc.ANGSTROM_TO_METER))
            return conversion
        else:
            return self._positions[index] * sc.ANGSTROM_TO_METER

    @classmethod
    def from_velocity(cls, velocities, change_in_time):
        """
        Return a displacement given a velocity and a time difference.

        Δx = v * Δt
        """
        if (type(velocities) is not Velocities
                or type(change_in_time) is not Time):
            raise TypeError(f"cannot create Positions object from"
                            f"'{type(velocities)}' "
                            f"and '{type(change_in_time)}'")

        displacement = cls()
        dt = change_in_time.as_second()
        for (x, y, z) in velocities.as_meter_per_sec():
            displacement.append(x * dt, y * dt, z * dt,
                                units=enums.DistanceUnits.METER)
        return displacement

    @classmethod
    def from_acceleration(cls, acceleration, change_in_time):
        """
        Return a displacement given an acceleration and a time difference.

        Δx = a*Δt^2
        """
        if (type(acceleration) is not Accelerations
                or type(change_in_time) is not Time):
            raise TypeError(f"cannot create Positions object from"
                            f"'{type(acceleration)}' "
                            f"and '{type(change_in_time)}'")

        displacement = cls()
        dt2 = change_in_time.as_second() ** 2
        for (x, y, z) in acceleration.as_meter_per_sec_sqrd():
            displacement.append(x * dt2, y * dt2, z * dt2,
                                units=enums.DistanceUnits.METER)
        return displacement

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for x, y, z in self._positions:
            strings.append(f"{x:10.6f}{y:10.6f}{z:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return f"<Positions object with {len(self._positions)} atoms.>"

    def __add__(self, other):
        """Add operator functionality."""
        if type(other) is not Positions:
            raise TypeError(f"unsupported operand type(s) for +: 'Positions' "
                            f"and '{type(other)}'")
        total = Positions()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_meter(),
                                              other.as_meter()):
            total.append(x1 + x2, y1 + y2, z1 + z2, enums.DistanceUnits.METER)
        return total

    def __sub__(self, other):
        """Subtraction operator functionality."""
        if type(other) is not Positions:
            raise TypeError(f"unsupported operand type(s) for -: 'Positions' "
                            f"and '{type(other)}'")
        difference = Positions()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_meter(),
                                              other.as_meter()):
            difference.append(x1 - x2, y1 - y2, z1 - z2,
                              enums.DistanceUnits.METER)
        return difference

    def __mul__(self, other):
        """Multiplication operator functionality."""
        if type(other) is not int and type(other) is not float:
            raise TypeError(f"can't multiply Positions object by non scalar "
                            f"of type '{type(other)}'")
        products = Positions()
        for (x1, y1, z1) in self.as_meter():
            products.append(x1 * other, y1 * other, z1 * other,
                            enums.DistanceUnits.METER)
        return products


class Velocities:
    """
    Extend a list to make x, y, z momentum data easy to use.

    Always store as meters per second
    """

    def __init__(self):
        """Create the internal list that holds velocities."""
        self._velocities = list()

    def __len__(self):
        """Return the length of _velocities."""
        return len(self._velocities)

    def append(self, x, y, z, units):
        """Append x, y, z as a tuple to the end of the list."""
        if units is enums.VelocityUnits.METER_PER_SEC:
            self._velocities.append((x, y, z))
        elif units is enums.VelocityUnits.ANGSTROM_PER_FS:
            x *= sc.ANGSTROM_TO_METER * (1 / sc.FEMTOSECOND_TO_SECOND)
            y *= sc.ANGSTROM_TO_METER * (1 / sc.FEMTOSECOND_TO_SECOND)
            z *= sc.ANGSTROM_TO_METER * (1 / sc.FEMTOSECOND_TO_SECOND)
            self._velocities.append((x, y, z))
        elif units is enums.VelocityUnits.ANGSTROM_PER_SEC:
            x *= sc.ANGSTROM_TO_METER
            y *= sc.ANGSTROM_TO_METER
            z *= sc.ANGSTROM_TO_METER
            self._velocities.append((x, y, z))
        else:
            raise ValueError(f"Unknown velocity unit: {units}")

    def as_meter_per_sec(self, index=None):
        """Return the entire list or specific index in m/s."""
        if index is None:
            return self._velocities
        else:
            return self._velocities[index]

    def as_angstrom_per_fs(self, index=None):
        """Return the entire list or specific index in angstrom/fs."""
        if index is None:
            conversion = list()
            conversion_factor = (sc.METER_TO_ANGSTROM *
                                (1 / sc.SECOND_TO_FEMTOSECOND))
            for (x, y, z) in self._velocities:
                conversion.append((x * conversion_factor,
                                   y * conversion_factor,
                                   z * conversion_factor))
            return conversion
        else:
            return (self._velocities[index] * sc.METER_TO_ANGSTROM *
                    (1 / sc.SECOND_TO_FEMTOSECOND))

    def as_angstrom_per_sec(self, index=None):
        """Return the entire list or specific index in angstrom/s."""
        if index is None:
            conversion = list()
            conversion_factor = (sc.METER_TO_ANGSTROM)
            for (x, y, z) in self._velocities:
                conversion.append((x * conversion_factor,
                                   y * conversion_factor,
                                   z * conversion_factor))
            return conversion
        else:
            return (self._velocities[index] * sc.METER_TO_ANGSTROM)

    @classmethod
    def from_acceleration(cls, acceleration, change_in_time):
        """
        Return a velocity given an acceleration and a time difference.

        Δv = a*Δt
        """
        if (type(acceleration) is not Accelerations
                or type(change_in_time) is not Time):
            raise TypeError(f"cannot create Velocities object from"
                            f"'{type(acceleration)}' "
                            f"and '{type(change_in_time)}'")

        displacement = cls()
        dt = change_in_time.as_second()
        for (x, y, z) in acceleration.as_meter_per_sec_sqrd():
            displacement.append(x * dt, y * dt, z * dt,
                                units=enums.VelocityUnits.METER_PER_SEC)
        return displacement

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for x, y, z in self._velocities:
            strings.append(f"{x:10.6f}{y:10.6f}{z:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return f"<Velocities object with {len(self._velocities)} atoms.>"

    def __add__(self, other):
        """Add operator functionality."""
        if type(other) is not Velocities:
            raise TypeError(f"unsupported operand type(s) for +: 'Velocities' "
                            f"and '{type(other)}'")
        total = Velocities()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_meter_per_sec(),
                                              other.as_meter_per_sec()):
            total.append(x1 + x2, y1 + y2, z1 + z2,
                         enums.VelocityUnits.METER_PER_SEC)
        return total

    def __sub__(self, other):
        """Subtraction operator functionality."""
        if type(other) is not Positions:
            raise TypeError(f"unsupported operand type(s) for -: 'Velocities' "
                            f"and '{type(other)}'")
        difference = Velocities()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_meter_per_sec(),
                                              other.as_meter_per_sec()):
            difference.append(x1 - x2, y1 - y2, z1 - z2,
                              enums.VelocityUnits.METER_PER_SEC)
        return difference

    def __mul__(self, other):
        """Multiplication operator functionality."""
        if type(other) is not int and type(other) is not float:
            raise TypeError(f"can't multiply Velocities object by non scalar "
                            f"of type '{type(other)}'")
        products = Velocities()
        for (x1, y1, z1) in self.as_meter_per_sec():
            products.append(x1 * other, y1 * other, z1 * other,
                            enums.VelocityUnits.METER_PER_SEC)
        return products


class Accelerations:
    """
    Extend a list to make x, y, z acceleration data easy to use.

    Always store as meters per second squared
    """

    def __init__(self):
        """Create the internal list that holds accelerations."""
        self._accelerations = list()

    def __len__(self):
        """Return the length of _accelerations."""
        return len(self._accelerations)

    def append(self, x, y, z, units):
        """Append x, y, z as a tuple to the end of the list."""
        if units is enums.AccelerationUnits.METER_PER_SEC_SQRD:
            self._accelerations.append((x, y, z))
        else:
            raise ValueError(f"Unknown Acceleration units: {units}")

    def as_meter_per_sec_sqrd(self, index=None):
        """Return the entire list or specific index in m/s^2."""
        if index is None:
            return self._accelerations
        else:
            return self._accelerations[index]

    @classmethod
    def from_forces(cls, forces, atoms):
        """
        Return an acceleration given forces and a masses.

        F = ma --> a = F / m
        """
        if type(forces) is not Forces or (type(atoms) is not list or
                                          type(atoms[0]) is not atom.Atom):
            raise TypeError(f"cannot create Positions object from"
                            f"'{type(forces)}' and '{type(atoms)}'")
        atom_masses = [element.mass * sc.AMU_TO_KG for element in atoms]
        atom_forces = forces.as_newton()

        accelerations = cls()
        for mass, (x, y, z) in zip(atom_masses, atom_forces):
            acceleration = (x / mass, y / mass, z / mass)
            accelerations.append(*acceleration,
                                 enums.AccelerationUnits.METER_PER_SEC_SQRD)
        return accelerations

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for x, y, z in self._accelerations:
            strings.append(f"{x:10.6f}{y:10.6f}{z:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return f"<Accelerations object with {len(self._accelerations)} atoms.>"

    def __add__(self, other):
        """Add operator functionality."""
        if type(other) is not Accelerations:
            return NotImplemented
        total = Accelerations()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_meter_per_sec_sqrd(),
                                              other.as_meter_per_sec_sqrd()):
            total.append(x1 + x2, y1 + y2, z1 + z2,
                         enums.AccelerationUnits.METER_PER_SEC_SQRD)
        return total

    def __sub__(self, other):
        """Subtraction operator functionality."""
        if type(other) is not Positions:
            return NotImplemented
        difference = Accelerations()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_meter_per_sec_sqrd(),
                                              other.as_meter_per_sec_sqrd()):
            difference.append(x1 - x2, y1 - y2, z1 - z2,
                              enums.AccelerationUnits.METER_PER_SEC_SQRD)
        return difference

    def __mul__(self, other):
        """Multiplication operator functionality."""
        if type(other) is not int and type(other) is not float:
            return NotImplemented
        products = Positions()
        for (x1, y1, z1) in self.as_meter_per_sec_sqrd():
            products.append(x1 * other, y1 * other, z1 * other,
                            enums.AccelerationUnits.METER_PER_SEC_SQRD)
        return products


class Forces:
    """
    Extend a list to make x, y, z forces data easy to use.

    Always store as Newtons
    """

    def __init__(self):
        """Create the internal list that holds forces."""
        self._forces = list()

    def __len__(self):
        """Return the length of _forces."""
        return len(self._forces)

    def append(self, x, y, z, units):
        """Append x, y, z as a tuple to the end of the list."""
        if units is enums.ForceUnits.NEWTON:
            self._forces.append((x, y, z))
        elif units is enums.ForceUnits.DYNE:
            x *= sc.DYNE_TO_NEWTON
            y *= sc.DYNE_TO_NEWTON
            z *= sc.DYNE_TO_NEWTON
            self._forces.append((x, y, z))
        elif units is enums.ForceUnits.MILLIDYNE:
            x *= sc.FROM_MILLI * sc.DYNE_TO_NEWTON
            y *= sc.FROM_MILLI * sc.DYNE_TO_NEWTON
            z *= sc.FROM_MILLI * sc.DYNE_TO_NEWTON
            self._forces.append((x, y, z))
        elif units is enums.ForceUnits.HARTREE_PER_BOHR:
            x *= sc.HARTREE_PER_BOHR_TO_NEWTON
            y *= sc.HARTREE_PER_BOHR_TO_NEWTON
            z *= sc.HARTREE_PER_BOHR_TO_NEWTON
            self._forces.append((x, y, z))
        else:
            raise ValueError(f"Unknown Force units: {units}")

    def as_newton(self, index=None):
        """Return the entire list or specific index in newtons."""
        if index is None:
            return self._forces
        else:
            return self._forces[index]

    def as_dyne(self, index=None):
        """Return the entire list or specific index in dyne."""
        if index is None:
            conversion = list()
            for (x, y, z) in self._forces:
                conversion.append((x * sc.NEWTON_TO_DYNE,
                                   y * sc.NEWTON_TO_DYNE,
                                   z * sc.NEWTON_TO_DYNE))
            return conversion
        else:
            return self._forces[index] * sc.NEWTON_TO_DYNE

    def as_millidyne(self, index=None):
        """Return the entire list or specific index in millidyne."""
        if index is None:
            conversion = list()
            for x, y, z in self._forces:
                x *= sc.NEWTON_TO_DYNE * sc.TO_MILLI
                y *= sc.NEWTON_TO_DYNE * sc.TO_MILLI
                z *= sc.NEWTON_TO_DYNE * sc.TO_MILLI
                conversion.append((x, y, z))
            return conversion
        else:
            return self._forces[index] * sc.NEWTON_TO_DYNE * sc.TO_MILLI

    def as_hartree_per_bohr(self, index=None):
        """Return the entire list or specific index in hartree/bohr."""
        if index is None:
            conversion = list()
            for (x, y, z) in self._forces:
                conversion.append((x * sc.NEWTON_TO_HARTREE_PER_BOHR,
                                   y * sc.NEWTON_TO_HARTREE_PER_BOHR,
                                   z * sc.NEWTON_TO_HARTREE_PER_BOHR))
            return conversion
        else:
            return self._forces[index] * sc.NEWTON_TO_HARTREE_PER_BOHR

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for x, y, z in self._forces:
            strings.append(f"{x:10.6f}{y:10.6f}{z:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return f"<Forces object with {len(self._forces)} atoms.>"

    def __add__(self, other):
        """Add operator functionality."""
        if type(other) is not Forces:
            return NotImplemented
        total = Forces()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_newton(),
                                              other.as_newton()):
            total.append(x1 + x2, y1 + y2, z1 + z2,
                         enums.ForceUnits.NEWTON)
        return total

    def __sub__(self, other):
        """Subtraction operator functionality."""
        if type(other) is not Positions:
            return NotImplemented
        difference = Forces()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_newton(),
                                              other.as_newton()):
            difference.append(x1 - x2, y1 - y2, z1 - z2,
                              enums.ForceUnits.NEWTON)
        return difference

    def __mul__(self, other):
        """Multiplication operator functionality."""
        if type(other) is not int or type(other) is not float:
            return NotImplemented
        products = Positions()
        for (x1, y1, z1) in self.as_newton():
            products.append(x1 * other, y1 * other, z1 * other,
                            enums.ForceUnits.NEWTON)
        return products


class Frequencies:
    """
    Extend a list to make frequency data easy to use.

    Stores modes and displacements.
    Always store as RECIP_CM
    """

    def __init__(self):
        """Create the internal list that holds forces."""
        self._frequencies = list()

    def __len__(self):
        """Return the length of _frequencies."""
        return len(self._frequencies)

    def append(self, frequency, units):
        """Append frequency to the end of the list."""
        if units is enums.FrequencyUnits.RECIP_CM:
            self._frequencies.append(frequency)
        else:
            raise ValueError(f"Unknown frequency units: {units}")

    def as_recip_cm(self, index=None):
        """Return the entire list or specific index in cm^-1."""
        if index is None:
            return self._frequencies
        else:
            return self._frequencies[index]

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for frequency in self._frequencies:
            strings.append(f"{frequency:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return f"<Frequency object with {len(self._frequencies)} " \
               f"vibrational modes.>"


class ForceConstants:
    """
    Extend a list to make force constant data easy to use.

    Always store as Newtons per meter
    """

    def __init__(self):
        """Create the internal list that holds forces."""
        self._force_constants = list()

    def __len__(self):
        """Return the length of _force_constants."""
        return len(self._force_constants)

    def append(self, force_constant, units):
        """Append force_constant to the end of the list."""
        if units is enums.ForceConstantUnits.NEWTON_PER_METER:
            self._force_constants.append(force_constant)
        elif units is enums.ForceConstantUnits.MILLIDYNE_PER_ANGSTROM:
            self._force_constants.append(force_constant * sc.FROM_MILLI
                                         * sc.DYNE_TO_NEWTON
                                         * (1 / sc.ANGSTROM_TO_METER))
        else:
            raise ValueError(f"Unknown Force Constant Units: {units}")

    def as_newton_per_meter(self, index=None):
        """Return the entire list or specific index in N/m."""
        if index is None:
            return self._force_constants
        else:
            return self._force_constants[index]

    def as_millidyne_per_angstrom(self, index=None):
        """Return the entire list or specific index in millidyne/angstrom."""
        if index is None:
            conversion = list()
            for force_constant in self._force_constants:
                conversion.append(force_constant * sc.TO_MILLI
                                  * sc.NEWTON_TO_DYNE
                                  * (1 / sc.METER_TO_ANGSTROM))
            return conversion
        else:
            return (self._force_constants[index] * sc.TO_MILLI
                    * sc.NEWTON_TO_DYNE * (1 / sc.METER_TO_ANGSTROM))

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for force_constant in self._force_constants:
            strings.append(f"{force_constant:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return (f"<ForceConstants object with {len(self._force_constants)} "
                f"vibrational modes.>")


class Masses:
    """
    Extend a list to make force constant data easy to use.

    Always store as amu
    """

    def __init__(self):
        """Create the internal list that holds forces."""
        self._masses = list()

    def __len__(self):
        """Return the length of _masses."""
        return len(self._masses)

    def append(self, mass, units):
        """Append mass as a tuple to the end of the list."""
        if units is enums.MassUnits.AMU:
            self._masses.append(mass)
        elif units is enums.MassUnits.KILOGRAMS:
            self._masses.append(mass * sc.KG_TO_AMU)
        elif units is enums.MassUnits.GRAMS:
            self._masses.append(mass * sc.TO_KILO * sc.KG_TO_AMU)
        else:
            raise ValueError(f"Unknown Mass Units: {units}")

    def as_amu(self, index=None):
        """Return the entire list or specific index in amu."""
        if index is None:
            return self._masses
        else:
            return self._masses[index]

    def as_kilogram(self, index=None):
        """Return the entire list or specific index in kg."""
        if index is None:
            conversion = list()
            for mass in self._masses:
                conversion.append(mass * sc.AMU_TO_KG)
            return conversion
        else:
            return self._masses[index] * sc.AMU_TO_KG

    def as_gram(self, index=None):
        """Return the entire list or specific index in grams."""
        if index is None:
            conversion = list()
            for mass in self._masses:
                conversion.append(mass * sc.AMU_TO_KG * sc.FROM_KILO)
            return conversion
        else:
            return self._masses[index] * sc.AMU_TO_KG * sc.FROM_KILO

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for mass in self._masses:
            strings.append(f"{mass:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return f"<Mass object with {len(self._masses)} masses.>"


class Time:
    """
    Extend an float for easy usage of time data.

    Always store in seconds
    """

    def __init__(self, time, unit=enums.TimeUnits.SECOND):
        """Create a time object."""
        if unit is enums.TimeUnits.SECOND:
            self.time = time
        elif unit is enums.TimeUnits.FEMTOSECOND:
            self.time = time * sc.FEMTOSECOND_TO_SECOND

    def as_second(self):
        """Return the value in seconds."""
        return self.time

    def as_femtosecond(self):
        """Return the value in femtoseconds."""
        return self.time * sc.SECOND_TO_FEMTOSECOND

    def __str__(self):
        """Return structure as string."""
        return f"{self.time} seconds"

    def __repr__(self):
        """Return object representation."""
        return f"Time({self.time})"


class Energies:
    """
    Extend a list to make energy data easy to use.

    Always store as Joules
    """

    def __init__(self):
        """Create the internal list that holds energies."""
        self._energies = list()

    def __len__(self):
        """Return the length of _energies."""
        return len(self._energies)

    def append(self, energy, units):
        """Append energy to the end of the list."""
        if units is enums.EnergyUnits.JOULE:
            self._energies.append(energy)
        elif units is enums.EnergyUnits.KCAL_PER_MOLE:
            self._energies.append(energy * sc.KCAL_PER_MOLE_TO_JOULE)
        elif units is enums.EnergyUnits.MILLIDYNE_ANGSTROM:
            self._energies.append(energy * sc.MILLIDYNE_ANGSTROM_TO_JOULE)
        elif units is enums.EnergyUnits.HARTREE:
            self._energies.append(energy * sc.HARTREE_TO_JOULE)
        else:
            raise ValueError(f"Unknown Energy Units: {units}")

    def alter_energy(self, index, energy, units):
        """Set the energy at index to new value."""
        if index >= len(self._energies):
            raise IndexError("Position to be altered is outside range.")
        if units is enums.EnergyUnits.JOULE:
            self._energies[index] = energy
        elif units is enums.EnergyUnits.KCAL_PER_MOLE:
            self._energies[index] = energy * sc.KCAL_PER_MOLE_TO_JOULE
        elif units is enums.EnergyUnits.MILLIDYNE_ANGSTROM:
            self._energies[index] = energy * sc.MILLIDYNE_ANGSTROM_TO_JOULE
        elif units is enums.EnergyUnits.HARTREE:
            self._energies[index] = energy * sc.HARTREE_TO_JOULE
        else:
            raise ValueError(f"Unknown Energy Units: {units}")

    def as_joules(self, index=None):
        """Return the entire list or specific index in joules."""
        if index is None:
            return self._energies
        else:
            return self._energies[index]

    def as_kcal_per_mole(self, index=None):
        """Return the entire list or specific index in kcal/mol."""
        if index is None:
            conversion = list()
            for energy in self._energies:
                conversion.append(energy * sc.JOULE_TO_KCAL_PER_MOLE)
            return conversion
        else:
            return self._energies[index] * sc.JOULE_TO_KCAL_PER_MOLE

    def as_millidyne_angstrom(self, index=None):
        """Return the entire list or specific index in mdyne*angstrom."""
        if index is None:
            conversion = list()
            for energy in self._energies:
                conversion.append(energy * sc.JOULE_TO_MILLIDYNE_ANGSTROM)
            return conversion
        else:
            return self._energies[index] * sc.JOULE_TO_MILLIDYNE_ANGSTROM

    def as_hartree(self, index=None):
        """Return the entire list or specific index in hartrees."""
        if index is None:
            conversion = list()
            for energy in self._energies:
                conversion.append(energy * sc.JOULE_TO_HARTREE)
            return conversion
        else:
            return self._energies[index] * sc.JOULE_TO_HARTREE

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for energy in self._energies:
            strings.append(f"{energy:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return f"<Energies object with {len(self._energies)} energies.>"
