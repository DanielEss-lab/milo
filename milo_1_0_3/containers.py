#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Define classes used to extend list functionality to make units easy."""

from milo_1_0_3 import enumerations as enums
from milo_1_0_3 import scientific_constants as sc
from milo_1_0_3 import atom


class Positions:
    """
    Extend a list to make xyz positional data easy to use.

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
        if units is enums.DistanceUnits.ANGSTROM:
            factor = 1
        elif units is enums.DistanceUnits.BOHR:
            factor = sc.BOHR_TO_ANGSTROM
        elif units is enums.DistanceUnits.METER:
            factor = sc.METER_TO_ANGSTROM
        else:
            raise ValueError(f"Unknown positions unit: {units}")
        self._positions[index] = tuple(i * factor for i in (x, y, z))

    def append(self, x, y, z, units):
        """Append x, y, z as a tuple to the end of the list."""
        if units is enums.DistanceUnits.ANGSTROM:
            factor = 1
        elif units is enums.DistanceUnits.BOHR:
            factor = sc.BOHR_TO_ANGSTROM
        elif units is enums.DistanceUnits.METER:
            factor = sc.METER_TO_ANGSTROM
        else:
            raise ValueError(f"Unknown positions unit: {units}")
        self._positions.append(tuple(i * factor for i in (x, y, z)))

    def as_angstrom(self, index=None):
        """Return the entire list or specific index in angstroms."""
        if index is None:
            return self._positions
        return self._positions[index]

    def as_bohr(self, index=None):
        """Return the entire list or specific index in bohr radii."""
        factor = sc.ANGSTROM_TO_BOHR
        if index is None:
            return [tuple(i * factor for i in j) for j in self._positions]
        return tuple(i * factor for i in self._positions[index])

    def as_meter(self, index=None):
        """Return the entire list or specific index in meters."""
        factor = sc.ANGSTROM_TO_METER
        if index is None:
            return [tuple(i * factor for i in j) for j in self._positions]
        return tuple(i * factor for i in self._positions[index])

    @classmethod
    def from_velocity(cls, velocities, change_in_time):
        """
        Return a displacement given a velocity and a time difference.

        Δx = v * Δt
        """
        if not (isinstance(velocities, Velocities) and
                isinstance(change_in_time, Time)):
            raise TypeError(f"Cannot create Positions object from "
                f"'{type(velocities)}' and '{type(change_in_time)}'.")
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
        if not (isinstance(acceleration, Accelerations) and
                isinstance(change_in_time, Time)):
            raise TypeError(f"cannot create Positions object from "
                f"'{type(acceleration)}' and '{type(change_in_time)}'")
        displacement = cls()
        dt2 = change_in_time.as_second() ** 2
        for (x, y, z) in acceleration.as_meter_per_sec_sqrd():
            displacement.append(x * dt2, y * dt2, z * dt2,
                                units=enums.DistanceUnits.METER)
        return displacement

    def __add__(self, other):
        """Add operator functionality."""
        if not isinstance(other, Positions):
            return NotImplemented
        total = Positions()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_angstrom(),
                                              other.as_angstrom()):
            total.append(x1 + x2, y1 + y2, z1 + z2,
                         enums.DistanceUnits.ANGSTROM)
        return total

    def __sub__(self, other):
        """Subtraction operator functionality."""
        if not isinstance(other, Positions):
            return NotImplemented
        difference = Positions()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_angstrom(),
                                              other.as_angstrom()):
            difference.append(x1 - x2, y1 - y2, z1 - z2,
                              enums.DistanceUnits.ANGSTROM)
        return difference

    def __mul__(self, other):
        """Multiplication operator functionality."""
        if not isinstance(other, (int, float)):
            return NotImplemented
        products = Positions()
        for (x, y, z) in self.as_angstrom():
            products.append(x * other, y * other, z * other,
                            enums.DistanceUnits.ANGSTROM)
        return products

    __rmul__ = __mul__

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for x, y, z in self._positions:
            strings.append(f"{x:10.6f} {y:10.6f} {z:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return f"<Positions object with {len(self._positions)} atoms.>"


class Velocities:
    """
    Extend a list to make xyz momentum data easy to use.

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
            factor = 1
        elif units is enums.VelocityUnits.ANGSTROM_PER_FS:
            factor = sc.ANGSTROM_TO_METER * (1 / sc.FEMTOSECOND_TO_SECOND)
        elif units is enums.VelocityUnits.ANGSTROM_PER_SEC:
            factor = sc.ANGSTROM_TO_METER
        else:
            raise ValueError(f"Unknown velocity unit: {units}")
        self._velocities.append(tuple(i * factor for i in (x, y, z)))

    def as_meter_per_sec(self, index=None):
        """Return the entire list or specific index in m/s."""
        if index is None:
            return self._velocities
        return self._velocities[index]

    def as_angstrom_per_fs(self, index=None):
        """Return the entire list or specific index in angstrom/fs."""
        factor = (sc.METER_TO_ANGSTROM * (1 / sc.SECOND_TO_FEMTOSECOND))
        if index is None:
            return [tuple(i * factor for i in j) for j in self._velocities]
        return tuple(i * factor for i in self._velocities[index])

    def as_angstrom_per_sec(self, index=None):
        """Return the entire list or specific index in angstrom/s."""
        factor = sc.METER_TO_ANGSTROM
        if index is None:
            return [tuple(i * factor for i in j) for j in self._velocities]
        return tuple(i * factor for i in self._velocities[index])

    @classmethod
    def from_acceleration(cls, acceleration, change_in_time):
        """
        Return a velocity given an acceleration and a time difference.

        Δv = a*Δt
        """
        if not (isinstance(acceleration, Accelerations) and
                isinstance(change_in_time, Time)):
            raise TypeError(f"Cannot create Velocities object from "
                f"'{type(acceleration)}' and '{type(change_in_time)}'.")
        displacement = cls()
        dt = change_in_time.as_second()
        for (x, y, z) in acceleration.as_meter_per_sec_sqrd():
            displacement.append(x * dt, y * dt, z * dt,
                                units=enums.VelocityUnits.METER_PER_SEC)
        return displacement

    def __add__(self, other):
        """Add operator functionality."""
        if not isinstance(other, Velocities):
            return NotImplemented
        total = Velocities()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_meter_per_sec(),
                                              other.as_meter_per_sec()):
            total.append(x1 + x2, y1 + y2, z1 + z2,
                         enums.VelocityUnits.METER_PER_SEC)
        return total

    def __sub__(self, other):
        """Subtraction operator functionality."""
        if not isinstance(other, Velocities):
            return NotImplemented
        difference = Velocities()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_meter_per_sec(),
                                              other.as_meter_per_sec()):
            difference.append(x1 - x2, y1 - y2, z1 - z2,
                              enums.VelocityUnits.METER_PER_SEC)
        return difference

    def __mul__(self, other):
        """Multiplication operator functionality."""
        if not isinstance(other, (int, float)):
            return NotImplemented
        products = Velocities()
        for (x, y, z) in self.as_meter_per_sec():
            products.append(x * other, y * other, z * other,
                            enums.VelocityUnits.METER_PER_SEC)
        return products

    __rmul__ = __mul__

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for x, y, z in self._velocities:
            strings.append(f"{x:10.6f} {y:10.6f} {z:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return f"<Velocities object with {len(self._velocities)} atoms.>"


class Accelerations:
    """
    Extend a list to make xyz acceleration data easy to use.

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
            factor = 1
        else:
            raise ValueError(f"Unknown Acceleration units: {units}")
        self._accelerations.append(tuple(i * factor for i in (x, y, z)))

    def as_meter_per_sec_sqrd(self, index=None):
        """Return the entire list or specific index in m/s^2."""
        if index is None:
            return self._accelerations
        return self._accelerations[index]

    @classmethod
    def from_forces(cls, forces, atoms):
        """
        Return an acceleration given forces and a masses.

        F = ma --> a = F / m
        """
        if not (isinstance(forces, Forces) and isinstance(atoms, list) and
                isinstance(atoms[0], atom.Atom)):
            raise TypeError(f"Cannot create Accelerations object from"
                f"'{type(forces)}' and '{type(atoms)}'.")
        accelerations = cls()
        atom_masses_kg = [element.mass * sc.AMU_TO_KG for element in atoms]
        for mass, (x, y, z) in zip(atom_masses_kg, forces.as_newton()):
            accelerations.append(x / mass, y / mass, z / mass,
                                 enums.AccelerationUnits.METER_PER_SEC_SQRD)
        return accelerations

    def __add__(self, other):
        """Add operator functionality."""
        if not isinstance(other, Accelerations):
            return NotImplemented
        total = Accelerations()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_meter_per_sec_sqrd(),
                                              other.as_meter_per_sec_sqrd()):
            total.append(x1 + x2, y1 + y2, z1 + z2,
                         enums.AccelerationUnits.METER_PER_SEC_SQRD)
        return total

    def __sub__(self, other):
        """Subtraction operator functionality."""
        if not isinstance(other, Accelerations):
            return NotImplemented
        difference = Accelerations()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_meter_per_sec_sqrd(),
                                              other.as_meter_per_sec_sqrd()):
            difference.append(x1 - x2, y1 - y2, z1 - z2,
                              enums.AccelerationUnits.METER_PER_SEC_SQRD)
        return difference

    def __mul__(self, other):
        """Multiplication operator functionality."""
        if not isinstance(other, (int, float)):
            return NotImplemented
        products = Accelerations()
        for (x, y, z) in self.as_meter_per_sec_sqrd():
            products.append(x * other, y * other, z * other,
                            enums.AccelerationUnits.METER_PER_SEC_SQRD)
        return products

    __rmul__ = __mul__

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for x, y, z in self._accelerations:
            strings.append(f"{x:10.6f} {y:10.6f} {z:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return f"<Accelerations object with {len(self._accelerations)} atoms.>"


class Forces:
    """
    Extend a list to make xyz forces data easy to use.

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
            factor = 1
        elif units is enums.ForceUnits.DYNE:
            factor = sc.DYNE_TO_NEWTON
        elif units is enums.ForceUnits.MILLIDYNE:
            factor = sc.FROM_MILLI * sc.DYNE_TO_NEWTON
        elif units is enums.ForceUnits.HARTREE_PER_BOHR:
            factor = sc.HARTREE_PER_BOHR_TO_NEWTON
        else:
            raise ValueError(f"Unknown Force units: {units}")
        self._forces.append(tuple(i * factor for i in (x, y, z)))

    def as_newton(self, index=None):
        """Return the entire list or specific index in newtons."""
        if index is None:
            return self._forces
        return self._forces[index]

    def as_dyne(self, index=None):
        """Return the entire list or specific index in dyne."""
        factor = sc.NEWTON_TO_DYNE
        if index is None:
            return [tuple(i * factor for i in j) for j in self._forces]
        return tuple(i * factor for i in self._forces[index])

    def as_millidyne(self, index=None):
        """Return the entire list or specific index in millidyne."""
        factor = sc.NEWTON_TO_DYNE * sc.TO_MILLI
        if index is None:
            return [tuple(i * factor for i in j) for j in self._forces]
        return tuple(i * factor for i in self._forces[index])

    def as_hartree_per_bohr(self, index=None):
        """Return the entire list or specific index in hartree/bohr."""
        factor = sc.NEWTON_TO_HARTREE_PER_BOHR
        if index is None:
            return [tuple(i * factor for i in j) for j in self._forces]
        return tuple(i * factor for i in self._forces[index])

    def __add__(self, other):
        """Add operator functionality."""
        if not isinstance(other, Forces):
            return NotImplemented
        total = Forces()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_newton(),
                                              other.as_newton()):
            total.append(x1 + x2, y1 + y2, z1 + z2,
                         enums.ForceUnits.NEWTON)
        return total

    def __sub__(self, other):
        """Subtraction operator functionality."""
        if not isinstance(other, Forces):
            return NotImplemented
        difference = Forces()
        for (x1, y1, z1), (x2, y2, z2) in zip(self.as_newton(),
                                              other.as_newton()):
            difference.append(x1 - x2, y1 - y2, z1 - z2,
                              enums.ForceUnits.NEWTON)
        return difference

    def __mul__(self, other):
        """Multiplication operator functionality."""
        if not isinstance(other, (int, float)):
            return NotImplemented
        products = Forces()
        for (x, y, z) in self.as_newton():
            products.append(x * other, y * other, z * other,
                            enums.ForceUnits.NEWTON)
        return products

    __rmul__ = __mul__

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for x, y, z in self._forces:
            strings.append(f"{x:10.6f}{y:10.6f}{z:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return f"<Forces object with {len(self._forces)} atoms.>"


class Frequencies:
    """
    Extend a list to make frequency data easy to use.

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
            factor = 1
        else:
            raise ValueError(f"Unknown frequency units: {units}")
        self._frequencies.append(frequency * factor)

    def as_recip_cm(self, index=None):
        """Return the entire list or specific index in cm^-1."""
        if index is None:
            return self._frequencies
        return self._frequencies[index]

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for frequency in self._frequencies:
            strings.append(f"{frequency:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return (f"<Frequency object with {len(self._frequencies)} "
                f"vibrational modes.>")


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
            factor = 1
        elif units is enums.ForceConstantUnits.MILLIDYNE_PER_ANGSTROM:
            factor = (sc.FROM_MILLI * sc.DYNE_TO_NEWTON
                      * (1 / sc.ANGSTROM_TO_METER))
        else:
            raise ValueError(f"Unknown Force Constant Units: {units}")
        self._force_constants.append(force_constant * factor)

    def as_newton_per_meter(self, index=None):
        """Return the entire list or specific index in N/m."""
        if index is None:
            return self._force_constants
        return self._force_constants[index]

    def as_millidyne_per_angstrom(self, index=None):
        """Return the entire list or specific index in millidyne/angstrom."""
        factor = (sc.FROM_MILLI * sc.DYNE_TO_NEWTON
                  * (1 / sc.ANGSTROM_TO_METER))
        if index is None:
            return [i * factor for i in self._force_constants]
        return self._force_constants[index] * factor

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for force_constant in self._force_constants:
            strings.append(f"{force_constant:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return (f"<ForceConstants object with {len(self._force_constants)} "
                f"force constants.>")


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
            factor = 1
        elif units is enums.MassUnits.KILOGRAMS:
            factor = sc.KG_TO_AMU
        elif units is enums.MassUnits.GRAMS:
            factor = sc.TO_KILO * sc.KG_TO_AMU
        else:
            raise ValueError(f"Unknown Mass Units: {units}")
        self._masses.append(mass * factor)

    def as_amu(self, index=None):
        """Return the entire list or specific index in amu."""
        if index is None:
            return self._masses
        return self._masses[index]

    def as_kilogram(self, index=None):
        """Return the entire list or specific index in kg."""
        factor = sc.AMU_TO_KG
        if index is None:
            return [i * factor for i in self._masses]
        return self._masses[index] * factor

    def as_gram(self, index=None):
        """Return the entire list or specific index in grams."""
        factor = sc.AMU_TO_KG * sc.FROM_KILO
        if index is None:
            return [i * factor for i in self._masses]
        return self._masses[index] * factor

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
        else:
            raise ValueError(f"Unknown Time units: {unit}")

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
            factor = 1
        elif units is enums.EnergyUnits.KCAL_PER_MOLE:
            factor = sc.KCAL_PER_MOLE_TO_JOULE
        elif units is enums.EnergyUnits.MILLIDYNE_ANGSTROM:
            factor = sc.MILLIDYNE_ANGSTROM_TO_JOULE
        elif units is enums.EnergyUnits.HARTREE:
            factor = sc.HARTREE_TO_JOULE
        else:
            raise ValueError(f"Unknown Energy Units: {units}")
        self._energies.append(energy * factor)

    def alter_energy(self, index, energy, units):
        """Set the energy at index to new value."""
        if units is enums.EnergyUnits.JOULE:
            factor = 1
        elif units is enums.EnergyUnits.KCAL_PER_MOLE:
            factor = sc.KCAL_PER_MOLE_TO_JOULE
        elif units is enums.EnergyUnits.MILLIDYNE_ANGSTROM:
            factor = sc.MILLIDYNE_ANGSTROM_TO_JOULE
        elif units is enums.EnergyUnits.HARTREE:
            factor = sc.HARTREE_TO_JOULE
        else:
            raise ValueError(f"Unknown Energy Units: {units}")
        self._energies[index] = energy * factor

    def as_joules(self, index=None):
        """Return the entire list or specific index in joules."""
        if index is None:
            return self._energies
        return self._energies[index]

    def as_kcal_per_mole(self, index=None):
        """Return the entire list or specific index in kcal/mol."""
        factor = sc.JOULE_TO_KCAL_PER_MOLE
        if index is None:
            return [i * factor for i in self._energies]
        return self._energies[index] * factor

    def as_millidyne_angstrom(self, index=None):
        """Return the entire list or specific index in mdyne*angstrom."""
        factor = sc.JOULE_TO_MILLIDYNE_ANGSTROM
        if index is None:
            return [i * factor for i in self._energies]
        return self._energies[index] * factor

    def as_hartree(self, index=None):
        """Return the entire list or specific index in hartrees."""
        factor = sc.JOULE_TO_HARTREE
        if index is None:
            return [i * factor for i in self._energies]
        return self._energies[index] * factor

    def __str__(self):
        """Return structure as multiline string."""
        strings = []
        for energy in self._energies:
            strings.append(f"{energy:10.6f}")
        return '\n'.join(strings)

    def __repr__(self):
        """Return object representation."""
        return f"<Energies object with {len(self._energies)} energies.>"
