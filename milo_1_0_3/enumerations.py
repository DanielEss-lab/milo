#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Declare all enumerations.

Quick reference: https://docs.python.org/3/library/enum.html
Enums declared in alphabetical order

Suggested usage:
from milo_1_0_3 import enumerations as enums
my_var = enums.ProgramID.GAUSSIAN_16
"""

from enum import Enum


class AccelerationUnits(Enum):
    """Define enum for different units used to measure acceleration."""

    METER_PER_SEC_SQRD = 1


class AngleUnits(Enum):
    """Define enum for different units to measure angles."""

    RADIANS = 1
    DEGREES = 2


class DistanceUnits(Enum):
    """Define enum for different units used to measure distance."""

    ANGSTROM = 1
    BOHR = 2
    METER = 3


class EnergyBoost(Enum):
    """Define enum for energy boost usage."""

    OFF = 1
    ON = 2


class EnergyUnits(Enum):
    """Define enum for different units to measure energy."""

    JOULE = 1
    KCAL_PER_MOLE = 2
    MILLIDYNE_ANGSTROM = 3
    HARTREE = 4


class ForceConstantUnits(Enum):
    """Define enum for different units to measure force constants."""

    NEWTON_PER_METER = 1
    MILLIDYNE_PER_ANGSTROM = 2


class ForceUnits(Enum):
    """Define enum for different units to measure force."""

    NEWTON = 1
    DYNE = 2
    MILLIDYNE = 3
    HARTREE_PER_BOHR = 4


class FrequencyUnits(Enum):
    """Define enum for different units to measure frequency."""

    RECIP_CM = 1


class GeometryDisplacement(Enum):
    """Define enum for initial structure equilibrium displacement."""

    NONE = 1
    EDGE_WEIGHTED = 2
    GAUSSIAN_DISTRIBUTION = 3
    UNIFORM = 4


class MassUnits(Enum):
    """Define enum for different units to measure mass."""

    AMU = 1
    KILOGRAM = 2
    GRAM = 3


class OscillatorType(Enum):
    """Define enum for how to treat harmonic oscillator."""

    QUASICLASSICAL = 1
    CLASSICAL = 2


class PhaseDirection(Enum):
    """Define enum for possible phase directions."""

    BRING_TOGETHER = 1
    PUSH_APART = 2
    RANDOM = 3


class ProgramID(Enum):
    """Define enum for electronic strucutre programs."""

    GAUSSIAN_16 = 1
    G16 = 1  # Creates an alias for GAUSSIAN_16
    GAUSSIAN_09 = 2
    G09 = 2


class PropagationAlgorithm(Enum):
    """Define enum for force propagation algorithms."""

    VERLET = 1
    VELOCITY_VERLET = 2


class RotationalEnergy(Enum):
    """Define enum for whether to add rotational energy."""

    YES = 1
    NO = 2


class TimeUnits(Enum):
    """Define enum for different units to measure time."""

    SECOND = 1
    FEMTOSECOND = 2


class VelocityUnits(Enum):
    """Define enum for different units to measure velocity."""

    METER_PER_SEC = 1
    ANGSTROM_PER_FS = 2
    ANGSTROM_PER_SEC = 3
