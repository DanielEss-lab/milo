#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Store data used throughout the simulation."""

from milo_1_0_1 import containers
from milo_1_0_1 import enumerations as enums
from milo_1_0_1 import random_number_generator as rng


class ProgramState:
    """Contain all the data used throughout the simulation."""

    def __init__(self):
        """Create all variables and sets some to default values."""
        # Basic job data
        self.job_name = None
        self.number_atoms = None
        self.spin = None
        self.charge = None
        self.atoms = list()
        self.current_step = 0  # current_step is 0 for the first step
        self.step_size = containers.Time(1.00, enums.TimeUnits.FEMTOSECOND)
        self.max_steps = None  # If None, no limit
        self.temperature = 298.15  # kelvin
        self.input_structure = containers.Positions()

        # Enums for job control
        self.propagation_algorithm = enums.PropagationAlgorithm.VERLET
        self.oscillator_type = enums.OscillatorType.QUASICLASSICAL
        self.add_rotational_energy = enums.RotationalEnergy.NO
        self.geometry_displacement_type = enums.GeometryDisplacement.NONE
        self.phase_direction = enums.PhaseDirection.RANDOM
        self.phase = None  # (1, 2)
        self.fixed_mode_directions = dict()  # maps mode number to 1 or -1
        self.fixed_vibrational_quanta = dict()  # maps mode number to quanta

        # Data added to at every time step
        self.structures = list()
        self.velocities = list()
        self.forces = list()
        self.accelerations = list()
        self.energies = list()

        # Frequency data
        self.frequencies = containers.Frequencies()
        self.mode_displacements = list()  # mode[frequency index][atom index]
        self.force_constants = containers.ForceConstants()
        self.reduced_masses = containers.Masses()
        self.zero_point_energy = None
        self.zero_point_correction = None

        # Energy Boost
        self.energy_boost = enums.EnergyBoost.OFF
        self.energy_boost_min = None
        self.energy_boost_max = None

        # Random Number Generation
        # self.random_seed = 0  # this is in self.random as self.random.seed
        self.random = rng.RandomNumberGenerator()

        # Program information
        self.program_id = enums.ProgramID.GAUSSIAN_16
        self.gaussian_header = None  # "m06 6-31g(d,p)"
        self.gaussian_footer = None
        self.processor_count = None
        self.memory_amount = None

        # Final output files
        self.output_xyz_file = True
