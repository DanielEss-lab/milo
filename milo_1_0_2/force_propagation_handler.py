#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculate the next position using force data from previous step."""

from milo_1_0_2 import containers
from milo_1_0_2 import enumerations as enums


def get_propagation_handler(program_state):
    """Return the correct propagation handler class based on program state."""
    if program_state.propagation_algorithm is \
            enums.PropagationAlgorithm.VERLET:
        return Verlet
    elif program_state.propagation_algorithm is \
            enums.PropagationAlgorithm.VELOCITY_VERLET:
        return VelocityVerlet
    else:
        raise ValueError(f"Unknown force propagation algorithm:"
                         f"{program_state.propagation_algorithm}")


class ForcePropagationHandler:
    """Template class for force propagation handlers."""

    @staticmethod
    def _calculate_acceleration(program_state):
        """Calculate acceleration with F = m*a."""
        acceleration = containers.Accelerations.from_forces(
            program_state.forces[-1], program_state.atoms)
        program_state.accelerations.append(acceleration)

    @staticmethod
    def _calculate_velocity(program_state):
        """Calculate velocity with v(n) = v(n-1) + 1/2*(a(n-1) + a(n))*dt."""
        velocity = (program_state.velocities[-1] +
                    (containers.Velocities.from_acceleration(
                        (program_state.accelerations[-1] +
                            program_state.accelerations[-2]),
                        program_state.step_size) * 0.5))
        program_state.velocities.append(velocity)


class Verlet(ForcePropagationHandler):
    """
    Calculates the next structure using the Verlet algorithm.

    Each pass calculates the acceleration and velocity of the previous time
    step, then calculates the new positions for the current time step.

    Velocities are calculated in this algorithm only to be output. They are
    not used in the actual force propagation algorithm.

    Verlet algorithm equations:
        a(n-1) = F*m
        x(n) = x(n-1) + v(n-1)*dt + 1/2*a(n-1)*(dt^2);    when n == 1
             = 2*x(n-1) - x(n-2) + a(n-1)*(dt^2);         when n >= 2

        v(n-1) = v(n-2) + 1/2*(a(n-2) + a(n-1))*dt;       only used for output
    """

    @classmethod
    def run_next_step(cls, program_state):
        """Calculate the next structure using the Verlet algorithm."""
        # Calculate a(n-1) from forces
        cls._calculate_acceleration(program_state)

        # Calculate v(n-1) using equation above (only used for output)
        if len(program_state.structures) > 1:
            cls._calculate_velocity(program_state)

        # Calculate x(n) using one of the equations above
        if len(program_state.structures) == 1:
            structure = (program_state.structures[-1] +
                         containers.Positions.from_velocity(
                             program_state.velocities[-1],
                             program_state.step_size) +
                         (containers.Positions.from_acceleration(
                             program_state.accelerations[-1],
                             program_state.step_size) * 0.5))
            program_state.structures.append(structure)
        elif len(program_state.structures) >= 2:
            structure = ((program_state.structures[-1] * 2) -
                         program_state.structures[-2] +
                         containers.Positions.from_acceleration(
                             program_state.accelerations[-1],
                             program_state.step_size))
            program_state.structures.append(structure)


class VelocityVerlet(ForcePropagationHandler):
    """
    Calculates the next structure using the Velocity Verlet algorithm.

    Each pass calculates the acceleration and velocity of the previous time
    step, then calculates the new positions for the current time step.

    On the first pass, the velocities are not calculated and instead come from
    initial_energy_sampler or the input file.

    Velocity Verlet algorithm equations:
        a(n-1) = F*m
        v(n-1) = v(n-2) + 1/2*(a(n-2) + a(n-1))*dt
        x(n) = x(n-1) + v(n-1)*dt + 1/2*a(n-1)*dt^2
    """

    @classmethod
    def run_next_step(cls, program_state):
        """Calculate the next structure using the Velocity Verlet algorithm."""
        # Calculate a(n-1) from forces
        cls._calculate_acceleration(program_state)

        # Calculate v(n-1) using equation above
        if len(program_state.structures) > 1:
            cls._calculate_velocity(program_state)

        # Calculate x(n) using equation above
        structure = (program_state.structures[-1] +
                     containers.Positions.from_velocity(
                         program_state.velocities[-1],
                         program_state.step_size) +
                     (containers.Positions.from_acceleration(
                         program_state.accelerations[-1],
                         program_state.step_size) * 0.5))
        program_state.structures.append(structure)
