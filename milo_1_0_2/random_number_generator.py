#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Provide extended functionality from standard random."""

import math
import os
import random


class RandomNumberGenerator():
    """
    Provide extended functionality from standard random.

    This was created to ensure that all random calls throughout the program
    used the same seed, which could then be output to allow reproducibility
    between different runs.
    """

    def __init__(self, seed=None):
        """
        Create object from seed, generating new seed if seed is 0 or not given.

        The seed it used to initialize a random.Random() object. If the seed is
        0 or not given, a new seed will be generated. It will first attempt to
        generate the seed from os.urandom (which is platform dependent but
        usually pulls from /dev/random/). If this fails due to a
        NotImplementedError, it will create a seed using system time and the
        process id.
        """
        self.seed = seed
        if self.seed is None:
            try:
                self.seed = int.from_bytes(os.urandom(5),
                                           byteorder='big', signed=False)
            except NotImplementedError:
                import time
                current_time = str(int(time.time()))
                process_id = str(os.getpid())
                self.seed = int(process_id
                                + current_time[-(12 - len(process_id)):])
        self.rng = random.Random(self.seed)

    def reset_seed(self, seed=None):
        """Reset the seed and create a new random.Random() object with it."""
        self.__init__(seed)

    def uniform(self):
        """Return a uniformly distributed random number between 0 and 1."""
        return self.rng.random()

    def edge_weighted(self):
        """
        Return an edge weighted random number between -1 and 1.

        Implemented by taking the sin of a random number from a uniform
        distribution between 0 and 2*pi.
        """
        return math.sin(2 * math.pi * self.rng.random())

    def gaussian(self):
        """
        Return a random number [-1, 1] from a modified normal distribution.

        The random number is returned from a normal distribution with mu = 0
        and, sigma = 1/sqrt(2), only sampled between -1 and 1, inclusive.

        Sigma was chosen to be 1/sqrt(2) so the probability of the distribution
        being greater than 1 or less than -1 would be the same as the
        probability of being outside the classically allowed region of the
        quantum harmonic oscillator. References: https://physicspages.com/pdf/G
        riffiths%20QM/Griffiths%20Problems%2002.15.pdf; https://phys.libretexts
        .org/Bookshelves/University_Physics/Book%3A_University_Physics_(OpenSta
        x)/Map%3A_University_Physics_III_-_Optics_and_Modern_Physics_(OpenStax)
        /07%3A_Quantum_Mechanics/7.06%3A_The_Quantum_Harmonic_Oscillator
        """
        mu = 0
        sigma = 0.7071067811865475  # 1 / sqrt(2)
        possible_number = self.rng.gauss(mu, sigma)
        while possible_number < -1 or possible_number > 1:
            possible_number = self.rng.gauss(mu, sigma)
        return possible_number

    def one_or_neg_one(self):
        """Return one or negative one randomly with equal probability."""
        if self.rng.random() >= 0.5:
            return 1
        return -1
