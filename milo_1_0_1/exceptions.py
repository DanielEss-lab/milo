#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Define Milo specific exceptions."""


class ElectronicStructureProgramError(Exception):
    """Raise when the ESP returns an error."""

    def __init__(self, message):
        """Create an exception with given message."""
        self.message = message


class InputError(Exception):
    """Raise when input file is not valid."""

    def __init__(self, message):
        """Create an exception with given message."""
        self.message = message
