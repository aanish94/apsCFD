#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Thermodynamic Functions/Properties of a gas mixture in chemical equilibrium
"""

import cantera as ct
import math


def initialize_gas_object(name):
    """Initialize Cantera solution object

    Inputs are listed:
        - https://github.com/Cantera/cantera/tree/master/data/inputs

    :param <str> name: Solution name

    :return <cantera.composite.Solution> object
    """

    name = name.lower().strip()

    if not name.endswith('.cti'):
        name += '.cti'

    # TODO: Validate name
    if name not in []:
        pass

    return ct.Solution(name)
    # return ct.ThermoPhase(name)


def speed_of_sound(gas, rtol=1.0e-6, maxiter=5000):
    """Compute the equilbrium speed of sound for a gas

    :param <cantera.composite.Solution> gas: Solution object

    :return <float> Gas speed of sound w/ an equilibrium composition (m/s)
    """

    # Set the gas to equilibrium at its current T and P
    gas.equilibrate('TP', rtol=rtol, maxiter=maxiter)

    # Store properties
    s0 = gas.s
    p0 = gas.P
    r0 = gas.density

    # Perturb the pressure
    p1 = p0 * 1.0001

    # Set a new state w/ the same entropy and composition but adjusted pressure
    gas.SP = s0, p1

    # Method #1 to calculate "frozen" speed of sound
    afrozen = math.sqrt((p1 - p0)/(gas.density - r0))

    # Equilibrate the gas while holding entropy & pressure constant
    gas.equilibrate('SP', rtol=rtol, maxiter=maxiter)

    # Calcluate the equilibrium speed of sound
    aequil = math.sqrt((p1 - p0)/(gas.density - r0))

    # Method #2 to calculate "frozen" a via the ideal gas expression
    gamma = gas.cp / gas.cv
    mw = gas.mean_molecular_weight
    afrozen2 = math.sqrt(gamma * ct.gas_constant * gas.T / mw)

    # Set pressure back to initial value (not needed if perturbation is small)
    # gas.SP = s0, p0
    # gas.equilibrate('SP', rtol=rtol, maxiter=maxiter)

    return aequil, afrozen, afrozen2
