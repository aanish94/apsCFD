#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Equilibrium Quasi 1-D Nozzle Flow
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

from gas import initialize_gas_object, speed_of_sound


def nozzle(gas, u1, T1, p1):
    """Equilibrium normal shock wave flow - calculate post-shock
    conditions (eq. 17.4)

    :param <cantera.composite.Solution> gas: Solution object
    :param <float> u1: Velocity pre-shock (m/s)
    :param <float> T1: Temperature pre-shock (K)
    :param <float> p1: Pressure pre-shock (Pa)

    :return <float> u2: Velocity post-shock (m/s)
    :return <float> T2: Temperature post-shock (K)
    """

    # Bring gas to equilibrium at specified temperature and pressure
    gas.TP = T1, p1
    gas.equilibrate('TP')

    # Extract pre-shock properties
    rho1 = gas.density
    h1 = gas.h
    a1, _, _ = speed_of_sound(gas)


# TODO: Solve 17.15, 17.16 and 17.18 numerically


if __name__ == "__main__":
    air = initialize_gas_object('airNASA9')
