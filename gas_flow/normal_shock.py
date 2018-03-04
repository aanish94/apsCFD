#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Normal Shock Relations
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

from gas import initialize_gas_object, speed_of_sound


def conservation_of_energy(h1, u1, rho1, rho2):
    """Solve conservation of energy for enthalpy

    :param <float> h1: Enthalpy state 1 (J/kg/K)
    :param <float> u1: Velocity state 1 (m/s)
    :param <float> rho1: Density state 1 (kg/m^3)
    :param <float> rho2: Density state 2 (kg/m^3)

    :return <float> h2: Enthalpy state 2 (Pa)
    """

    return h1 + 0.5 * u1 ** 2 * (1 - (rho1 / rho2) ** 2)


def conservation_of_momentum(p1, u1, rho1, rho2):
    """Solve conservation of momentum for pressure

    :param <float> p1: Pressure state 1 (Pa)
    :param <float> u1: Velocity state 1 (m/s)
    :param <float> rho1: Density state 1 (kg/m^3)
    :param <float> rho2: Density state 2 (kg/m^3)

    :return <float> p2: Pressure state 2 (Pa)    
    """

    return p1 + rho1 * u1 ** 2 * (1 - rho1 / rho2)


def normal_shock(gas, u1, T1, p1):
    """Equilibrium normal shock wave flow - calculate post-shock conditions (eq. 17.4)

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

    # print ("Initial Mach #: {}".format(u1 / a1))

    def solve_for_rho2(rho):
        """Helper function - solved iteratively for rho2
        """
        # Calculate pressure and enthalpy from conservation laws
        p = conservation_of_momentum(p1, u1, rho1, rho)
        h = conservation_of_energy(h1, u1, rho1, rho)
        
        # Bring gas to equilibrium at specified enthalpy and pressure
        air.HP = h, p
        air.equilibrate('HP')

        # Determine new density and return delta
        rho2 = air.density
        return rho - rho2
    
    # Solve iteratively for rho2
    rho2_guess = 10 * rho1  # assue an initial value
    rho2 = fsolve(solve_for_rho2, x0=rho2_guess)[0]

    # Calculate final properties
    T2 = air.T
    p2 = conservation_of_momentum(p1, u1, rho1, rho2)
    u2 = (rho1 / rho2) * u1
    a2, _, _ = speed_of_sound(air)

    # print ("Final Mach #: {}".format(u2 / a2))
    
    return u2, T2, p2


if __name__ == "__main__":
    air = initialize_gas_object('airNASA9')

    p1 = 0.0001 * 101325
    T1 = 225.0

    u1 = np.linspace(1000, 14000)
    T2_T1 = np.zeros(u1.shape)
    
    for idx, u in enumerate(u1):
        air = initialize_gas_object('airNASA9')

        _, T2, _ = normal_shock(air, u, T1, p1)
        T2_T1[idx] = T2 / T1        

    # plt.figure()
    plt.plot(u1, T2_T1)
