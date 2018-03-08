#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Equilibrium Normal Shock Relations
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

from atmosphere import Atmosphere
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
    """Equilibrium normal shock wave flow - calculate post-shock
    conditions (eq. 17.4)

    :param <cantera.composite.Solution> gas: Solution object
    :param <float> u1: Velocity pre-shock (m/s)
    :param <float> T1: Temperature pre-shock (K)
    :param <float> p1: Pressure pre-shock (Pa)

    :return <float> u2: Velocity post-shock (m/s)
    :return <float> rho2: Density post-shock (kg/m^3)
    :return <float> T2: Temperature post-shock (K)
    :return <float> p2: Pressure post-shock (Pa)
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

    return u2, rho2 / rho1, T2, p2


if __name__ == "__main__":
    pass
    # Unit conversions
    atm_to_pa = 101325
    ft_to_m = 0.3048

    # Gas model
    air = initialize_gas_object('airNASA9')

    """Recreate Fig. 14.3 from Anderson Hypersonic & High Temperature Gas Dynamics
    """
    exact_x = []
    exact_y = []

    fig0, ax0 = plt.subplots(figsize=(12, 9))
    T1 = 225.0

    pressures = np.array([1e-1, 1e-2, 1e-4]) * atm_to_pa
    u1 = np.linspace(1000, 14000)

    # Calculate shock properties for each pressure
    for p1 in pressures:
    
        T2_T1 = np.zeros(u1.shape)

        for idx, u in enumerate(u1):

            _, _, T2, _ = normal_shock(air, u, T1, p1)
            T2_T1[idx] = T2 / T1

        ax0.plot(u1, T2_T1, label='Pressure: {} atm'.format(p1 / atm_to_pa))

    # Plot exact values
    for idx, x in enumerate(exact_x):
        ax0.plot(x, exact_y[idx], marker='x', color='r', markersize=10)

    # Plot formatting
    ax0.legend(loc='best')
    ax0.set_xlabel(r'$u_1 (\frac{km}{s})$', fontsize=18)
    ax0.set_ylabel(r'$\frac{T_2}{T_1}$', fontsize=18)
    ax0.set_title('Influence of Pressure on Normal Shock Temperature in Equilibrium Air', fontsize=16)

    fig0.tight_layout()
    fig0.savefig('results/topic_2b_figure_14_3.png', bbox_inches='tight')

    """Recreate Fig. 14.4 & 14.5 from Anderson Hypersonic & High Temperature Gas Dynamics
    """

    altitudes = np.array([35900, 59800, 82200, 100000, 120300, 154800, 173500,
                          200100, 230400, 259700, 294800, 322900]) * ft_to_m

    exact_x = []
    exact_y_1 = []
    exact_y_2 = []

    fig1, ax1 = plt.subplots(figsize=(12, 9))
    fig2, ax2 = plt.subplots(figsize=(12, 9))

    T1 = 225.0

    # Convert altitudes into ambient pressures
    pressures = np.zeros_like(altitudes)
    for idx, altitude in enumerate(altitudes):
        a = Atmosphere(altitude)
        pressures[idx] = a.P

    u1 = np.linspace(4000, 46000) * ft_to_m

    # Calculate shock properties for each pressure
    for pidx, p1 in enumerate(pressures):
        T2_arr = np.zeros_like(u1)
        rho2_rho1_arr = np.zeros_like(u1)

        for uidx, u in enumerate(u1):

            _, rho2_rho1, T2, _ = normal_shock(air, u, T1, p1)
            T2_arr[uidx] = T2
            rho2_rho1_arr[uidx] = rho2_rho1

        ax1.plot(u1 / ft_to_m, T2_arr, label='Altitude: {} ft'.format(altitudes[pidx] / ft_to_m))
        ax2.plot(u1 / ft_to_m, rho2_rho1_arr, label='Altitude: {} ft'.format(altitudes[pidx] / ft_to_m))

    # Plot exact values
    for idx, x in enumerate(exact_x):
        ax1.plot(x, exact_y_1[idx], marker='x', color='r', markersize=10)
        ax2.plot(x, exact_y_2[idx], marker='x', color='r', markersize=10)

    # Plot formatting
    ax1.legend(loc='best')
    ax1.set_xlabel(r'$u_1 (\frac{ft}{s})$', fontsize=18)
    ax1.set_ylabel(r'$T_2 (K)$', fontsize=18)
    ax1.set_title('Variation of Normal Shock Temperature with Velocity & Altitude', fontsize=16)
    ax1.set_xticks(np.arange(min(u1 / ft_to_m), max(u1 / ft_to_m) + 1, 2000))
    fig1.tight_layout()
    fig1.savefig('results/topic_2b_figure_14_4.png', bbox_inches='tight')

    # Plot formatting
    ax2.legend(loc='best')
    ax2.set_xlabel(r'$u_1 (\frac{ft}{s})$', fontsize=18)
    ax2.set_ylabel(r'$\frac{\rho_2}{\rho_1}$', fontsize=18)
    ax2.set_title('Variation of Normal Shock Density with Velocity & Altitude', fontsize=16)
    ax2.set_xticks(np.arange(min(u1 / ft_to_m), max(u1 / ft_to_m) + 1, 2000))    
    fig2.tight_layout()
    fig2.savefig('results/topic_2b_figure_14_5.png', bbox_inches='tight')
