#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Isentropic Flow Relations
"""

import numpy as np

import isentropic
import normal_shock


class Nozzle(object):
    """Represents a Nozzle with Quasi 1-D flow
    """

    def __init__(self):
        pass

    def mass_flow_rate(rho, u, A):
        """Calculate flow through the nozzle given stagnation
        properties and back pressure

        :param <float> rho: Density (kg/m^3)
        :param <float> u: Flow speed (m/s)
        :param <float> A: Nozzle cross-sectional area (m^2)

        :return <float> mdot (kg/s)
        """

        return rho * u * A

    def solve_flow_given_stagnation(self, p0, T0, p_back, Ae_Astar, gamma):
        """Calculate flow through the nozzle given stagnation
        properties and back pressure

        :param <float> p0: Stagnation pressure (Pa)
        :param <float> T0: Stagnation temperature (K)
        :param <float> p_back: Back pressure at the exit (Pa)
        :param <float> Ae_Astar: Exit area ratio
        :param <float> gamma: Specific heat ratio

        :return <float> ...
        """

        p_ratio = p_back / p0

        if p_ratio > 1.0:
            raise ValueError("Back pressure must be less than stagnation")

        # Calculate exit mach #
        Me1, Me2 = isentropic.mach_from_area_ratio(Ae_Astar, gamma)

        # Subsonic case
        pe_p0_1 = 1 / isentropic.p0_p(Me1, gamma)

        # Supersonic case
        pe_p0_2 = 1 / isentropic.p0_p(Me2, gamma)

        # Shock at the nozzle exit case
        pe_p0_3 = pe_p0_2 * normal_shock.p2_p1(Me2, gamma)

        # Check which case matches back pressure
        M = None
        if p_ratio >= pe_p0_1:
            # Fully subsonic flow
            M = Me1

        elif p_ratio >= pe_p0_3:
            # Shock inside the nozzle
            # TODO: Compute shock location & properties
            pass

        else:
            # Fully isentropic subsonic - supersonic flow
            M = Me2

        return M

    def visualize_flow(self, Ae_Astar, gamma):
        """Visualize property variations through the nozzle

        :param <float> Ae_Astar: Exit area ratio
        :param <float> gamma: Specific heat ratio

        """

        # Create area array from throat to exit
        A_arr = np.linspace(1, Ae_Astar, num=100)

        M_arr = np.empty_like(A_arr)

        # Calculate Mach # at each area ratio
        for idx, A_Astar in enumerate(A_arr):
            M_arr[idx] = isentropic.mach_from_area_ratio(A_Astar, gamma)

        p_p0_arr = 1 / isentropic.p0_p(M_arr, gamma)
        T_T0_arr = 1 / isentropic.T0_T(M_arr, gamma)

        # rho_arr = ideal_density(p_arr, R/MW, T_arr)
        # c_arr = ideal_sound(p_arr, rho_arr, gamma)
        # u_arr = M_arr * c_arr

        return M_arr, p_p0_arr, T_T0_arr


if __name__ == "__main__":
    pass
