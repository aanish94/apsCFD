#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Isentropic Flow Relations
"""

import math
from scipy import optimize


def A_Astar(M, gamma):
    """Area-Mach # relation for isentropic, quasi-1D flow (eq. 5.20)

    :param <float> M: Mach # at area A
    :param <float> gamma: Specific heat ratio

    :return <float> Area ratio A:A_star
    """

    t1 = (gamma + 1.0) / (gamma - 1.0)
    m2 = M ** 2
    t2 = 1.0 / m2 * (2.0 / (gamma + 1.0) * (1.0 +(gamma - 1.0) * 0.5 * m2)) ** t1

    return math.sqrt(t2)


def mach_angle(M):
    """Calculate the Mach angle

    :param <float> M: Mach #

    :return <float> Mach angle, mu
    """

    return math.asin(1 / M)


def mach_from_area_ratio(A_Astar_in, gamma):
    """Invert the area-Mach # relation for isentropic, quasi-1D flow (eq. 5.20)

    :param <float> A_Astar_in: Area ratio A:A_star
    :param <float> gamma: Specific heat ratio

    :return <float> M_sub: Subsonic Mach # at area A
    :return <float> M_super: Supersonic Mach # at area A    
    """

    def f_to_solve(M):
        return A_Astar(M, gamma) - A_Astar_in

    M_sub = optimize.bisect(f_to_solve, 0.01, 1.0)
    M_super = optimize.newton(f_to_solve, x0=2.0)

    return M_sub, M_super


def p0_p(M, gamma):
    """Ratio of total to static pressure for isentropic flow (eq. 3.30)

    :param <float> M: Mach # at area A
    :param <float> gamma: Specific heat ratio

    :return <float> Pressure ratio p0/p
    """

    return (T0_T(M, gamma)) ** (gamma / (gamma - 1.0))


def r0_r(M, gamma):
    """Ratio of stagnation to free-stream density for isentropic flow (eq. 3.31)

    :param <float> M: Mach # at area A
    :param <float> gamma: Specific heat ratio

    :return <float> Density ratio r0/r
    """

    return (T0_T(M, gamma)) ** (1.0 / (gamma - 1.0))


def T0_T(M, gamma):
    """Ratio of total to static temperature for adiabatic flow (eq. 3.28)

    :param <float> M: Mach # at area A
    :param <float> gamma: Specific heat ratio

    :return <float> Temperature ratio T0/T
    """

    return 1.0 + 0.5 * (gamma - 1.0) * M ** 2
