#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Rayleigh Flow (1-D flow w/ heat addition)

Note:
    - Star denotes conditions achieved if sufficient heat was added to achieve sonic conditions.
"""

import math
from scipy import optimize


def mach(T0T0star, gamma):
    """Static pressure ratio for flow with heat addition (eq. 3.85)

    :param <float> T0T0star: Total temperature ratio T0/T0star
    :param <float> gamma: Specific heat ratio

    :return <float> M: Initial Mach #
    """

    def f_to_solve(M):
        return T0_T0star(M, gamma) - T0T0star

    return optimize.newton(f_to_solve, 2.5)


def p_pstar(M, gamma):
    """Static pressure ratio for flow with heat addition (eq. 3.85)

    :param <float> M: Initial Mach #
    :param <float> gamma: Specific heat ratio

    :return <float> Static pressure ratio P/Pstar
    """

    return (1.0 + gamma) / (1.0 + gamma * M ** 2)


def p0_p0star(M, gamma):
    """Stagnation pressure ratio for flow with heat addition (eq. 3.88)

    :param <float> M: Initial Mach #
    :param <float> gamma: Specific heat ratio

    :return <float> Stagnation pressure ratio P0/P0star
    """

    t1 = (2.0 + (gamma - 1.0) * M ** 2) / (gamma + 1.0)
    t2 = gamma / (gamma - 1.0)

    return (1.0 + gamma) / (1.0 + gamma * M ** 2) * t1 * t2


def rho_rhostar(M, gamma):
    """Static density ratio for flow with heat addition (eq. 3.87)

    :param <float> M: Initial Mach #
    :param <float> gamma: Specific heat ratio

    :return <float> Static density ratio r/rstar
    """

    return 1.0 / M ** 2 / (1.0 + gamma) * (1.0 + gamma * M ** 2)


def T_Tstar(M, gamma):
    """Static temperature ratio for flow with heat addition (eq. 3.86)

    :param <float> M: Initial Mach #
    :param <float> gamma: Specific heat ratio

    :return <float> Static temperature ratio T/Tstar
    """

    return M ** 2  * ((1.0 + gamma) / (1.0 + gamma * M ** 2)) ** 2


def T0_T0star(M, gamma):
    """Total temperature ratio for flow with heat addition (eq. 3.89)

    :param <float> M: Initial Mach #
    :param <float> gamma: Specific heat ratio

    :return <float> Total temperature ratio T0/T0star
    """

    t1 = (gamma + 1) * M ** 2
    t2 = (1.0 + gamma * M ** 2) ** 2
    t3 = 2.0 + (gamma - 1.0) * M ** 2

    return t1 / t2 * t3
