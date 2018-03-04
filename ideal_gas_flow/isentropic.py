#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Isentropic Flow Relations
"""

import math


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


def p0_T(M, gamma):
    """Ratio of total to static pressure for isentropic flow (eq. 3.30)

    :param <float> M: Mach # at area A
    :param <float> gamma: Specific heat ratio

    :return <float> Pressure ratio p0/p
    """

    return (T0_T(M, gamma)) ** (gamma / (gamma - 1.0))


def r0_T(M, gamma):
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
