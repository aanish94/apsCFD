#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
1-D and 2-D steady flow of an ideal gas

    - Isentropic flow relations
    - 1D (Normal) shock relations
    - 1D flow with heat addition (Rayleigh)
    - 1D flow with friction (Fanno)
    - Prandtl-Meyer functions
    - Oblique shock relations
"""

import math

"""Isentropic Flow
"""


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

"""1-D Normal Shock Relations
"""


def shock_density(M1, gamma):
    """Density ratio across a normal shock (eq. 3.53)

    :param <float> M1: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Density ratio r2/r1
    """

    n1 = (g + 1.0) * M1 ** 2
    d1 = 2.0 + (g - 1.0) * M1 ** 2

    return n1 / d1


def shock_mach(M1, gamma):
    """Mach # after a normal shock (eq. 3.51)

    :param <float> M1: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Mach # after the chock
    """

    n1 = 1.0 + (gamma - 1.0) * 0.5 * M1 ** 2
    d1 = gamma * M1 ** 2 - (gamma - 1.0) * 0.5

    return math.sqrt(n1 / d1)


def shock_pressure(M1, gamma):
    """Pressure ratio across a normal shock (eq. 3.57)

    :param <float> M1: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Pressure ratio p2/p1
    """

    return 1.0 + 2.0 * gamma / (gamma + 1.0) * (M1 ** 2 - 1.0)


def shock_stag_pressure(M1, gamma):
    """Stagnation pressure ratio across a normal shock

    :param <float> M1: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Stagnation pressure ratio p02/p01
    """

    t1 = (gamma + 1.0) / (2.0 * gamma * M1 ** 2 - (gamma - 1.0))
    t2 = (gamma + 1.0) * M1 ** 2 / (2.0 + (gamma - 1.0) * M1 ** 2)

    return t1 ** (1.0 / (gamma - 1.0)) * t2 ** (gamma / (gamma - 1.0))


def shock_temperature(M1, gamma):
    """Temperature ratio across a normal shock (eq. 3.59)

    :param <float> M1: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Temperature ratio T2/T1
    """

    return shock_pressure(M1, gamma) / shock_density(M1, gamma)


def shock_velocity(M1, gamma):
    """Velocity ratio across a normal shock (eq. 3.53)

    :param <float> M1: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Velocity ratio u2/u1
    """

    return 1 / shock_rho(M1, gamma)


"""Rayleigh Flow (1-D flow w/ heat addition)

Note:
    - Star denotes conditions achieved if sufficient heat was added to achieve sonic conditions.
"""


def heat_p_pstar(M, gamma):
    """Static pressure ratio for flow with heat addition (eq. 3.85)

    :param <float> M: Initial Mach #
    :param <float> gamma: Specific heat ratio

    :return <float> Static pressure ratio P/Pstar
    """

    return (1.0 + gamma) / (1.0 + gamma * M ** 2)


def heat_p0_p0star(M, gamma):
    """Stagnation pressure ratio for flow with heat addition (eq. 3.88)

    :param <float> M: Initial Mach #
    :param <float> gamma: Specific heat ratio

    :return <float> Stagnation pressure ratio P0/P0star
    """

    t1 = (2.0 + (gamma - 1.0) * M ** 2) / (gamma + 1.0)
    t2 = gamma / (gamma - 1.0)

    return (1.0 + gamma) / (1.0 + gamma * M ** 2) * t1 * t2


def heat_r_rstar(M, gamma):
    """Static density ratio for flow with heat addition (eq. 3.87)

    :param <float> M: Initial Mach #
    :param <float> gamma: Specific heat ratio

    :return <float> Static density ratio r/rstar
    """

    return 1.0 / M ** 2 / (1.0 + gamma) * (1.0 + gamma * M ** 2)


def heat_t_tstar(M, gamma):
    """Static temperature ratio for flow with heat addition (eq. 3.86)

    :param <float> M: Initial Mach #
    :param <float> gamma: Specific heat ratio

    :return <float> Static temperature ratio T/Tstar
    """

    return M ** 2  * ((1.0 + gamma) / (1.0 + gamma * M ** 2)) ** 2


def heat_t0_t0star(M, gamma):
    """Total temperature ratio for flow with heat addition (eq. 3.89)

    :param <float> M: Initial Mach #
    :param <float> gamma: Specific heat ratio

    :return <float> Total temperature ratio T0/T0star
    """

    t1 = (gamma + 1) * M ** 2
    t2 = (1.0 + gamma * M ** 2) ** 2
    t3 = 2.0 + (gamma - 1.0) * M ** 2

    return t1 / t2 * t3


"""Prandtl-Meyer Expansion Waves
"""

def prandtl_meyer(M , gamma):
    """Prandtl-Meyer function (eq. 4.44)

    :param <float> M: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Prandtl-Meyer function (radians)
    """

    # Supersonic
    if M > 1.0:
        t1 = M ** 2 - 1.0
        t2 = math.sqrt(( gamma - 1.0) / ( gamma + 1.0) * t1 )
        t3 = math.sqrt( t1 )
        t4 = math.sqrt(( gamma + 1.0) / ( gamma - 1.0) )
        nu = t4 * math.atan( t2 ) - math.atan( t3 )
    # Subsonic
    else :
        nu = 0.0

    return nu


def prandtl_meyer(M , gamma):
    """Invert the Prandtl-Meyer function (eq. 4.44)

    :param <float> Prandtl-Meyer function (radians)
    :param <float> gamma: Specific heat ratio

    :return <float> M: Mach #
    """

    pass


"""Oblique shock relations
"""

