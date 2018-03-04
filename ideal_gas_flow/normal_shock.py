#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Normal Shock Relations
"""

import math


def mach(M1, gamma):
    """Mach # after a normal shock (eq. 3.51)

    :param <float> M1: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Mach # after the chock
    """

    n1 = 1.0 + (gamma - 1.0) * 0.5 * M1 ** 2
    d1 = gamma * M1 ** 2 - (gamma - 1.0) * 0.5

    return math.sqrt(n1 / d1)


def p2_p1(M1, gamma):
    """Pressure ratio across a normal shock (eq. 3.57)

    :param <float> M1: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Pressure ratio p2/p1
    """

    return 1.0 + 2.0 * gamma / (gamma + 1.0) * (M1 ** 2 - 1.0)


def p02_p01(M1, gamma):
    """Stagnation pressure ratio across a normal shock

    :param <float> M1: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Stagnation pressure ratio p02/p01
    """

    t1 = (gamma + 1.0) / (2.0 * gamma * M1 ** 2 - (gamma - 1.0))
    t2 = (gamma + 1.0) * M1 ** 2 / (2.0 + (gamma - 1.0) * M1 ** 2)

    return t1 ** (1.0 / (gamma - 1.0)) * t2 ** (gamma / (gamma - 1.0))


def rho2_rho1(M1, gamma):
    """Density ratio across a normal shock (eq. 3.53)

    :param <float> M1: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Density ratio r2/r1
    """

    n1 = (gamma + 1.0) * M1 ** 2
    d1 = 2.0 + (gamma - 1.0) * M1 ** 2

    return n1 / d1


def T2_T1(M1, gamma):
    """Temperature ratio across a normal shock (eq. 3.59)

    :param <float> M1: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Temperature ratio T2/T1
    """

    return p2_p1(M1, gamma) / rho2_rho1(M1, gamma)


def u2_u1(M1, gamma):
    """Velocity ratio across a normal shock (eq. 3.53)

    :param <float> M1: Mach # before the shock
    :param <float> gamma: Specific heat ratio

    :return <float> Velocity ratio u2/u1
    """

    return 1 / rho2_rho1(M1, gamma)
