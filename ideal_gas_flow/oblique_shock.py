#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Oblique Shock Relations

Note:
    - Beta is the shock angle (radians) w.r.t incoming stream direction
    - Theta is the flow deflection (radians) w.r.t incoming stream
"""

import math
from scipy import optimize


def beta(M, theta, gamma):
    """Oblique shock wave angle

    :param <float> M: Mach # upstream
    :param <float> theta: Flow deflection angle (radians)
    :param <float> gamma: Specific heat ratio

    :return <float> Shock angle w.r.t initial flow direction (radians)
    """

    # Initial guess
    b1 = math.asin(1.0 / M)

    def f_to_solve(beta):
        return theta(M, beta, gamma) - theta

    return optimize.newton(f_to_solve, x0=b1)


def mach(M, beta, theta, gamma):
    """Mach # after an oblique shock (eq. 4.12)

    :param <float> M: Mach # upstream
    :param <float> Beta: Shock angle w.r.t initial flow direction (radians)
    :param <float> theta: Flow deflection angle (radians)
    :param <float> gamma: Specific heat ratio

    :return <float> Mach # in flow after the shock
    """

    m1sb = M * math.sin(beta)
    
    n1 = 1.0 + (gamma - 1.0) * 0.5 * m1sb ** 2
    d1 = gamma * m1sb ** 2 - (gamma - 1.0) * 0.5

    m2 = math.sqrt(n1 / d1 / (math.sin(beta - theta)) ** 2)

    return m2


def p2_p1(M, beta, gamma):
    """Pressure ratio across an olique shock (eq. 4.9)

    :param <float> M: Mach # upstream
    :param <float> Beta: Shock angle w.r.t initial flow direction (radians)
    :param <float> gamma: Specific heat ratio

    :return <float> Pressure ratio p2/p1
    """

    m1sb = M * math.sin(beta)

    return 1.0 + 2.0 * gamma / (gamma + 1.0) * (m1sb ** 2 - 1.0)


def p02_p01(M, beta, gamma):
    """Stagnation pressure ratio across an olique shock (eq. 4.9)

    :param <float> M: Mach # upstream
    :param <float> Beta: Shock angle w.r.t initial flow direction (radians)
    :param <float> gamma: Specific heat ratio

    :return <float> Total pressure ratio p02/p01
    """

    m1sb = M * math.sin(beta)

    t1 = (gamma + 1.0) / (2.0 * gamma * m1sb ** 2 - (gamma - 1.0))
    t2 = (gamma + 1.) * m1sb ** 2 / (2.0 + (gamma - 1.0) * m1sb ** 2)

    return t1 ** (1.0 / (gamma - 1.0)) * t2 ** (gamma / (gamma - 1.0))


def r2_r1(M, beta, gamma):
    """Density ratio across an olique shock (eq. 4.8)

    :param <float> M: Mach # upstream
    :param <float> Beta: Shock angle w.r.t initial flow direction (radians)
    :param <float> gamma: Specific heat ratio

    :return <float> Density ratio r2/r1
    """

    m1sb = M * math.sin(beta)

    n1 = (gamma + 1.0) * m1sb ** 2
    d1 = 2.0 + (gamma - 1.0) * m1sb ** 2

    return n1 / d1


def theta(M, beta, gamma):
    """Flow deflection angle (eq. 4.17)

    :param <float> M: Mach # upstream
    :param <float> Beta: Shock angle w.r.t initial flow direction (radians)
    :param <float> gamma: Specific heat ratio

    :return <float> Flow deflection angle (radians)
    """

    m1sb = M * math.sin(beta)

    t1 = 2.0 / math.tan(beta) * (m1sb ** 2 - 1.0)
    t2 = M ** 2 * (gamma + math.cos(2.0 * beta)) + 2.0
    
    return math.atan(t1 / t2)


def t2_t1(M, beta, gamma):
    """Temperature ratio across an olique shock (eq. 4.11)

    :param <float> M: Mach # upstream
    :param <float> Beta: Shock angle w.r.t initial flow direction (radians)
    :param <float> gamma: Specific heat ratio

    :return <float> Pressure ratio T2/T1
    """

    return p2_p1(M, beta, gamma) / r2_r1(M, beta, gamma)


def u2_u1(M, beta, gamma):
    """Flow speed ratio across an olique shock (eq. 4.8)

    :param <float> M: Mach # upstream
    :param <float> Beta: Shock angle w.r.t initial flow direction (radians)
    :param <float> gamma: Specific heat ratio

    :return <float> Velocity ratio u2/u1
    """

    t1 = math.sqrt((math.sin(beta) / r2_r1(M, beta, gamma)) ** 2 + (math.cos(beta)) ** 2)

    return t1
