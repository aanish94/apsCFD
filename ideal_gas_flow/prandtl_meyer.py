#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Prandtl-Meyer Expansion Waves
"""

import math
from scipy import optimize


def nu(M , gamma):
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


def mach_from_nu(nu , gamma):
    """Invert the Prandtl-Meyer function (eq. 4.44)

    :param <float> Prandtl-Meyer function (radians)
    :param <float> gamma: Specific heat ratio

    :return <float> M: Mach #
    """

    def f_to_solve(M):
        return prandtl_meyer(M, gamma)

    return optimize.newton(f_to_solve, x0=2.0)
