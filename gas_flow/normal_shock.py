#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Normal Shock Relations
"""

import math
from scipy.optimize import fsolve

from gas import initialize_gas_object, speed_of_sound


air = initialize_gas_object('airNASA9')

# Given u1, p1, T1

u1 = 10000
p1 = 0.01 * 101325
T1 = 225
#
#air.TP = T1, p1
#air.equilibrate('TP')
#
#rho1 = air.density
#h1 = air.h
#a1, _, _ = speed_of_sound(air)
#
#print ("Initial Mach #: {}".format(u1 / a1))
#
## Assume a value for rho1 / rho2
#rho1_rho2 = 0.1
#
#for count in range(20):
#    rho1_rho2_old = rho1_rho2
#
#    # Calculate associated p2 and h2
#    p2 = p1 + rho1 * u1 ** 2 * (1 - rho1_rho2)
#    h2 = h1 + 0.5 * u1 ** 2 * (1 - rho1_rho2 ** 2)        
#    air.HP = h2, p2
#    air.equilibrate('HP')
#    
#    # Calculate new rho1 / rho2
#    rho1_rho2 = rho1 / air.density
#    
#    if (abs(rho1_rho2 - rho1_rho2_old)) < 1.0e-6:
#        break
#
#T2 = air.T
#u2 = rho1_rho2 * u1
#a2, _, _ = speed_of_sound(air)
#
#print ("Final Mach #: {}".format(u2 / a2))
#print ("T2/T1 = {}".format(T2 / T1))


def normal_shock(gas, u1, T1, p1):
    """Equilibrium normal shock wave flow - calculate post-shock conditions

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

    print ("Initial Mach #: {}".format(u1 / a1))
    
    def solve_for_rho2(rho):
        """Helper function - solved iteratively for rho2
        """

        # Calculate pressure and enthalpy
        p = p1 + rho1 * u1 ** 2 * (1 - rho1 / rho)
        h = h1 + 0.5 * u1 ** 2 * (1 - (rho1 / rho) ** 2)        
        
        # Bring gas to equilibrium at specified enthalpy and pressure
        air.HP = h, p
        air.equilibrate('HP')

        # Determine new density and return delta
        rho2 = air.density
        return rho - rho2
    
    # Solve iteratively for rho2
    rho2_guess = 10 * rho1
    rho2 = fsolve(solve_for_rho2, x0=rho2_guess)

    T2 = air.T
    u2 = (rho1 / rho2) * u1
    a2, _, _ = speed_of_sound(air)

    print ("Final Mach #: {}".format(u2 / a2))
    print ("T2/T1 = {}".format(T2 / T1))
    
    return u2, T2

normal_shock(air, u1, T1, p1)
