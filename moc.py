# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 17:32:36 2018

@author: asikora
"""

import math
import numpy as np

#from fluid import Fluid
#
#MEDIUM = Fluid(['Air'])
#
#def calculate_stagnation_h_s(T0, P0):
#    """
#    Calculate stagnation enthalpy and entropy from stagnation temperature and pressure
#    """
#    h0 = MEDIUM.enthalpy(T0, P0)
#    s0 = MEDIUM.entropy(T0, P0)
#    
#    return h0, s0

#
#def speed_of_sound_from_h_s(h, s):
#    """
#    Calculate speed of sound from enthalpy and entropy
#    """
#    
#    return MEDIUM.custom('A', 'H', h, 'S', s)

import cantera as ct


def nu(m, gamma):
    """
    compute prandtl meyer angles given by input mach numbers
    
    input:
        m - mach numbers
        gamma - specific heat ratio
    
    output:
        n - prandtl meyer angles (radians)
    """
    return np.sqrt((gamma + 1) / (gamma - 1)) * np.arctan(np.sqrt((gamma - 1) / (gamma + 1) * (m**2 - 1))) - np.arctan(np.sqrt(m**2 - 1))


def m_nu(n, gamma):
    """
    compute mach numbers given by set of prandtl meyer angles
    
    input:
        n - pm angles (radians)
        gamma - specific heat ratio
    
    output:
        m - mach numbers
    """
    if isinstance(n, (int, float)):
        n = np.array([n])

    m = np.zeros(n.shape)
    numax = (np.sqrt((gamma + 1) / (gamma - 1)) - 1) * np.pi / 2

    for i in range(len(n)):
        if n[i] < 0 or n[i] > numax:
            m[i] = np.nan
        elif n[i] == 0:
            m[i] = 1
        else:
            mnew = 2
            m[i] = 0
            while (abs(mnew - m[i]) > 0.00001):
                m[i] = mnew
                fm = nu(m[i], gamma) - n[i]
                fdm = np.sqrt(m[i]**2 - 1) / (1 + 0.5 * (gamma - 1) * m[i] ** 2) / m[i]
                mnew = m[i] - fm / fdm
            m[i] = mnew

    return m

    
def degrees_to_radians(deg):
    return deg * math.pi / 180

    
def prandtl_meyer(M, dtheta):
    """Calculate the resultant Mach number from input Mach number and expansion fan angle
    """
    
    gamma = 1.4
    
    nu_M = nu(M, gamma)
    
    nu_M2 = nu_M + dtheta
    
    M2 = m_nu(nu_M2, gamma)
    
    return M2[0]
    

def enthalpy_from_h_stag_u(h0, u):
    """Calculate enthalpy from stagnation enthalpy and velocity
    """
    
    return h0 - 0.5 * (u ** 2)
    

def prandtl_meyer_real_gas(gas, u1, h0, s0, dtheta):
    """
    Calculate velocity and speed of sound behind an expansion fan of angle dtheta
    """
    
    # Calculate enthalpy, entropy and speed of sound before the expansion
    h1 = enthalpy_from_h_stag_u(h0, u1)
    s1 = s0  # isentropic

    # print h1, s1
    a1 = equilSoundSpeeds(gas, h1)
    print a1
    print 'Mach Number: {}'.format(u1 / a1)
    # Guess velocity after the expansion based on perfect gas result
    u2 = 1300 # prandtl_meyer(u1 / a1, dtheta) * a1
    print u2

    # Calculate enthalpy and entropy based off this velocity
    h2 = enthalpy_from_h_stag_u(h0, u2)
    # print h2
    s2 = s0

    a2 = equilSoundSpeeds(gas, h2)
    # a2 = speed_of_sound_from_h_s(h2, s2)

    # print a2
    # return
    
    # return
    def I(v, a):
        v = float(v)
        a = float(a)
        
        return ((1 / v) * math.sqrt(v ** 2 / a ** 2 - 1))    

    N = 6

    dtheta_numerical = 0.5 * I(u1, a1)    
    for i in range(N):
        ui = u1 + i * (u2 - u1) / (N + 1)
        hi = h1 + 0.5 * (u1 ** 2) - 0.5 * (ui ** 2)
        si = s1
        # ai = speed_of_sound_from_h_s(hi, si)
        ai = equilSoundSpeeds(gas, hi)

        dtheta_numerical += (I(ui, ai) - I(u2, a2))

    print 'Numerical: {}'.format(dtheta_numerical)
    print 'Error: {}'.format(dtheta - dtheta_numerical)
    
    return


def equilSoundSpeeds(gas, h, rtol=1.0e-6, maxiter=5000):
    """
    www.cantera.org/docs/sphinx/html/cython/examples/thermo_sound_speed.html
    """
    
    gas.HP = h, gas.P
    gas.equilibrate('TP', rtol=rtol, maxiter=maxiter)
    
    s0 = gas.s
    p0 = gas.P
    r0 = gas.density
    
    p1 = p0 * 1.0001
    
    gas.SP = s0, p1
    gas.equilibrate('SP', rtol=rtol, maxiter=maxiter)
    
    a_equil = math.sqrt((p1 - p0) / (gas.density - r0))
    
    # gamma = gas.cp / gas.cv
    # a_frozen = math.sqrt(gamma * ct.gas_constant * gas.T / gas.mean_molecular_weight)
    
    # print a_equil
    # print a_frozen

    return a_equil
    
    
if __name__ == "__main__":
    pass

    # Combustion Chamber Properties
    Tc = 900  # K
    Pc = 1.013e7  # Pa

    AIR = ct.Solution('air.cti')
    AIR.TP = Tc, Pc

    h0 = AIR.h
    s0 = AIR.s
    # h0, s0 = calculate_stagnation_h_s(Tc, Pc)

    prandtl_meyer_real_gas(AIR, 700, h0, s0, degrees_to_radians(20))
