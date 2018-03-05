# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 17:32:36 2018

@author: asikora
"""

import math
import numpy as np

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


def enthalpy_from_h_stag_u(h0, u):
    """Calculate enthalpy from stagnation enthalpy and velocity
    """
    
    return h0 - 0.5 * (u ** 2)


def equilSoundSpeeds(gas, h=None, rtol=1.0e-6, maxiter=5000):
    """
    www.cantera.org/docs/sphinx/html/cython/examples/thermo_sound_speed.html
    """
    
    if h:
        gas.HP = h, gas.P
        gas.equilibrate('HP', rtol=rtol, maxiter=maxiter)
    else:
        gas.equilibrate('TP', rtol=rtol, maxiter=maxiter)

    s0 = gas.s
    p0 = gas.P
    r0 = gas.density
    
    p1 = p0 * 1.0001
    
    gas.SP = s0, p1
    gas.equilibrate('SP', rtol=rtol, maxiter=maxiter)
    
    a_equil = math.sqrt((p1 - p0) / (gas.density - r0))
    
#    gamma = gas.cp / gas.cv
#    a_frozen = math.sqrt(gamma * ct.gas_constant * gas.T / gas.mean_molecular_weight)
    
    # print a_equil
    # print a_frozen

    return a_equil# , a_frozen
#    return a_frozen


def eta(gas):
    
    gamma = gas.cp / gas.cv
    eta_frozen = (1.0 + gamma) / (gamma - 1.0)
    
    h = gas.h
    a, _ = equilSoundSpeeds(gas)
    eta_equil = 1 + 2 * h / (a ** 2)
    
    return eta_equil, eta_frozen
    

#def main(gas, h, M0, dtheta):
#    
#    # Get speed of sound
#    a, _ = equilSoundSpeeds(gas, h)
#    
#    u = M0 * a
#
#    eta = 1 + 2 * h / (a ** 2)
#    v_theta = a
#    v_r = math.sqrt(u ** 2 - a ** 2)
#    c = math.sqrt(2 * (h + (v_theta ** 2 + v_r ** 2)/2))
#    
#    phi = math.acos(a * math.sqrt(eta) / c)
#    # print(a, eta, c)
#    phi = 1
#    
#    M = math.sqrt(1 + eta * math.tan(phi) ** 2)
#    
#    nu = dtheta - math.acos(1 / M)
#    
#    print(M, nu)
    

def integrand(F, h0):
    """
    """

    h = h0 - F
    a = equilSoundSpeeds(AIR, h)
    
#    print(a)
#    print(F)
#    print('----')
    
    try:
        fn = (1.0 / F) * math.sqrt((2 * F / a) - 1)
    except ValueError:
        return np.inf
        
    return fn

from scipy.integrate import quad
from scipy.optimize import fsolve


def func(F, Fa, h0, dtheta):

    y, err = quad(integrand, Fa, F, args=(h0))
    
    print(y  -dtheta)
    
    return y - dtheta


def main(ua, h0, dtheta):
    
    dtheta = degrees_to_radians(dtheta)
    
    Fa = 0.5 * ua ** 2
    
    sol = fsolve(func, 0, args=(Fa, h0, dtheta))
    
    sol = sol[0]
    
    u = math.sqrt(2 * sol)
    h = h0 - 0.5 * u ** 2
    
    a = equilSoundSpeeds(AIR, h)

    print(u/a)

    # print('Solution: {}'.format(sol[0]))
    
    
if __name__ == "__main__":
    pass

    # Combustion Chamber Properties
    Tc = 6140  # K
    Pc = 1.2 * 101325  # Pa

    AIR = ct.Solution('air.cti')
    AIR.TP = Tc, Pc
    AIR.equilibrate('TP')
    
    main(1700, AIR.h, 20)

    # sol = fsolve(func, 3.0, args=())
    # print('Solution: {}'.format(sol[0]))
