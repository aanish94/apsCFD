#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Equilibrium Quasi 1-D Nozzle Flow
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

from gas import initialize_gas_object, speed_of_sound


def nozzle(gas, T0, p0, A_Astar):
    """Equilibrium normal shock wave flow - calculate post-shock
    conditions (eq. 17.4)

    :param <cantera.composite.Solution> gas: Solution object
    :param <float> T0: Stagnation Temperature (K)
    :param <float> p0: Stagnation Pressure (Pa)
    :param <float> A_Astar: Area ratio to throat

    :return <float> u2: Velocity post-shock (m/s)
    :return <float> T2: Temperature post-shock (K)
    """

    # Bring gas to equilibrium at specified temperature and pressure
    gas.TP = T0, p0
    gas.equilibrate('TP')

    # Extract pre-shock properties
    # rho = gas.density
    h0 = gas.h
    a0, _, _ = speed_of_sound(gas)
    
    
def calculate_mole_mass_ratios(gas):
    """Calculate mole-mass ratio
    """

    species = gas.species()
    species = [x.name for x in species]

    mole_fractions = gas.mole_fraction_dict()
    molecular_weights = gas.molecular_weights

    M = 0
    mole_mass_ratios = {}    
    for idx, spec in enumerate(species):
        print(spec)

        cur_x = mole_fractions[spec]
        cur_weight = molecular_weights[idx]
        print (cur_weight)
        
        M += cur_x * cur_weight
        
        mole_mass_ratios[spec] = cur_x
        
    for key, value in mole_mass_ratios.items():
        mole_mass_ratios[key] = value / M
        
    return mole_mass_ratios

# TODO: Solve 17.15, 17.16 and 17.18 numerically


if __name__ == "__main__":
    air = initialize_gas_object('air')

    p0 = 100 * 101325
    T0 = 8000
    
    A_Astars = np.linspace(1, 2)
    
    air.TP = T0, p0
    air.equilibrate('TP')
    
    h0 = air.h
    rho0 = air.density

    def solve_for_p(p_in, h0):

        air.SP = air.s, p_in
        air.equilibrate('SP')

        u_calc = np.sqrt(2 * (h0 - air.h))
        a_calc, _, _ = speed_of_sound(air)

        print(u_calc - a_calc)
        return u_calc - a_calc
    
    pstar = fsolve(solve_for_p, p0 / 1.8888, args=(h0), xtol=1e-5)
    
    air.SP = air.s, pstar
    air.equilibrate('SP')

    rhostar = air.density
    ustar = np.sqrt(2 * (h0 - air.h))
    
    for A_Astar in A_Astars:
        pass
    
        def solve_for_h(h, h0, rhostar, ustar, A_Astar):
            
            air.HP = h, air.P
            air.equilibrate('HP')
            
            rho_cur = air.density
            
            rhs = rhostar * ustar / A_Astar
            lhs = rho_cur * np.sqrt(2 * (h0 - h))
            
            return lhs - rhs
    
        # h_cur = fsolve(solve_for_h, h0 * 0.5, args=(h0, rhostar, ustar, A_Astar))
