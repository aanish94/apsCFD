# -*- coding: utf-8 -*-
"""
Created on Thu Mar 08 14:42:14 2018

@author: asikora
"""

from air_equilibrium_properties import a_from_p_s, h_from_p_rho, rho_from_p_s, T_from_p_rho

from gas import initialize_gas_object, speed_of_sound


T = 8000
P = 100 * 101325

air = initialize_gas_object('air')
air.TP = T, P
air.equilibrate('TP')

# print air.density
# print air.h
#
print T_from_p_rho(P, air.density) / T
print h_from_p_rho(P, air.density) / air.h

print speed_of_sound(air)


from scipy.optimize import fsolve, newton


def rho_from_h_p(h, p):
    def f_to_solve(rho, h, p):
        
        h1 = h_from_p_rho(p, rho)
        
        return h - h1
    

    rho = fsolve(f_to_solve, x0=1, args=(h, p))
    
    return rho

def rho_from_T_p(T, p):
    
    def f_to_solve(rho, T, p):
    
        if rho <= 0:
            return 1000 + abs(rho)
        T1 = T_from_p_rho(p, rho)
        
        return T - T1
    

    rho = fsolve(f_to_solve, x0=1, args=(T, p), xtol=1e-6)

    return rho
    
    
def s_from_p_rho(p, rho):

    def f_to_solve(s, p, rho):
    
        if s <= 0:
            return 1000 + abs(s)
        
        rho1 = rho_from_p_s(p, s)
        
        return rho - rho1
    
    s = fsolve(f_to_solve, x0=10000, args=(p, rho), xtol=1e-6)

    return s
    
    
    
Ttest = 225.0
Ptest = 22906.9277467

#print rho_from_T_p(Ttest, Ptest)
#
#
#invert_t_p_rho = lambda x: T_from_p_rho(Ptest, x) - Ttest
#
#r1 = newton(invert_t_p_rho, 1)
