#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Equilibrium Quasi 1-D Nozzle Flow
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

from gas import initialize_gas_object

from air_equilibrium_properties import a_from_p_s, h_from_p_rho, rho_from_p_s, T_from_p_rho

import isentropic

# Unit conversions
atm_to_pa = 101325
ft_to_m = 0.3048
si_to_btu_lb = 1 / 2326.
k_to_rankine = 9.0/5.0


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
    # gas.TP = T0, p0
    # gas.equilibrate('TP')

    # Extract pre-shock properties
    # rho = gas.density
    # h0 = gas.h
    # a0, _, _ = speed_of_sound(gas)
    pass


def rho_from_T_p(T, p):

    def f_to_solve(rho, T, p):

        if rho <= 0:
            return 1000 + abs(rho)
        T1 = T_from_p_rho(p, rho)

        return T - T1

    rho = fsolve(f_to_solve, x0=1, args=(T, p), xtol=1e-6)

    return rho


def rho_from_h_p(h, p):

    def f_to_solve(rho, h, p):

        if rho <= 0:
            return 1000 + abs(rho)
        h1 = h_from_p_rho(p, rho)

        return h - h1

    rho = fsolve(f_to_solve, x0=1, args=(h, p), xtol=1e-6)

    return rho


def s_from_p_rho(p, rho):

    def f_to_solve(s, p, rho):

        if s <= 0:
            return 1000 + abs(s)

        rho1 = rho_from_p_s(p, s)

        return rho - rho1

    s = fsolve(f_to_solve, x0=10000, args=(p, rho), xtol=1e-6)

    return s


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

        cur_x = mole_fractions.get(spec, 0)
        cur_weight = molecular_weights[idx]
        print (cur_weight)

        M += cur_x * cur_weight

        mole_mass_ratios[spec] = cur_x

    for key, value in mole_mass_ratios.items():
        mole_mass_ratios[key] = value / M

    return mole_mass_ratios


if __name__ == "__main__":
    pass

#    fig6, ax6 = plt.subplots(figsize=(12, 9))
#    fig7, ax7 = plt.subplots(figsize=(12, 9))

    # all_h0s = [400, 600, 8000, 1000, 1200, 1400, 1600, 1800, 2000, 3000, 4000, 6000, 8000, 10000]
    all_h0s = [8000]

    # Area ratios
    A_Astars = np.linspace(1, 16)

    # Combustion Chamber properties
    p0 = 100 * 101325.0
    # T0 = 8000.0
    # Calculate stagnation properties
    # rho0 = rho_from_T_p(T0, p0)
    # h0 = h_from_p_rho(p0, rho0)
    # s0 = s_from_p_rho(p0, rho0)

    h0 = 8000
    for h0 in all_h0s:
        h0 = h0 / si_to_btu_lb
        rho0 = rho_from_h_p(h0, p0)
        T0 = T_from_p_rho(p0, rho0)
        s0 = s_from_p_rho(p0, rho0)

        # Calculate throat conditions
        def solve_for_p(p_in, h0, s0):

            rhotemp = rho_from_p_s(p_in, s0)
            htemp = h_from_p_rho(p_in, rhotemp)
            utemp = np.sqrt(2 * (h0 - htemp))
            atemp = a_from_p_s(p_in, s0)

            return utemp - atemp

        pstar = fsolve(solve_for_p, p0 / 1.8888, args=(h0, s0))
        rhostar = rho_from_p_s(pstar, s0)
        hstar = h_from_p_rho(pstar, rhostar)
        ustar = np.sqrt(2 * (h0 - hstar))
        astar = a_from_p_s(pstar, s0)
        Tstar = T_from_p_rho(pstar, rhostar)

        p_arr = []
        rho_arr = []
        u_arr = []
        m_arr = []
        T_arr = []

        n2 = []
        o = []
        o2 = []
        n = []
        no = []

        p_ideal = []
        rho_ideal = []
        T_ideal = []

    #    air = initialize_gas_object('air')

        # Iterate through area ratios
        for A_Astar in A_Astars:

            # Solve nozzle conditions at this area ratio
            def solve_for_p(p1, h0, s0, A_Astar):

                if p1 <= 0:
                    return 1000 + abs(p1)

                rhotemp = rho_from_p_s(p1, s0)
                htemp = h_from_p_rho(p1, rhotemp)

                utemp = np.sqrt(2 * (h0 - htemp))

                rhocomp = rhostar * ustar / (A_Astar * utemp)

                return rhotemp - rhocomp

            p1 = fsolve(solve_for_p, p0 / 20, args=(h0, s0, A_Astar))
            rho1 = rho_from_p_s(p1, s0)
            h1 = h_from_p_rho(p1, rho1)
            u1 = np.sqrt(2 * (h0 - h1))
            a1 = a_from_p_s(p1, s0)
            T1 = T_from_p_rho(p1, rho1)

            p_arr.append(p1/p0)
            rho_arr.append(rho1/rho0)
            u_arr.append(u1)
            m_arr.append(u1/a1)
            T_arr.append(T1/T0)

            gamma_air = 1.4
            p_ideal.append(1/isentropic.p0_p(u1/a1, gamma_air))
            rho_ideal.append(1/isentropic.r0_r(u1/a1, gamma_air))
            T_ideal.append(1/isentropic.T0_T(u1/a1, gamma_air))


#        T_converted = np.array(T_arr) * k_to_rankine * T0
#        ax6.semilogy(m_arr, T_converted, label=r'h0 = {} $\frac{{m^2}}{{s^2}}$'.format(h0))
#
#        ax7.plot(m_arr, np.array(u_arr) / ft_to_m, label=r'h0 = {} $\frac{{m^2}}{{s^2}}$'.format(h0))

#        air.TP = T1, p1
#        air.equilibrate('TP')

#        cdict = calculate_mole_mass_ratios(air)
#        n2.append(cdict['N2'])
#        o.append(cdict['O'])
#        o2.append(cdict['O2'])
#        n.append(cdict['N'])
#        no.append(cdict['NO'])

    # Plot p1, rho1, u1 and T1
    fig1, ax1 = plt.subplots(figsize=(12, 9))
    fig2, ax2 = plt.subplots(figsize=(12, 9))
    fig3, ax3 = plt.subplots(figsize=(12, 9))
    fig4, ax4 = plt.subplots(figsize=(12, 9))
#    fig5, ax5 = plt.subplots(figsize=(12, 9))

    ax1.plot(A_Astars, p_arr, label='Real Gas', lw=2.5, color="#0768C6")
    ax1.plot(A_Astars, p_ideal, label='Ideal Gas', lw=2.5, color="#ED4728", ls='--')

    ax2.plot(A_Astars, rho_arr, label='Real Gas', lw=2.5, color="#0768C6")
    ax2.plot(A_Astars, rho_ideal, label='Ideal Gas', lw=2.5, color="#ED4728", ls='--')

    ax3.plot(A_Astars, u_arr, label='Real Gas', lw=2.5, color="#0768C6")

    ax4.plot(A_Astars, T_arr, label='Real Gas', lw=2.5, color="#0768C6")
    ax4.plot(A_Astars, T_ideal, label='Ideal Gas', lw=2.5, color="#ED4728", ls='--')

#    ax5.plot(A_Astars, n2, label='N2')
#    ax5.plot(A_Astars, o, label='O')
#    ax5.plot(A_Astars, o2, label='O2')
#    ax5.plot(A_Astars, n, label='N')
#    ax5.plot(A_Astars, no, label='NO')
#
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    ax3.legend(loc='best')
    ax4.legend(loc='best')
    # ax5.legend(loc='best')

    # Plot formatting
    # ax1.set_xlabel(r'$\frac{A}{A^*}$', fontsize=18)
    ax1.set_xlabel(r'$Mach Number$', fontsize=18)
    # ax2.set_xlabel(r'$\frac{A}{A^*}$', fontsize=18)
    ax2.set_xlabel(r'$Mach Number$', fontsize=18)
    ax3.set_xlabel(r'$\frac{A}{A^*}$', fontsize=18)
    # ax4.set_xlabel(r'$\frac{A}{A^*}$', fontsize=18)
    ax4.set_xlabel(r'$Mach Number$', fontsize=18)
#    ax5.set_xlabel(r'$\frac{A}{A^*}$', fontsize=18)

    ax1.set_ylabel(r'$\frac{p}{p_0}$', fontsize=18)
    ax2.set_ylabel(r'$\frac{\rho}{\rho_0}$', fontsize=18)
    ax3.set_ylabel(r'$u (m/s)$', fontsize=18)
    ax4.set_ylabel(r'$\frac{T}{T_0}$', fontsize=18)
#    ax5.set_ylabel(r'$\eta$', fontsize=18)

    ax1.set_title('Variation of Nozzle Pressure Ratio for T0: {} K & p0: {} atm'.format(round(T0, 2), p0/101325), fontsize=16)
    ax2.set_title('Variation of Nozzle Density Ratio for T0: {} K & p0: {} atm'.format(round(T0, 2), p0/101325), fontsize=16)
    ax3.set_title('Variation of Nozzle Velocity for T0: {} K & p0: {} atm'.format(round(T0, 2), p0/101325), fontsize=16)
    ax4.set_title('Variation of Nozzle Temperature Ratio for T0: {} K & p0: {} atm'.format(round(T0, 2), p0/101325), fontsize=16)
#    ax5.set_title('Variation of Nozzle Chemical Composition', fontsize=16)

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()
#    fig5.tight_layout()

    fig1.savefig('p_profile.png', bbox_inches='tight')
    fig2.savefig('rho_profile.png', bbox_inches='tight')
    fig3.savefig('u_profile.png', bbox_inches='tight')
    fig4.savefig('T_profile.png', bbox_inches='tight')

    """COMPARISION PLOTS"""


#    import csv
#
#    with open('topic_3a/reference_T_M.csv','r') as csvfile:
#        exacts = csv.reader(csvfile, delimiter=',')
#        for row in exacts:
#            ax6.plot(float(row[0]), float(row[1]), marker='x', color='r', markersize=5)
#
#    with open('topic_3a/reference_u_M.csv','r') as csvfile:
#        exacts = csv.reader(csvfile, delimiter=',')
#        for row in exacts:
#            ax7.plot(float(row[0]), float(row[1]), marker='x', color='r', markersize=5)
#
#    ax6.set_xlim(1, 3.5)
#    ax7.set_xlim(1, 3.5)
#
#    ax6.legend(loc='best', fontsize='x-small')
#    ax7.legend(loc='best', fontsize='x-small')
#
#    ax6.set_xlabel(r'$M$', fontsize=18)
#    ax7.set_xlabel(r'$M$', fontsize=18)
#
#    ax6.set_ylabel(r'$T (\degree R)$', fontsize=18)
#    ax7.set_ylabel(r'$u (\frac{ft}{s})$', fontsize=18)
#
#    ax6.set_title('Variation of Nozzle Temperature with Mach # for p0 = {} atm'.format(p0/101325), fontsize=16)
#    ax7.set_title('Variation of Nozzle Velocity with Mach # for p0 = {} atm'.format(p0/101325), fontsize=16)
#
#    fig6.tight_layout()
#    fig7.tight_layout()
#
#    fig6.savefig('T_vs_M.png', bbox_inches='tight')
#    fig7.savefig('u_vs_M.png', bbox_inches='tight')
