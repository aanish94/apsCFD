#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Equilibrium Normal Shock Relations
"""

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

from ideal_oblique_shock import beta, theta

# from atmosphere import Atmosphere
# from air_equilibrium_properties import a_from_p_s, h_from_p_rho, rho_from_p_s, T_from_p_rho
from normal_shock import normal_shock
from gas import initialize_gas_object

# Unit conversions
atm_to_pa = 101325
ft_to_m = 0.3048


def deg_to_rad(d):
    """Degrees to radians
    """
    return d / 180.0 * math.pi


def rad_to_deg(r):
    """Radians to degrees
    """
    return r * 180.0 / math.pi


if __name__ == "__main__":
    pass
    fig, ax = plt.subplots(figsize=(12, 9))

    air = initialize_gas_object('air')

    u_arr = np.array([5000, 10000, 15000, 20000, 25000]) * ft_to_m

    T1 = 226
    p1 = 1114.309329279853

    for u1 in u_arr:
        beta_arr = deg_to_rad(np.linspace(0, 90))

        theta_arr = []
        beta_arr_plot = []
        for beta_in in beta_arr:
            def solve_for_theta(t1):

                lhs = math.tan(beta_in - t1)
                vn1 = u1 * math.sin(beta_in)

                res = normal_shock(air, vn1, T1, p1)
                if not res:
                    return 1000 + abs(t1)
                vn2, _, _, _ = res

                rhs = math.tan(beta_in) * vn2 / vn1

                return lhs - rhs

            theta_guess = theta(10, beta_in, 1.4)

            try:
                theta_real = fsolve(solve_for_theta, x0=theta_guess)
            except ValueError:
                continue

            theta_arr.append(theta_real)
            beta_arr_plot.append(beta_in)
        # print rad_to_deg(theta_real)

        ax.plot(np.degrees(theta_arr), np.degrees(beta_arr_plot), label=r'$u_1 = {} \frac{{m}}{{s}}$'.format(u1), lw=2)

    import csv
    with open('topic_6/reference.csv','r') as csvfile:
        exacts = csv.reader(csvfile, delimiter=',')
        for row in exacts:
            ax.plot(float(row[0]), float(row[1]), marker='x', color='r', markersize=5)
        
    ax.set_xlabel(r'Flow deflection, $\theta \degree$', fontsize=18)
    ax.set_ylabel(r'Shock angle, $\beta \degree$', fontsize=18)
    ax.set_title('Deflection angle-wave angle-velocity diagram for oblique shocks at 100,000 ft', fontsize=16)
    ax.legend(loc='best')

    fig.tight_layout()
    fig.savefig('theta_beta.png', bbox_inches='tight')
