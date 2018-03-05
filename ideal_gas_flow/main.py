#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

from isentropic import mach_angle
from oblique_shock import theta


def plot_shock_curve(ax, M, gamma, N=1000):
    """Plot the theta-beta curve corresponding to an oblique shock with
    given M
    """

    mu = mach_angle(M)
    beta_arr = np.linspace(mu, np.pi / 2, N)

    theta_arr = np.zeros_like(beta_arr)

    for i in range(N):
        theta_arr[i] = theta(M, beta_arr[i], gamma)

    ax.plot(np.degrees(theta_arr), np.degrees(beta_arr), label=r'$M = {}$'.
        format(M if np.isfinite(M) else "\infty"))

    return


def create_oblique_shock_graph(M_list, gamma):
    """Construct an oblique shock graph
    """

    fig, ax = plt.subplots(figsize=(12, 9))

    for M in M_list:
        plot_shock_curve(ax, M, gamma)

    ax.grid(True)
    ax.set_ylim(0, 90)
    ax.set_xlim(0, 50)
    ax.set_xticks(np.arange(11) * 5)
    ax.set_xlabel(r"Deflection Angle $\theta$ [deg]")
    ax.set_ylabel(r"Shock Wave Angle $\beta$ [deg]")
    ax.legend(loc=4)

    ax.set_title(r"Oblique Shock Graph for $\gamma$ = {}".format(gamma))
    fig.savefig("oblique_shock_graph.png", bbox_inches="tight", dpi=100)

    return


if __name__ == "__main__":
    pass

    M_list =[1.5, 2.0, 3.0, 5.0, np.inf]
    gamma = 1.4

    create_oblique_shock_graph(M_list, gamma)
