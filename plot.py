#!/usr/bin/env python
# -*- coding: utf-8 -*-

from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np

from utils import calculate_color_scale
"""
Plot Utilities
"""


def plot_polygon(ax, x, y, facecolor, edgecolor='black'):
    """Plot a polygon given array of x & coordinates

    :param <matplotlib.AxesSubplot> ax: Axis handle
    :param <np.array> x: X-coordinates
    :param <np.array> y: Y-coordinates
    :param <np.array> facecolor: Polygon face color
    :param <np.array> edgecolor: Polygon border color
    """

    poly = Polygon(np.c_[x, y], closed=True, edgecolor=edgecolor, facecolor=facecolor)
    ax.add_patch(poly)

    return


def plot_complex_region(ax, theta, m, x, y, gamma, cl, ch):
    """Plot Mach # contours for a complex wave region

    :param <matplotlib.AxesSubplot> ax: Axis handle
    :param <np.array> theta: Flow angles measured from the x-axis (radians)
    :param <np.array> m: Mach #'s
    :param <np.array> x: X-coordinates along the outgoing wave boundary
    :param <np.array> y: Y-coordinates along the outgoing wave boundary
    :param <float> gamma: Specific heat ratio
    :param <float> cl: Minimum Mach #
    :param <float> ch: Maximum Mach #
    """

    # Count the number of waves
    nwaves = theta.shape[1]

#    # Calculate Mach # from PM angles
#    m = np.zeros(pm.shape)
#    for i in range(len(pm)):
#        m[i, :] = m_nu(pm[i, :], gamma)

    # Iterate through each grid point
    for i in range(1, nwaves):
        for j in range(i, nwaves):
            # Check if point is on the wall
            if i == j:
                xp = (x[i, j], x[i - 1, j - 1], x[i - 1, j])
                yp = (y[i, j], y[i - 1, j - 1], y[i - 1, j])

                # Use Mach # to determine color scale
                color = calculate_color_scale([m[i, j], m[i - 1, j - 1], m[i - 1, j]], cl, ch)
            # Normal grid point
            else:
                xp = (x[i, j], x[i, j - 1], x[i - 1, j - 1], x[i - 1, j])
                yp = (y[i, j], y[i, j - 1], y[i - 1, j - 1], y[i - 1, j])
                # Use Mach # to determine color scale
                color = calculate_color_scale([m[i, j], m[i, j - 1], m[i - 1, j - 1], m[i - 1, j]], cl, ch)

            plot_polygon(ax, xp, yp, color)

    return


def plot_simple_region(ax, theta, m, x, y, gamma, cl, ch):
    """Plot Mach # contours for a simple wave region

    :param <matplotlib.AxesSubplot> ax: Axis handle
    :param <np.array> theta: Flow angles measured from the x-axis (radians)
    :param <np.array> m: Mach #'s
    :param <np.array> x: X-coordinates along the outgoing wave boundary
    :param <np.array> y: Y-coordinates along the outgoing wave boundary
    :param <float> gamma: Specific heat ratio
    :param <float> cl: Minimum Mach #
    :param <float> ch: Maximum Mach #
    """

    # Count the number of waves
    nwaves = theta.shape[1]

#    # Calculate Mach # from PM angles
#    m = np.zeros(pm.shape)
#    for i in range(len(pm)):
#        m[i, :] = m_nu(pm[i, :], gamma)

    # Iterate through each wave
    for j in range(1, nwaves):
        # Define x & y coordinates of wave region
        xp = (x[1, j], x[1, j - 1], x[0, j - 1], x[0, j])
        yp = (y[1, j], y[1, j - 1], y[0, j - 1], y[0, j])

        # Use Mach # to determine color scale
        color = calculate_color_scale([m[1, j], m[1, j - 1], m[0, j - 1], m[0, j]], cl, ch)

        # Plot polygon
        plot_polygon(ax, xp, yp, color)

    return


def plot_uniform_region(ax, theta, m, x, y, gamma, cl, ch):
    """Plot Mach # contours for a uniform wave region

    :param <matplotlib.AxesSubplot> ax: Axis handle
    :param <np.array> theta: Flow angles measured from the x-axis (radians)
    :param <np.array> m: Mach #'s
    :param <np.array> x: X-coordinates along the outgoing wave boundary
    :param <np.array> y: Y-coordinates along the outgoing wave boundary
    :param <float> gamma: Specific heat ratio
    :param <float> cl: Minimum Mach #
    :param <float> ch: Maximum Mach #
    """

    # Calculate Mach # from flow angles
#    m = np.zeros(pm.shape)
#    for i in range(len(pm)):
#        m[i, :] = m_nu(pm[i, :], gamma)

    # Determine color and plot polygon
    color = calculate_color_scale(m[:], cl, ch)
    plot_polygon(ax, x, y, color)

    return
