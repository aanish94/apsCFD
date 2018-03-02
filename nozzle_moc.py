#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import matplotlib.pyplot as plt
import numpy as np

from plot import plot_complex_region, plot_simple_region, plot_uniform_region
from utils import intercept
"""
Method of Characteristics implementation to design a nozzle for shock free flowflow
"""

# Constants
THROAT_MACH_NUM = 1


def calculate_mach_angles(mach_numbers):
    """Calculate Mach wave angle. A Mach wave propagates across the flow at the
    Mach angle, mu, which is the angle formed between the Mach wave wavefront and
    a vector pointing opposite to thedirection of motion.

    :param <np.array> mach_numbers: Mach numbers

    :return <np.array> mu: Mach angles
    """
    
    mach_angles = np.arcsin(1 / mach_numbers)
    
    return mach_angles


def calculate_mach_numbers(pm_angles, gamma):
    """Compute the Mach # given by a set of Prandtl-Meyer angles

    :param <np.array> pm_angles: PM angles (radians)
    :param <float> gamma: Specific heat ratio

    :return <np.array> mach_numbers: Mach #'s
    """
    # Track input type - Convert to 2-D array if not already
    flag_1d = False

    # Input is a number
    if isinstance(pm_angles, (int, float)):
        pm_angles = np.matrix(pm_angles)
    else:
        # Input is a 1-D array
        if pm_angles.ndim == 1:
            flag_1d = True
            pm_angles = np.matrix(pm_angles).transpose()
    
    m = np.zeros(pm_angles.shape)
    numax = (np.sqrt((gamma + 1) / (gamma - 1)) - 1) * np.pi / 2

    for i in range(pm_angles.shape[0]):
        for j in range(pm_angles.shape[1]):
            if pm_angles[i, j] < 0 or pm_angles[i, j] > numax:
                m[i, j] = np.nan
            elif pm_angles[i, j] == 0:
                m[i, j] = 1
            else:
                mnew = 2
                m[i, j] = 0
                while (abs(mnew - m[i, j]) > 0.00001):
                    m[i, j] = mnew
                    fm = nu(m[i, j], gamma) - pm_angles[i, j]
                    fdm = np.sqrt(m[i, j]**2 - 1) / (1 + 0.5 * (gamma - 1) * m[i, j] ** 2) / m[i, j]
                    mnew = m[i, j] - fm / fdm
                m[i, j] = mnew

    if flag_1d:
        return np.squeeze(np.asarray(m))
    else:
        return m


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


def complex_wave(thetai, pmi, xi, yi, bc, gamma):
    """
    Implement MOC in a complex wave region formed by the reflection of a simple wave
    from a plane wall. MOC solves for the flow & wave geometry in this triangular complex wave
    region.

    Inputs are of size 1xN where N is the number of wavelets to be computed. They are ordered
    with the upstream wavelet first.

    Outputs are of size NxN where:
        - 1st row represents the incoming wave boundary (replicates the inputs)
        - Last column represents the outgoing wave boundary
        - Main diagonal represents values on the wall

    :param <np.array> thetai: Flow angles measured from the x-axis (radians)
    :param <np.array> pmi: Prandtl-Meyer angles (radians)
    :param <np.array> xi: X-coordinates along the incoming wave boundary
    :param <np.array> yi: Y-coordinates along the incoming wave boundary
    :param <float> bc: Angle of the wall from which the wave is reflecting (radians)
    :param <float> gamma: Specific heat ratio

    :return <np.array> theta: Flow angles measured from the x-axis (radians)
    :return <np.array> pm: Prandtl-Meyer angles (radians)
    :return <np.array> x: X-coordinates along the outgoing wave boundary
    :return <np.array> y: Y-coordinates along the outgoing wave boundary
    """
    # Calculate Mach # and angle
    mach_numbers_initial = calculate_mach_numbers(pmi, gamma)
    mach_angles_initial = calculate_mach_angles(mach_numbers_initial)

    # Determine characteristic family i.e C- or C+
    slope = math.atan2(yi[1] - yi[0], xi[1] - xi[0]) - thetai[0]
    if abs(slope - mach_angles_initial[0]) < abs(slope + mach_angles_initial[0]):
        # C- characteristic
        c = -1
    else:  # C+ characteristic
        c = 1

    # Initialize outputs
    nwaves = len(thetai)
    theta = np.zeros((nwaves, nwaves))
    pm = theta.copy()
    x = theta.copy()
    y = theta.copy()
    mach_angles = theta.copy()

    theta[0, :] = thetai
    pm[0, :] = pmi
    x[0, :] = xi
    y[0, :] = yi
    mach_angles[0, :] = mach_angles_initial

    # Iterate through each grid point
    for i in range(1, nwaves):
        for j in range(i, nwaves):
            # Check if point is on the wall
            if i == j:
                # Apply boundary condition
                theta[i, j] = bc
                # Calculate prandtl-meyer and mach angle
                pm[i, j] = pm[i - 1, j] + c * (theta[i, j] - theta[i - 1, j])
                mach_angles[i, j] = calculate_mach_angles(calculate_mach_numbers(pm[i, j], gamma))
                # Calculate coordinates
                cur_angle = (theta[i - 1, j] + c * mach_angles[i - 1,j] + theta[i, j] + c * mach_angles[i, j]) / 2
                prev_angle = (theta[i - 1, j - 1] + theta[i, j]) / 2
                x[i, j], y[i, j] = intercept(x[i - 1, j], y[i - 1, j], cur_angle, x[i - 1, j - 1], y[i - 1, j - 1], prev_angle)

            # Normal grid point
            else:
                # Determine theta from previous point
                theta[i, j] = 0.5 * (theta[i, j - 1] + theta[i - 1,j] + c * (pm[i, j - 1] - pm[i - 1, j]))
                pm[i, j] = pm[i, j - 1] + c * (theta[i, j - 1] - theta[i, j])
                mach_angles[i, j] = calculate_mach_angles(calculate_mach_numbers(pm[i, j], gamma))
                # Calculate coordinates
                cur_angle = (theta[i, j - 1] - c * mach_angles[i,j - 1] + theta[i, j] - c * mach_angles[i, j]) / 2
                prev_angle = (theta[i - 1, j] + c * mach_angles[i - 1,j] + theta[i, j] + c * mach_angles[i, j]) / 2
                x[i, j], y[i, j] = intercept(x[i, j - 1], y[i, j - 1], cur_angle, x[i - 1, j], y[i - 1, j], prev_angle)

    return theta, pm, x, y


def simple_wave(thetai, pmi, xi, yi, leading_length, gamma):
    """
    Implement MOC in a simple wave region where only straight waves of one characteristic
    family exist. MOC solves for the flow & wave geometry in simple wave region terminated by
    a characteristic of the opposite sign.

    Inputs are of size 1xN where N is the number of wavelets to be computed. They are ordered
    with the upstream wavelet first.

    Outputs are of size 2xN where the first row replicates the inputs and second row provides
    the values at the outgoing wave boundary.

    :param <np.array> thetai: Flow angles measured from the x-axis (radians)
    :param <np.array> pmi: Prandtl-Meyer angles (radians)
    :param <np.array> xi: X-coordinates along the incoming wave boundary
    :param <np.array> yi: Y-coordinates along the incoming wave boundary
    :param <float> leading_length: Leading wavelet length (sign gives characteristic family, C- or C+)
    :param <float> gamma: Specific heat ratio

    :return <np.array> theta: Flow angles measured from the x-axis (radians)
    :return <np.array> pm: Prandtl-Meyer angles (radians)
    :return <np.array> x: X-coordinates along the outgoing wave boundary
    :return <np.array> y: Y-coordinates along the outgoing wave boundary
    """

    # Determine characteristic family i.e C- or C+
    characteristic_family = np.sign(leading_length)  # 1 or -1
    leading_length = abs(leading_length)

    # Initialize outputs
    nwaves = len(thetai)
    theta = np.zeros((2, nwaves))

    pm = theta.copy()
    pm[:, :] = pmi

    x = theta.copy()
    x[0, :] = xi

    y = theta.copy()
    y[0, :] = yi

    theta[:, :] = thetai

    # Calculate Mach # and angle
    mach_numbers = calculate_mach_numbers(pmi, gamma)
    mach_angles = calculate_mach_angles(mach_numbers)

    # Calculate coordinates of 1st centerline point
    x[1, 0] = xi[0] + leading_length * np.cos(thetai[0] + characteristic_family * mach_angles[0])
    y[1, 0] = yi[0] + leading_length * np.sin(thetai[0] + characteristic_family * mach_angles[0])

    x = np.real(x)
    y = np.real(y)

    # Calculate coordinates of all other points
    for j in range(1, nwaves):
        # Find the current and previous point angles
        cur_angle = thetai[j] + characteristic_family * mach_angles[j]
        prev_angle = (thetai[j - 1] + thetai[j] - characteristic_family * (mach_angles[j - 1] + mach_angles[j])) / 2
        # Determine the intersection point
        x[1, j], y[1, j] = intercept(xi[j], yi[j], cur_angle, x[1, j - 1], y[1, j - 1], prev_angle)

    return theta, pm, x, y


def simple_wave_boundary(thetai, pmi, xi, yi, c, a0, x0, y0, gamma):
    """
    Implement MOC in a simple wave region along with the wall shape needed to cancel its
    reflection. MOC solves for the flow & wave geometry.

    Inputs are of size 1xN where N is the number of wavelets to be computed. They are ordered
    with the upstream wavelet first.

    Outputs are of size 2xN where the first row replicates the inputs and second row provides
    the values at the outgoing wave boundary.

    :param <np.array> thetai: Flow angles measured from the x-axis (radians)
    :param <np.array> pmi: Prandtl-Meyer angles (radians)
    :param <np.array> xi: X-coordinates along the incoming wave boundary
    :param <np.array> yi: Y-coordinates along the incoming wave boundary
    :param <int> c: Characteristic family (C+ or C-)
    :param <float> theta0: Flow angle of the last point on the flow boundary (radians)
    :param <float> x0: X-coordinate of the last point on the flow boundary
    :param <float> y0: Y-coordinate of the last point on the flow boundary
    :param <float> gamma: Specific heat ratio

    :return <np.array> theta: Flow angles measured from the x-axis (radians)
    :return <np.array> pm: Prandtl-Meyer angles (radians)
    :return <np.array> x: X-coordinates along the outgoing wave boundary
    :return <np.array> y: Y-coordinates along the outgoing wave boundary
    """

    # Initialize outputs
    nwaves = len(thetai)
    theta = np.zeros((2, nwaves))

    pm = theta.copy()
    pm[:, :] = pmi

    x = theta.copy()
    x[0, :] = xi

    y = theta.copy()
    y[0, :] = yi

    theta[:, :] = thetai

    # Calculate Mach # and angle
    mach_numbers = calculate_mach_numbers(pmi, gamma)
    mach_angles = calculate_mach_angles(mach_numbers)

    # Calculate coordinates of 1st point
    x[1, 0], y[1, 0] = intercept(xi[0], yi[0], thetai[0] + c * mach_angles[0], x0, y0, (a0 + thetai[0]) / 2)

    # Calculate coordinates of all other points
    for j in range(1, nwaves):
        # Find the current and previous point angles
        cur_angle = thetai[j] + c * mach_angles[j]
        prev_angle = (thetai[j - 1] + thetai[j]) / 2
        # Determine the intersection point
        x[1, j], y[1, j] = intercept(xi[j], yi[j], cur_angle, x[1, j - 1], y[1, j - 1], prev_angle)

    return theta, pm, x, y


def moc_minimum_length_nozzle(mexit, throat_height, gamma, nwaves, plot_fig):
    """
    compute the shape of and flow through a minimum length nozzle given a design
    exit mach number using method of charateristics (moc), unit process

    involves solving compatibility equations step by step along the characteristic
    lines (along which flow variables are continuous, but deritivatics are indeterminate)

    moc allows for designing supersonic nozzle contour for shock-free, istentropic flow
    that accounts for multi-dimensional nature. subsonic flow is accelerated to sonic 
    speed at the throat (sonic line usually curved due to 2-d nature, but is assumed
    straight). 3 distinct flow regimes form = simple wave, complex wave and simple wave
    with boundary

    rocket nozzles must be as short as possible to minimize weight - aka minimum length
    nozzle - wherein expansion section is shrunk to a point and expansion takes
    place through a centered prandtl-meyer wave emanating from a sharp corner throat 
    with angle theta_wmax. moc determines length L such that shock-free flow results and
    a nozzle shorter than L will have shocks within the nozzle.

    for a minimum length nozzle, the expansion angle of the wall downstream of the throat
    is equal to 1/2 the prandtl meyer function for the design exit mach number

    input:
        mexit - exit mach number (design parameter)
        throat_height - height of throat (m)
        gamma - specific heat ratio
        nwaves - # of increments for corner expansion fan (higher number - more accurate calculation)

    output:
        ... - nozzle countour
        nozzle_length - nozzle length
        nozzle_height - nozzle height at exit
    """
    astart = nu(mexit, gamma) / 2  # expansion angle of wall downstream of throat
    ai_1 = [astart / nwaves / 1000, astart / nwaves /
            100, astart / nwaves / 10, astart / nwaves / 2]
    ai_2 = list(np.linspace(astart / (nwaves - 1), astart, nwaves))

    ai = np.array(ai_1 + ai_2)
    xi = np.zeros(ai.shape)
    yi = np.ones(ai.shape) * throat_height
    ni = ai.copy()
    le = -throat_height

    cl = 1  # throat mach number
    ch = mexit

    # Simple Wave Region
    a1, n1, x1, y1 = simple_wave(ai, ni, xi, yi, le, gamma)
    # Complex Wave Reflection
    a2, n2, x2, y2 = complex_wave(a1[-1, :], n1[-1, :], x1[-1, :], y1[-1, :], 0, gamma)
    # simple wave boundary
    a3, n3, x3, y3 = simple_wave_boundary(a2[:, -1], n2[:, -1], x2[:, -1], y2[:, -1], 1, ai[-1], xi[-1], yi[-1],  gamma)

    # plot solution
    if plot_fig:
        fig, ax = plt.subplots(figsize=(12, 9))

        m1 = calculate_mach_numbers(n1, gamma)
        plot_simple_region(ax, a1, m1, x1, y1, gamma, cl, ch)  # Initial wave front
        
        m2 = calculate_mach_numbers(n2, gamma)                         
        plot_complex_region(ax, a2, m2, x2, y2, gamma, cl, ch)  # Reflection from wall
        
        m3 = calculate_mach_numbers(n3, gamma)                           
        plot_simple_region(ax, a3, m3, x3, y3, gamma, cl, ch)

        # Region between 1st simple region and complex region
        theta_uniform_1 = np.array([[a1[0, -1], a1[1, -1], a3[1, 0]]])
        pm_uniform_1 = np.array([[n1[0, -1], n1[1, -1], n3[1, 0]]])
        m_uniform_1 = calculate_mach_numbers(pm_uniform_1, gamma)
        x_uniform_1 = np.array([x1[0, -1], x1[1, -1], x3[1, 0]])
        y_uniform_1 = np.array([y1[0, -1],  y1[1, -1], y3[1, 0]])
        plot_uniform_region(ax, theta_uniform_1, m_uniform_1, x_uniform_1, y_uniform_1, gamma, cl, ch)

        # Region between complex region and exit
        theta_uniform_2 = np.array([[a3[0, -1], a3[1, -1], a3[0, -1]]])
        pm_uniform_2 = np.array([[n3[0, -1], n3[1, -1], n3[0, -1]]])
        m_uniform_2 = calculate_mach_numbers(pm_uniform_2, gamma)
        x_uniform_2 = np.array([x3[0, -1], x3[1, -1], x3[1, -1]])
        y_uniform_2 = np.array([y3[0, -1],  y3[1, -1], y3[0, -1]])
        plot_uniform_region(ax, theta_uniform_2, m_uniform_2, x_uniform_2, y_uniform_2, gamma, cl, ch)

        # Format the figure
        plt.autoscale(axis='both', tight=False)
        plt.margins(0.01, 0.01)
        plt.xlabel('$x$ [m]', fontsize=14)
        plt.ylabel('$y$ [m]', fontsize=14)
        plt.title('$M_{{exit}}={}$'.format(mexit), fontsize=20)
        plt.subplots_adjust(top=0.925)
        # plt.savefig("moc_nozzle.png", bbox_to_inches=True)

    # Nozzle Design Summary
    nozzle_length = x3[-1, -1]
    nozzle_height = y3[-1, -1]
    mach_number_exit = calculate_mach_numbers(n3[-1, :], gamma)[-1]

    print("nozzle length: {} [m]".format(nozzle_length))
    print("nozzle height: {} [m]".format(nozzle_height / throat_height))

    # Confirm exit Mach # is achieved
    assert(round(mach_number_exit, 1) == round(mexit, 1))

    return nozzle_length, nozzle_height / throat_height, np.insert(x3[-1, :], 0, x1[0, 0]), np.insert(y3[-1, :], 0, y1[0, 0])


if __name__ == "__main__":
    pass
    mach_number_exit = 2.4
    throat_height = 1.0
    gamma = 1.4
    number_waves = 20
    length, ae_as, x_vec, y_vec = moc_minimum_length_nozzle(mach_number_exit,
                                                            throat_height,
                                                            gamma,
                                                            number_waves,
                                                            True)
