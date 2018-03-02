# -*- coding: utf-8 -*-
"""
Created on Tue Mar 07 21:12:53 2017

@author: asikora
"""

import math
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np

"""
Method of Characteristics implementation to design a nozzle for shock free flowflow
"""

# Constants
MAX_RGB = 255
MIN_RGB = 0
THROAT_MACH_NUM = 1


def rgb(value, minimum=MIN_RGB, maximum=MAX_RGB):
    """Calculate an RGB color scaled between the min and max limits to use in colormaps

    :param <float> value: Scaled value
    :param <int> minimum: Minimum RGB value
    :param <int> maximum: Maximum RGB value

    :return <str> hex: RGB values scaled between min/max and converted to HEX code
    """
    # Convert to floats
    minimum, maximum = float(minimum), float(maximum)

    # Determine scale of value
    ratio = 2 * (value - minimum) / (maximum - minimum)

    # Calculate RGB accordingly
    b = int(max(MIN_RGB, MAX_RGB * (1 - ratio)))
    r = int(max(MIN_RGB, MAX_RGB * (ratio - 1)))
    g = MAX_RGB - b - r

    # Convert to HEX
    hex_color = '#%02x%02x%02x' % (r, g, b)

    return hex_color


def intercept(x1, y1, t1, x2, y2, t2):
    """Calculate the intercept coordinates of two lines defined by an angle and start coordinates.

    :param <int> x1: Point A x-coordinate
    :param <int> y1: Point A y-coordinate
    :param <int> theta1: Point A angle
    :param <int> x2: Point B x-coordinate
    :param <int> y2: Point B y-coordinate
    :param <int> theta2: Point B angle

    :return <int> x: Intersection x-coordinate
    :return <int> y: Intersection y-coordinate
    """

    # Lines with the same angle do not intersect
    if t2 == t1:
        x = y = None
    # Calculate intersection
    else:
        x = ((y1 - y2) * np.cos(t1) * np.cos(t2) + x2 * np.sin(t2) *
             np.cos(t1) - x1 * np.sin(t1) * np.cos(t2)) / np.sin(t2 - t1)
        y = ((x1 - x2) * np.sin(t1) * np.sin(t2) + y2 * np.cos(t2) *
             np.sin(t1) - y1 * np.cos(t1) * np.sin(t2)) / np.sin(t1 - t2)

    return np.real(x), np.real(y)


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
                fdm = np.sqrt(m[i]**2 - 1) / (1 + 0.5 *
                                              (gamma - 1) * m[i] ** 2) / m[i]
                mnew = m[i] - fm / fdm
            m[i] = mnew

    return m


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
    mach_numbersi = m_nu(pmi, gamma)
    mach_anglesi = np.arcsin(1 / mach_numbersi)

    # Determine characteristic family i.e C- or C+
    slope = math.atan2(yi[1] - yi[0], xi[1] - xi[0]) - thetai[0]
    if abs(slope - mach_anglesi[0]) < abs(slope + mach_anglesi[0]):
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
    mach_angles[0, :] = mach_anglesi

    # Iterate through each grid point
    for i in range(1, nwaves):
        for j in range(i, nwaves):
            # Check if point is on the wall
            if i == j:
                # Apply boundary condition
                theta[i, j] = bc
                # Calculate prandtl-meyer and mach angle
                pm[i, j] = pm[i - 1, j] + c * (theta[i, j] - theta[i - 1, j])
                mach_angles[i, j] = np.arcsin(1 / m_nu(pm[i, j], gamma))
                # Calculate coordinates
                cur_angle = (theta[i - 1, j] + c * mach_angles[i - 1,
                                                               j] + theta[i, j] + c * mach_angles[i, j]) / 2
                prev_angle = (theta[i - 1, j - 1] + theta[i, j]) / 2
                x[i, j], y[i, j] = intercept(
                    x[i - 1, j], y[i - 1, j], cur_angle, x[i - 1, j - 1], y[i - 1, j - 1], prev_angle)

            # Normal grid point
            else:
                # Determine theta from previous point
                theta[i, j] = 0.5 * (theta[i, j - 1] + theta[i - 1,
                                                             j] + c * (pm[i, j - 1] - pm[i - 1, j]))
                pm[i, j] = pm[i, j - 1] + c * (theta[i, j - 1] - theta[i, j])
                mach_angles[i, j] = np.arcsin(1 / m_nu(pm[i, j], gamma))
                # Calculate coordinates
                cur_angle = (theta[i, j - 1] - c * mach_angles[i,
                                                               j - 1] + theta[i, j] - c * mach_angles[i, j]) / 2
                prev_angle = (theta[i - 1, j] + c * mach_angles[i - 1,
                                                                j] + theta[i, j] + c * mach_angles[i, j]) / 2
                x[i, j], y[i, j] = intercept(
                    x[i, j - 1], y[i, j - 1], cur_angle, x[i - 1, j], y[i - 1, j], prev_angle)

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
    mach_numbers = m_nu(pmi, gamma)
    mach_angles = np.arcsin(1 / mach_numbers)

    # Calculate coordinates of 1st centerline point
    x[1, 0] = xi[0] + leading_length * \
        np.cos(thetai[0] + characteristic_family * mach_angles[0])
    y[1, 0] = yi[0] + leading_length * \
        np.sin(thetai[0] + characteristic_family * mach_angles[0])

    x = np.real(x)
    y = np.real(y)

    # Calculate coordinates of all other points
    for j in range(1, nwaves):
        # Find the current and previous point angles
        cur_angle = thetai[j] + characteristic_family * mach_angles[j]
        prev_angle = (thetai[j - 1] + thetai[j] - characteristic_family *
                      (mach_angles[j - 1] + mach_angles[j])) / 2
        # Determine the intersection point
        x[1, j], y[1, j] = intercept(xi[j], yi[j], cur_angle, x[
                                     1, j - 1], y[1, j - 1], prev_angle)

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
    mach_numbers = m_nu(pmi, gamma)
    mach_angles = np.arcsin(1 / mach_numbers)

    # Calculate coordinates of 1st point
    x[1, 0], y[1, 0] = intercept(
        xi[0], yi[0], thetai[0] + c * mach_angles[0], x0, y0, (a0 + thetai[0]) / 2)

    # Calculate coordinates of all other points
    for j in range(1, nwaves):
        # Find the current and previous point angles
        cur_angle = thetai[j] + c * mach_angles[j]
        prev_angle = (thetai[j - 1] + thetai[j]) / 2
        # Determine the intersection point
        x[1, j], y[1, j] = intercept(xi[j], yi[j], cur_angle, x[
                                     1, j - 1], y[1, j - 1], prev_angle)

    return theta, pm, x, y


def calculate_color_scale(values, scale_min, scale_max):
    """Determine color based on input scale

    :param <list> values: List of values to scale
    :param <float> scale_min: Minimum scale value
    :param <float> scale_max: Maximum scale value

    :return <str> hex_color: HEX color code
    """

    # Scale number
    color = np.round(1 + MAX_RGB * (np.mean(values) -
                                    scale_min) / (scale_max - scale_min))

    # Handle too low/high
    if color < MIN_RGB:
        color = MIN_RGB
    elif color > MAX_RGB:
        color = MAX_RGB

    # Convert to HEX
    hex_color = rgb(color)

    return hex_color


def plot_polygon(ax, x, y, facecolor, edgecolor='black'):
    """Plot a polygon given array of x & coordinates

    :param <matplotlib.AxesSubplot> ax: Axis handle
    :param <np.array> x: X-coordinates
    :param <np.array> y: Y-coordinates
    :param <np.array> facecolor: Polygon face color
    :param <np.array> edgecolor: Polygon border color
    """

    poly = Polygon(np.c_[x, y], closed=True,
                   edgecolor=edgecolor, facecolor=facecolor)
    ax.add_patch(poly)

    return


def plot_complex_region(ax, theta, pm, x, y, gamma, cl, ch):
    """Plot Mach # contours for a complex wave region

    :param <matplotlib.AxesSubplot> ax: Axis handle
    :param <np.array> theta: Flow angles measured from the x-axis (radians)
    :param <np.array> pm: Prandtl-Meyer angles (radians)
    :param <np.array> x: X-coordinates along the outgoing wave boundary
    :param <np.array> y: Y-coordinates along the outgoing wave boundary
    :param <float> gamma: Specific heat ratio
    :param <float> cl: Minimum Mach #
    :param <float> ch: Maximum Mach #
    """

    # Count the number of waves
    nwaves = theta.shape[1]

    # Calculate Mach # from PM angles
    m = np.zeros(pm.shape)
    for i in range(len(pm)):
        m[i, :] = m_nu(pm[i, :], gamma)

    # Iterate through each grid point
    for i in range(1, nwaves):
        for j in range(i, nwaves):
            # Check if point is on the wall
            if i == j:
                xp = (x[i, j], x[i - 1, j - 1], x[i - 1, j])
                yp = (y[i, j], y[i - 1, j - 1], y[i - 1, j])

                # Use Mach # to determine color scale
                color = calculate_color_scale(
                    [m[i, j], m[i - 1, j - 1], m[i - 1, j]], cl, ch)
            # Normal grid point
            else:
                xp = (x[i, j], x[i, j - 1], x[i - 1, j - 1], x[i - 1, j])
                yp = (y[i, j], y[i, j - 1], y[i - 1, j - 1], y[i - 1, j])
                # Use Mach # to determine color scale
                color = calculate_color_scale(
                    [m[i, j], m[i, j - 1], m[i - 1, j - 1], m[i - 1, j]], cl, ch)

            plot_polygon(ax, xp, yp, color)

    return


def plot_simple_region(ax, theta, pm, x, y, gamma, cl, ch):
    """Plot Mach # contours for a simple wave region

    :param <matplotlib.AxesSubplot> ax: Axis handle
    :param <np.array> theta: Flow angles measured from the x-axis (radians)
    :param <np.array> pm: Prandtl-Meyer angles (radians)
    :param <np.array> x: X-coordinates along the outgoing wave boundary
    :param <np.array> y: Y-coordinates along the outgoing wave boundary
    :param <float> gamma: Specific heat ratio
    :param <float> cl: Minimum Mach #
    :param <float> ch: Maximum Mach #
    """

    # Count the number of waves
    nwaves = theta.shape[1]

    # Calculate Mach # from PM angles
    m = np.zeros(pm.shape)
    for i in range(len(pm)):
        m[i, :] = m_nu(pm[i, :], gamma)

    # Iterate through each wave
    for j in range(1, nwaves):
        # Define x & y coordinates of wave region
        xp = (x[1, j], x[1, j - 1], x[0, j - 1], x[0, j])
        yp = (y[1, j], y[1, j - 1], y[0, j - 1], y[0, j])

        # Use Mach # to determine color scale
        color = calculate_color_scale(
            [m[1, j], m[1, j - 1], m[0, j - 1], m[0, j]], cl, ch)

        # Plot polygon
        plot_polygon(ax, xp, yp, color)

    return


def plot_uniform_region(ax, theta, pm, x, y, gamma, cl, ch):
    """Plot Mach # contours for a uniform wave region

    :param <matplotlib.AxesSubplot> ax: Axis handle
    :param <np.array> theta: Flow angles measured from the x-axis (radians)
    :param <np.array> pm: Prandtl-Meyer angles (radians)
    :param <np.array> x: X-coordinates along the outgoing wave boundary
    :param <np.array> y: Y-coordinates along the outgoing wave boundary
    :param <float> gamma: Specific heat ratio
    :param <float> cl: Minimum Mach #
    :param <float> ch: Maximum Mach #
    """

    # Calculate Mach # from flow angles
    m = np.zeros(pm.shape)
    for i in range(len(pm)):
        m[i, :] = m_nu(pm[i, :], gamma)

    # Determine color and plot polygon
    color = calculate_color_scale(m[:], cl, ch)
    plot_polygon(ax, x, y, color)

    return


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
    astart = nu(mexit, gamma) / \
        2  # expansion angle of wall downstream of throat
    ai_1 = [astart / nwaves / 1000, astart / nwaves /
            100, astart / nwaves / 10, astart / nwaves / 2]
    ai_2 = list(np.linspace(astart / (nwaves - 1), astart, nwaves))

    ai = np.array(ai_1 + ai_2)
    xi = np.zeros(ai.shape)
    yi = np.ones(ai.shape) * throat_height
    ni = ai.copy()
    le = -throat_height
    mmax = m_nu(max(ai) * 2, gamma)

    cl = 1  # throat mach number
    ch = mmax

    # Simple Wave Region
    a1, n1, x1, y1 = simple_wave(ai, ni, xi, yi, le, gamma)
    # Complex Wave Reflection
    a2, n2, x2, y2 = complex_wave(
        a1[-1, :], n1[-1, :], x1[-1, :], y1[-1, :], 0, gamma)
    # simple wave boundary
    a3, n3, x3, y3 = simple_wave_boundary(
        a2[:, -1], n2[:, -1], x2[:, -1], y2[:, -1], 1, ai[-1], xi[-1], yi[-1],  gamma)

    # plot solution
    if plot_fig:
        fig, ax = plt.subplots(figsize=(12, 9))

        plot_simple_region(ax, a1, n1, x1, y1, gamma,
                           cl, ch)  # Initial wave front
        plot_complex_region(ax, a2, n2, x2, y2, gamma, cl,
                            ch)  # Reflection from wall
        plot_simple_region(ax, a3, n3, x3, y3, gamma, cl, ch)

        # Region between 1st simple region and complex region
        plot_uniform_region(ax, np.array([[a1[0, -1], a1[1, -1], a3[1, 0]]]), np.array([[n1[0, -1], n1[1, -1], n3[1, 0]]]), np.array(
            [x1[0, -1], x1[1, -1], x3[1, 0]]), np.array([y1[0, -1],  y1[1, -1], y3[1, 0]]), gamma, cl, ch)
        # plot_uniform_region(ax, np.array([[a1[0, 0], a1[1, 0], 0]]), np.array([[n1[0, 0], n1[1, 0], 0]]), np.array([x1[0, 0], x1[1, 0], x1[0, 0]]), np.array([y1[0, 0],  y1[1, 0], y1[1, 0]]), gamma, cl, ch)
        plot_uniform_region(ax, np.array([[a3[0, -1], a3[1, -1], a3[0, -1]]]), np.array([[n3[0, -1], n3[1, -1], n3[0, -1]]]), np.array(
            [x3[0, -1], x3[1, -1], x3[1, -1]]), np.array([y3[0, -1],  y3[1, -1], y3[0, -1]]), gamma, cl, ch)

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
    mach_number_exit = m_nu(n3[-1, :], gamma)[-1]

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
