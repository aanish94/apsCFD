#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import matplotlib.pyplot as plt
import numpy as np

"""
General Utilities
"""

# Constants
MAX_RGB = 255
MIN_RGB = 0
THROAT_MACH_NUM = 1


def calculate_color(value, minimum=MIN_RGB, maximum=MAX_RGB):
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
    hex_color = calculate_color(color)

    return hex_color


def deg_to_rad(d):
    """Degrees to radians
    """
    return d / 180.0 * math.pi


def rad_to_deg(r):
    """Radians to degrees
    """
    return r * 180.0 / math.pi


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
