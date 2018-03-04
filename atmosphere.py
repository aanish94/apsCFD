#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Implements the mathematical representation of the 1976 Committee on Extension to the Standard Atmosphere (COESA) United States standard lower atmospheric values for absolute temperature, pressure, density, and speed of sound for the input geopotential altitude.

Below 32,000 meters (approximately 104,987 feet), the U.S. Standard Atmosphere is identical with the Standard Atmosphere of the International Civil Aviation Organization (ICAO).

Matlab version: https://www.mathworks.com/help/aeroblks/coesaatmospheremodel.html

Notes
    - Based on the U.S. 1976 Standard Atmosphere.
    https://ntrs.nasa.gov/search.jsp?R=19770009539
"""

import math
import numpy as np
from scipy import interpolate, constants

# Constants

# Earth Radius [m]
EARTH_RADIUS = 6356.7660e3

# Gas Constant [J/m/K]
R = 8314.4598

# Atmospheric Conditions [N/m/m]
P0 = 101325.0

# Standard Gravity [m/s/s]
G0 = 9.80665

# Molecular Weight of Air [g/mol]
MW_AIR = 28.9644
GAMMA_AIR = 1.4

# Standard Properties

# Height [m]
H_STD = [
    0, 11E3, 20E3, 32E3, 47E3, 51E3, 71E3, 84.85205E3
]

# Temperature Gradient [K/m]
T_GRAD = [
    -6.5E-3, 0, 1E-3, 2.8E-3, 0, -2.8E-3, -2E-3, 0,
]

# Temperature [K]
T_STD = [
    288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946
]

# Pressure [Pa]
P_STD = [
    101325, 22632.06397346291, 5474.8886696777745, 868.0186847552279, 110.90630555496608,
    66.93887311868738, 3.956420428040732, 0.3733835899762159
]


def geometric_to_geopotential(z):
    """Calculates geopotential altitude from geometric altitude.

    :param <float> z: Geometric altitude (m)
    
    :return <float> h: Geopotential altitude (m)
    """

    return z * EARTH_RADIUS / (z + EARTH_RADIUS)


def geopotential_to_geometric(h):
    """Calculates geometric altitude from geopotential altitude.

    :param <float> h: Geopotential altitude (m)

    :return <float> z: Geometric altitude (m)    
    """

    return h * EARTH_RADIUS / (EARTH_RADIUS - h)


class Atmosphere(object):
    """
    US Standard Atmosphere 1976 class to calculate temperature, pressure,
    density, etc... as functions of altitude above sea level.

    Reference: [1] NOAA, NASA, and USAF. "U.S. Standard Atmosphere, 1976"
               October 15, 1976. http://ntrs.nasa.gov/search.jsp?R=19770009539.
    """

    def __init__(self, Z, kind='geometric', dT=0):
        """
        :param <float> Z: elevation (m)
        :param <float> dT: temperature delta from standard conditions (K)
        """

        self.Z = Z

        # Height
        if kind == 'geopotential':
            self.H = Z
        elif kind == 'geometric':
            self.H = geometric_to_geopotential(Z)
        else:
            raise ValueError("Unsupported kind: {}".format(kind))

        # Determine atmospheric layer
        idx = self.get_idx_from_height(self.H)
        self.T_layer = T_STD[idx]
        self.T_increase = T_GRAD[idx]
        self.P_layer = P_STD[idx]
        self.H_layer = H_STD[idx]

        self.H_above_layer = self.H - self.H_layer
        self.T = self.T_layer + self.T_increase * self.H_above_layer

        # determine corresponding pressure
        if self.T_increase == 0:
            self.P = self.P_layer * \
                math.exp(-G0 * MW_AIR * (self.H_above_layer) / R / self.T_layer)
        else:
            self.P = self.P_layer * \
                (self.T_layer / self.T) ** (G0 * MW_AIR / R / self.T_increase)

        if dT:
            self.T += dT

        # determine remaining properties
        self.rho = self.density(self.T, self.P)
        self.v_sonic = self.sonic_velocity(self.T)
        self.mu = self.viscosity(self.T)
        self.k = self.thermal_conductivity(self.T)
        self.g = self.gravity(self.Z)

    @staticmethod
    def get_idx_from_height(H):
        """Determine the index of the layer a specified elevation is above.
        Levels are fixed and defined in H_STD.

        :param <float> H: height (m)

        :return <int> Atmospheric layer index
        """

        if H <= 0:
            return 0

        for idx, Hi in enumerate(H_STD):
            if Hi >= H:
                return idx - 1

        return len(H_STD) - 1

    @staticmethod
    def density(T, P):
        """Density of air as a function of T and P.

        :param <float> T: temperature (K)
        :param <float> P: pressure (Pa)

        :return <flooat> Density (kg / m^3)
        """

        return MW_AIR / R * P / T

    @staticmethod
    def gravity(Z):
        """Gravitational acceleration above Earth as a function of elevationonly.

        :param <float> Z: elevation (m)

        :return <float> Gravitational acceleration (m/s/s)
        """

        return G0 * (EARTH_RADIUS / (EARTH_RADIUS + Z)) ** 2

    @staticmethod
    def sonic_velocity(T):
        """Speed of sound as a function of T only.

        :param <float> T: temperature (K)

        :return <float> Speed of sound (m/s)
        
        """
        return ((GAMMA_AIR * R / MW_AIR) * T) ** 0.5

    @staticmethod
    def thermal_conductivity(T):
        """Thermal conductivity of air [W/m/K] as a function of T only

        :param <float> T: temperature (K)

        :return <float> Thermal conductivity (W/m/K)
        """

        k = 2.64638E-3 * T ** 1.5 / (T + 245.4 * 10 ** (-12. / T))
        return k

    @staticmethod
    def viscosity(T):
        """Viscosity of air [Pa-s] as a function of T only

        :param <float> T: temperature [K]

        :return <float> Viscosity (Pa-s)
        """

        mu = 1.458E-6 * T ** 1.5 / (T + 110.4)
        return mu
