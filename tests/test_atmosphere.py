#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test atmosphere model
"""

import numpy as np
import unittest

from atmosphere import Atmosphere


class TestAtmosphere(unittest.TestCase):
    """
    """

    def test_sea_level(self):
        """Test Sea Level properties
        """

        h = 0
        a1 = Atmosphere(h)

        self.assertAlmostEqual(a1.T, 273 + 15, places=0)
        self.assertAlmostEqual(a1.P, 101325, places=0)
        self.assertAlmostEqual(a1.rho, 1.225, places=4)


    def test_under_1000m(self):
        """Tests for altitude values under 1000 meters
        """

        z = np.array([50.0, 550.0, 850.0])

        expected_T = np.array([287.825, 284.575, 282.626])
        expected_p = np.array([100720.0, 94890.0, 91523.0])
        expected_rho = np.array([1.2191, 1.1616, 1.1281])

        for idx, z_current in enumerate(z):
            a = Atmosphere(z_current)

            self.assertAlmostEqual(a.T, expected_T[idx], places=3)
            self.assertAlmostEqual(a.P, expected_p[idx], places=-2)
            self.assertAlmostEqual(a.rho, expected_rho[idx], places=4)


    def test_under_11km(self):
        """Tests for altitude values between 1 and 11 km
        """

        z = np.array([500.0, 2500.0, 6500.0, 9000.0, 11000.0])

        expected_T = np.array([284.900, 271.906, 245.943, 229.733, 216.774])
        expected_p = np.array([95461.0, 74691.0, 44075.0, 30800.0, 22699.0])
        expected_rho = np.array([1.1673, 0.95695, 0.62431, 0.46706, 0.36480])

        for idx, z_current in enumerate(z):
            a = Atmosphere(z_current)

            self.assertAlmostEqual(a.T, expected_T[idx], places=3)
            self.assertAlmostEqual(a.P, expected_p[idx], places=-2)
            self.assertAlmostEqual(a.rho, expected_rho[idx], places=4)


    def test_under_35km(self):
        """Tests for altitude values between 11 and 35 km
        """

        z = np.array([15000.0, 25000.0, 35000.0])
        expected_T = np.array([216.65, 221.552, 236.513])
        expected_p = np.array([12111.0, 2549.2, 574.59])
        expected_rho = np.array([0.19476, 0.040084, 0.0084634])

        for idx, z_current in enumerate(z):
            a = Atmosphere(z_current)

            self.assertAlmostEqual(a.T, expected_T[idx], places=3)
            self.assertAlmostEqual(a.P, expected_p[idx], places=-2)
            self.assertAlmostEqual(a.rho, expected_rho[idx], places=4)


    def test_under_86km(self):
        """Tests for altitude values between 35 and 86 km
        """

        z = np.array([50000.0, 70000.0, 86000.0])
        expected_T = np.array([270.65, 219.585, 186.87])
        expected_p = np.array([79.779, 5.2209, 0.37338])
        expected_rho = np.array([0.0010269, 0.000082829, 0.000006958])

        for idx, z_current in enumerate(z):
            a = Atmosphere(z_current)

            self.assertAlmostEqual(a.T, expected_T[idx], places=0)
            self.assertAlmostEqual(a.P, expected_p[idx], places=-2)
            self.assertAlmostEqual(a.rho, expected_rho[idx], places=4)


if __name__ == "__main__":
    unittest.main()
