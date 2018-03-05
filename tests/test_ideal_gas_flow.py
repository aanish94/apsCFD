#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Test ideal gas flow calculations
"""

import unittest

from utils.utils import deg_to_rad
from ideal_gas_flow import (isentropic, normal_shock, rayleigh,
                            prandtl_meyer, oblique_shock)


class TestIdealGasFlow(unittest.TestCase):
    """
    """

    def test_isentropic(self):
        """Test Isentropic Flow Relations
        """

        M = 2.0
        gamma = 1.4

        self.assertAlmostEqual(isentropic.A_Astar(M, gamma), 1.687, places=2)
        self.assertAlmostEqual(isentropic.T0_T(M, gamma), 1.80, places=2)
        self.assertAlmostEqual(isentropic.p0_p(M, gamma), 7.824, places=2)
        self.assertAlmostEqual(isentropic.r0_r(M, gamma), 4.347, places=2)

        # Inverse calculation
        self.assertAlmostEqual(isentropic.mach_from_area_ratio(1.687, gamma)[1], 2.0, places=2)

    def test_normal_shock(self):
        """Test Normal Shock Relations
        """

        M = 2.0
        gamma = 1.4

        self.assertAlmostEqual(normal_shock.mach(M, gamma), 0.5774, places=2)
        self.assertAlmostEqual(normal_shock.T2_T1(M, gamma), 1.687, places=2)
        self.assertAlmostEqual(normal_shock.p2_p1(M, gamma), 4.50, places=2)
        self.assertAlmostEqual(normal_shock.rho2_rho1(M, gamma), 2.667, places=2)

    def test_rayleigh(self):
        """Test Rayleigh Flow Relations
        """

        M = 2.0
        gamma = 1.4

        self.assertAlmostEqual(rayleigh.T0_T0star(M, gamma), 0.7934, places=2)
        self.assertAlmostEqual(rayleigh.T_Tstar(M, gamma), 0.5289, places=2)
        self.assertAlmostEqual(rayleigh.p_pstar(M, gamma), 0.3636, places=2)
        self.assertAlmostEqual(rayleigh.rho_rhostar(M, gamma), 0.6875, places=2)

        # Inverse calculations
        self.assertAlmostEqual(rayleigh.mach(0.7934, gamma), 2, places=2)
        
    def test_prandtl_meyer(self):
        """Test Prandtl Meyer Relations
        """

        M = 2.0
        gamma = 1.4

        self.assertAlmostEqual(prandtl_meyer.nu(M, gamma), 0.4604, places=2)
        
        # Inverse calculations
        self.assertAlmostEqual(prandtl_meyer.mach_from_nu(1.1481, gamma), 4, places=2)
        self.assertAlmostEqual(prandtl_meyer.mach_from_nu(0.4604, gamma), M, places=2)

    def test_oblique_shock(self):
        """Test Oblique Shock Relations
        """

        M = 2.0
        gamma = 1.4

        # From charts, for M = 2
        beta = deg_to_rad(44.0)
        theta = deg_to_rad(14.0)

        # Test beta <-> theta map
        self.assertAlmostEqual(oblique_shock.beta(M, theta, gamma), beta, places=2)
        self.assertAlmostEqual(oblique_shock.theta(M, beta, gamma), theta, places=2)
        
        # Test conditions behind the shock
        self.assertAlmostEqual(oblique_shock.mach(M, beta, theta, gamma), 1.482, places=1)
        self.assertAlmostEqual(oblique_shock.T2_T1(M, beta, gamma), 1.249, places=2)
        self.assertAlmostEqual(oblique_shock.p2_p1(M, beta, gamma), 2.088, places=2)
        self.assertAlmostEqual(oblique_shock.rho2_rho1(M, beta, gamma), 1.673, places=2)
        self.assertAlmostEqual(oblique_shock.u2_u1(M, beta, gamma), 0.8304, places=2)


if __name__ == "__main__":
    unittest.main()
