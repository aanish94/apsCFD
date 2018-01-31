#!/usr/bin/env python

""" Tests fluid properties.
"""

import unittest

from fluid import Fluid


def check_near_equality(f1, f2, resolution=2):
    """Check if two numbers are 'nearly' equal.

    :param <float> f1: Number 1
    :param <float> f2: Number 2
    :param <int> resolution: Number of decimal points

    :return <bool> flag: True/False
    """

    return round(f1, resolution) == round(f2, resolution)


class TestFluid(unittest.TestCase):
    """
    """

    def test_helium(self):
        """
        """

        helium = Fluid(['Helium'])

        self.assertEqual(helium.name, 'Helium')

    def test_mixture(self):
        """
        """

        air_mix = Fluid(['Nitrogen', 'Oxygen'], [0.79, 0.21])
        air = Fluid(['Air'])

        self.assertEqual(air_mix.name, 'Nitrogen[0.79]&Oxygen[0.21]')

        d1 = air_mix.density(300, 101325)
        d2 = air.density(300, 101325)

        self.assertTrue(check_near_equality(d1, d2, 1))


if __name__ == "__main__":
    unittest.main()
