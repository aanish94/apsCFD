#!/usr/bin/env python

""" Tests geometry generation logic.
"""

import unittest

import geometry


class TestGeometry(unittest.TestCase):
    """
    """

    def test_spacecraft(self):
        """
        """

        sc = geometry.create_spacecraft_geometry()

        # Confirm array lengths
        self.assertEqual(len(sc.x_vector), len(sc.y_vector_top))
        self.assertEqual(len(sc.x_vector), len(sc.y_vector_bottom))
        self.assertEqual(len(sc.x_vector), len(sc.y_vector_top_reverse))
        self.assertEqual(len(sc.x_vector), len(sc.y_vector_bottom_reverse))

        # Confirm array criteria
        self.assertGreater(min(sc.y_vector_top), max(sc.y_vector_bottom))
        self.assertGreaterEqual(min(sc.y_vector_bottom),
                                max(sc.y_vector_bottom_reverse))
        self.assertLess(max(sc.y_vector_top_reverse),
                        min(sc.y_vector_bottom_reverse))

    def test_pressure_vessel(self):
        """
        """

        pv = geometry.create_pressure_vessel_geometry()

        # Confirm array lengths
        self.assertEqual(len(pv.x_vector), len(pv.y_vector_top))
        self.assertEqual(len(pv.x_vector), len(pv.y_vector_bottom))
        self.assertEqual(len(pv.x_vector), len(pv.y_vector_top_reverse))
        self.assertEqual(len(pv.x_vector), len(pv.y_vector_bottom_reverse))

        # Confirm array criteria
        self.assertGreater(min(pv.y_vector_top), max(pv.y_vector_bottom))
        self.assertGreaterEqual(min(pv.y_vector_bottom),
                                max(pv.y_vector_bottom_reverse))
        self.assertLess(max(pv.y_vector_top_reverse),
                        min(pv.y_vector_bottom_reverse))


if __name__ == "__main__":
    unittest.main()
