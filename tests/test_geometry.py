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

        sc = geometry.create_sc_geometry()

        # Confirm array lengths
        self.assertEqual(sc.x_vector, sc.y_vector_top)
        self.assertEqual(sc.x_vector, sc.y_vector_bottom)
        self.assertEqual(sc.x_vector, sc.y_vector_top_reverse)
        self.assertEqual(sc.x_vector, sc.y_vector_bottom_reverse)

        # Confirm array criteria
        self.assertGreater(min(sc.y_vector_top), max(sc.y_vector_bottom))
        self.assertGreaterEqual(min(sc.y_vector_bottom), max(sc.y_vector_bottom_reverse))
        self.assertLess(max(sc.y_vector_top), min(sc.y_vector_bottom))

    def test_pressure_vessel(self):
        """
        """

        pv = geometry.create_sc_geometry()

        # Confirm array lengths
        self.assertEqual(pv.x_vector, pv.y_vector_top)
        self.assertEqual(pv.x_vector, pv.y_vector_bottom)
        self.assertEqual(pv.x_vector, pv.y_vector_top_reverse)
        self.assertEqual(pv.x_vector, pv.y_vector_bottom_reverse)

        # Confirm array criteria
        self.assertGreater(min(pv.y_vector_top), max(pv.y_vector_bottom))
        self.assertGreaterEqual(min(pv.y_vector_bottom), max(pv.y_vector_bottom_reverse))
        self.assertLess(max(pv.y_vector_top), min(pv.y_vector_bottom))


if __name__ == "__main__":
    unittest.main()
