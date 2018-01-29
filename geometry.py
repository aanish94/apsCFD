#!/usr/bin/env python

""" Provides Geometry class to specify the geometry of a symmetrical shape.

More Information:
"""

import matplotlib.pyplot as plt
import numpy as np


class Geometry(object):
    """Define a geometric domain using upper and lower y arrays.

    Geometry must be symmetric about the y-axis.
    """

    def __init__(self, x_max, x_min=None, resolution=1000,
                 bounds_top=None, funcs_top=None,
                 bounds_bottom=None, funcs_bottom=None):
        """Instantiate Geometry class.

        :param <int> x_max: Maximum x-value [meters]
        :param <int> x_min: Minimum x-value [meters]
        :param <int> resolution: Size of the x-array
        :param <list> bounds_top: y upper domain bounds
        :param <list> funcs_top: y upper domain functions
        :param <list> bounds_bottom: y lower domain bounds
        :param <list> funcs_bottom: y lower domain functions
        """

        self.x_max = x_max
        self.x_min = x_min or 0
        resolution = resolution or 1000

        if self.x_max <= self.x_min:
            raise ValueError("x_max : {} must be larger than x_min : {}".
                             format(self.x_max, self.x_min))

        # Construct x vector
        self.x_vector = np.linspace(self.x_min, self.x_max, num=resolution)

        # Use input functions to determine upper and lower y vectors
        self.y_vector_top = self.construct_y_vector(bounds_top, funcs_top, 1)
        self.y_vector_bottom = self.construct_y_vector(
            bounds_bottom, funcs_bottom, 0)

        # Rotate y vectors about the y-axis to complete the geometry symmetry
        self.y_vector_top_reverse = -1 * self.y_vector_top
        self.y_vector_bottom_reverse = -1 * self.y_vector_bottom

    def construct_y_vector(self, bounds, funcs, default=0):
        """Create a piecewise function given domain bounds and functions.

        :param <list> bounds: y domain bounds
        :param <list> funcs: y domain functions
        :param <int> default: Default y-value

        :return <np.array> y_vector: Resultant array of y values
        """

        # Use default if no input is provided
        if not bounds and not funcs:
            return np.ones(self.x_vector.shape) * default

        if isinstance(funcs, (int, float)):
            return np.ones(self.x_vector.shape) * funcs

        if len(bounds) != len(funcs):
            raise ValueError("The number of bounds and functions must match!")

        # Create simple functions from the specified bounds
        conditional_list = []

        # Iterate through each bound
        for idx, bound in enumerate(bounds):

            # Handle first boundary
            if idx == 0:
                x_left = self.x_min
                x_right = bound
            # Handle last boundary
            elif idx == len(bounds) - 1:
                x_left = bound
                x_right = self.x_max
            # Handle middle boundaries
            else:
                x_left = bounds[idx - 1]
                x_right = bound

            # Construct conditional list using boundaries
            conditional_list.append(
                (self.x_vector >= x_left) & (self.x_vector < x_right))

        # Evaluate piece-wise defined function for a given x array
        y_vector = np.piecewise(self.x_vector, conditional_list, funcs)

        if min(y_vector) < 0:
            raise ValueError(
                "y-vector minimum : {} is below 0!".format(min(y_vector)))

        return y_vector

    def visualize(self, file_name=None):
        """Output a visualization of the geometry.

        :param <str> file_name: Save plot as...
        """

        # Instantiate figure and axes objects
        fig, ax = plt.subplots()

        # Plot the upper and bottom boundaries
        ax.plot(self.x_vector, self.y_vector_top,
                label='Top', color='r', lw=3)

        ax.plot(self.x_vector, self.y_vector_top_reverse,
                label='Top - Reverse', color='r', lw=3)

        # Plot vertical lines at the begining & end of the domain
        plt.vlines(x=self.x_min,
                   ymin=self.y_vector_top_reverse[0], ymax=self.y_vector_top[0],
                   color='r', lw=3)
        plt.vlines(x=self.x_max,
                   ymin=self.y_vector_top_reverse[-1], ymax=self.y_vector_top[-1],
                   color='r', lw=3)

        # Color in the space surrounding the geometry
        ax.fill_between(self.x_vector,
                        self.y_vector_top_reverse, self.y_vector_bottom_reverse,
                        facecolor='black')
        ax.fill_between(self.x_vector, self.y_vector_top, self.y_vector_bottom,
                        facecolor='black')

        # Turn off axes
        plt.axis('off')
        # Save figure to file
        file_name = file_name or 'geometry.png'
        plt.savefig(file_name, bbox_inches='tight', pad_inches=0)

        return


def create_spacecraft_geometry():
    """Create a Geometry class representing a spacecraft.

    :return <Geometry> spacecraft: Geometry class
    """

    bounds_lower = [3, 7, 33]
    funcs_lower = [0, lambda y: y ** 1.5, 0]

    bounds_upper = None
    funcs_upper = 100

    x_max = 10
    x_min = 0
    resolution = 10000

    spacecraft = Geometry(x_max, x_min, resolution,
                          bounds_upper, funcs_upper,
                          bounds_lower, funcs_lower)

    return spacecraft


def create_rocket_engine_geometry():
    """Create a Geometry class representing a rocket engine.
    """

    bounds_lower = [3, 7, 33]
    funcs_lower = [0, lambda y: y ** 1.5, 0]

    bounds_upper = None
    funcs_upper = 100

    x_max = 10
    x_min = 0
    resolution = 10000

    rocket_engine = Geometry(x_max, x_min, resolution,
                             bounds_upper, funcs_upper,
                             bounds_lower, funcs_lower)

    return rocket_engine


def create_pressure_vessel_geometry():
    """Create a Geometry class representing a pressure vessel.
    """

    bounds_lower = [3, 7, 33]
    funcs_lower = [0, lambda y: y ** 1.5, 0]

    bounds_upper = None
    funcs_upper = 100

    x_max = 10
    x_min = 0
    resolution = 10000

    pressure_vessel = Geometry(x_max, x_min, resolution,
                               bounds_upper, funcs_upper,
                               bounds_lower, funcs_lower)

    return pressure_vessel


if __name__ == "__main__":
    pass

#    x = np.linspace(-2.5, 2.5, 600)
#    y_top = np.piecewise(x, [x < 0, x >= 0], [lambda x: -x, lambda x: x**2])
#    y_bottom = np.zeros(x.shape)
#
#    plt.plot(x, y_top, 'r')
#    plt.plot(x, y_bottom, 'b')

    sc = create_spacecraft_geometry()
    sc.visualize()
