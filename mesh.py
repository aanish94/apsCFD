#!/usr/bin/env python

""" Defines MeshGrid and MeshCell classes to sub-divide a geometry
into a computional domain via grid generation.

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

from geometry import create_spacecraft_geometry


class MeshGrid(object):
    """
    """

    def __init__(self):
        """
        """
        pass

#        resolution_x
#        resolution_y
#        # Specify boundary cells
#        cells_top
#        cells_bottom
#        cells_left
#        cells_right
#
#        # Specify interior cells
#        cells_interior
#
#        construct_grid_from_geometry()
#        visualize()
#        visualize_colormap(u or T or p or etc..)


class MeshCell(object):
    """
    """

#    def __init__(self):
#        """
#        """
#        pass
#
#    x, y, ID
#
#    rho, u, T, p, etc...
#    rho_old, u_old, T_old, p_old, etc...
#
#    get_left_neighbor()
#    get_right_neighbor()
#    get_top_neighbor()
#    get_bottom_neighbor()


def create_structured_mesh_via_algebraic_method(geometry, num_cells=None, cell_size=None):
    """Construct a structured mesh for the given physical domain via the algebraic method.

    :param <Geometry> geometry: A Geometry object representing the physical domain
    :param <list> num_cells: The number of cells in the x and y dimensions
    :param <float> cell_size: The desired size of each square cell

    :return <np.array> x_computational: 2-D X matrix representing a computational domain
    :return <np.array> y_computational: 2-D Y matrix representing a computational domain
    :return <tuple> stats: Statistics about the generated mesh

    """

    # Extract the x vector and associated max/min x values
    x_vector = geometry.x_vector
    x_max = max(x_vector)
    x_min = min(x_vector)

    # Extract the top & bottom y vectors and determine the max/min y values
    y_vector_top = geometry.y_vector_top
    y_vector_bottom = geometry.y_vector_bottom
    y_max = max(y_vector_top)
    y_min = min(y_vector_bottom)

    # Determine the number of cells in the x, y dimensions for either input
    # type
    if not num_cells and cell_size:
        x_num = int((x_max - x_min) / cell_size)
        y_num = int((y_max - y_min) / (2 * cell_size))

        x_cell_size = y_cell_size = cell_size
    elif num_cells:
        x_num = int(num_cells[0])
        y_num = int(num_cells[1] / 2)
        # Calculate associated cell size
        x_cell_size = (x_max - x_min) / x_num
        y_cell_size = (y_max - y_min) / y_num
    else:
        raise ValueError("Invalid inputs!")

    # Alert user if the desired grid resolution is too coarse
    if y_num < 0.5 * len(y_vector_top) or x_num < 0.5 * len(x_vector):
        print("WARNING! Input grid size {} may be too small.".format(
            (x_num, y_num * 2)))

    # Instantiate the x and y computational domains - the y domain will be
    # doubled later
    x_computational = np.zeros((x_num, y_num * 2))
    y_computational = np.zeros((x_num, y_num))

    # Scale the x vector to the size of the computational domain (xi)
    x_vector_scaled = np.linspace(x_min, x_max, num=x_num)

    # Define an equally spaced grid in x
    for jdx in np.arange(0, y_num * 2):
        x_computational[:, jdx] = x_vector_scaled

    # Scale the y vectors to the size of the computational domain (eta)
    scale_function_top = interpolate.interp1d(x_vector, y_vector_top)
    scale_function_bottom = interpolate.interp1d(x_vector, y_vector_bottom)

    y_vector_top_scaled = scale_function_top(x_vector_scaled)
    y_vector_bottom_scaled = scale_function_bottom(x_vector_scaled)

    # Define equally an equally spaced grid in y
    for idx in np.arange(0, x_num):
        for jdx in np.arange(0, y_num):
            y_height = y_vector_top_scaled[idx] - y_vector_bottom_scaled[idx]
            y_computational[idx, jdx] = y_vector_bottom_scaled[
                idx] + (y_height * jdx) / (y_num - 1)

    # Construct complete y grid by utilizing the symmetric nature of the
    # geometry
    y_computational = np.concatenate(
        (np.fliplr(-1 * y_computational), y_computational), axis=1)

    # Basic error checking
    if x_computational.shape != y_computational.shape:
        raise ValueError("X {} & Y {} shapes do not match!".format(
            x_computational.shape, y_computational.shape))

    return x_computational, y_computational, (x_num, y_num, x_cell_size, y_cell_size)


def visualize_mesh(geometry, x, y):
    """Visualize a mesh created on a physical domain.

    :param <Geometry> geometry: A Geometry object representing the physical domain
    :param <np.array> x_computational: 2-D X matrix representing a computational domain
    :param <np.array> y_computational: 2-D Y matrix representing a computational domain

    """

    fig, ax = plt.subplots(figsize=(12, 9))

    plt.plot(x, y)

    plt.margins(0.1, 0.1)


if __name__ == "__main__":
    pass
    sc = create_spacecraft_geometry()
    # x1, y1 = algebraic_structured_mesh(sc)

    # plot_meshgrid(x1, y1)
    x, y, stats = create_structured_mesh_via_algebraic_method(sc, num_cells=[
                                                              100, 100])
    # x, y = create_structured_mesh_via_algebraic_method(sc, cell_size=1)

    visualize_mesh(sc, x, y)
